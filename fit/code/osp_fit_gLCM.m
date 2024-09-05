function [MRSCont] = osp_fit_gLCM(MRSCont)
%% [MRSCont] = osp_fit_gLCM(MRSCont)
% If you are using osp_fit_gLCM please cite the following paper in addition
% to the original Osprey paper:
%   Zöllner HJ, Davies-Jenkins C, Simicic D, Tal A, Sulam J, Oeltzschner G. 
%   Simultaneous multi-transient linear-combination modeling of MRS data improves uncertainty estimation. 
%   Magn Reson Med. 2024 Sep;92(3):916-925. doi: 10.1002/mrm.30110
%
%
%   USAGE:
%       [MRSCont] = osp_fit_gLCM(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2024-08-22)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2024-08-22: First version of the code.


% Loop over all the datasets here
metFitTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end


% We want a loop over the extra dimension for separate fitting
SeparateExtraDims = 1;
if MRSCont.processed.metab{1}.dims.extras > 0
    SeparateExtraDims = MRSCont.processed.metab{1}.sz(MRSCont.processed.metab{1}.dims.extras);
end
ModelProcedure = jsonToStruct(MRSCont.opts.fit.ModelProcedure.metab{1,1});
if isstruct(ModelProcedure.Steps)
    ModelProcedureCell = cell(size(ModelProcedure.Steps));
    for ss = 1 : size(ModelProcedure.Steps,1)
        ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
    end
    ModelProcedure.Steps = ModelProcedureCell;
end
if ModelProcedure.Steps{1, 1}.extra.flag                                    % Model is actually 2D along extra dimension so we don't want to loop 
    SeparateExtraDims = 1;
end

for ex = 1 : SeparateExtraDims
    % We want a loop over all the model procedures which indicate the
    % subspectra to be modelled. This way it is possible to change the
    % subspectra in the model procedure. 
    for ss = 1 : size(MRSCont.opts.fit.ModelProcedure.metab,2)
        % Read model procedure 
        ModelProcedure = jsonToStruct(MRSCont.opts.fit.ModelProcedure.metab{1,ss});
        if isstruct(ModelProcedure.Steps)
            ModelProcedureCell = cell(size(ModelProcedure.Steps));
            for ss = 1 : size(ModelProcedure.Steps,1)
                ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
            end
            ModelProcedure.Steps = ModelProcedureCell;
        end
        if ~isfield(ModelProcedure,'basisset') || ~isfield(ModelProcedure.basisset, 'file') || ... 
            isempty(ModelProcedure.basisset.file)
            if ~iscell(MRSCont.fit.basisSet)
                ModelProcedure.basisset.file = {MRSCont.fit.basisSet};
            else
                ModelProcedure.basisset.file = MRSCont.fit.basisSet;
            end
        end
        if SeparateExtraDims > 1
            ModelProcedure.basisset.opts.index = ex;
        end
        if (ss == 1) && (ex == 1)
            [MRSCont.fit.results.metab(1,:,ss,ex)] = Osprey_gLCM(MRSCont.processed.metab,ModelProcedure)';
        else
            if isprop(MRSCont.fit.results.metab{1,1,1,1}, 'scale')
                scale = [];
                for kk = 1:MRSCont.nDatasets(1)
                    scale = [scale MRSCont.fit.results.metab{1,kk,1,1}.scale];
                end
            else
                scale = 0;
            end
            [MRSCont.fit.results.metab(1,:,ss,ex)] = Osprey_gLCM(MRSCont.processed.metab,ModelProcedure,0,0,scale)';
        end        
    end   
end
%% Model MM spectra
if MRSCont.flags.hasMM == 1
    ModelProcedure = jsonToStruct(MRSCont.opts.fit.ModelProcedure.mm{1});
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    if ~isfield(ModelProcedure,'basisset') || ~isfield(ModelProcedure.basisset, 'file') || ... 
        isempty(ModelProcedure.basisset.file)
        ModelProcedure.basisset.file = {MRSCont.fit.basisSet};
    end
    [MRSCont.fit.results.mm(1,:)] = Osprey_gLCM(MRSCont.processed.mm,ModelProcedure)';

    % Model conventional MRS with spline model of clean MM data
    for kk = 1:MRSCont.nDatasets(1)
        specs = MRSCont.fit.results.mm{1,kk}.Model{1, end}.fit.baseline;
        fids=ifft(fftshift(specs,1),[],1);
        t=MRSCont.fit.results.mm{1, kk}.Data.t';
        f = MRSCont.fit.results.mm{1,kk}.Model{1, end}.parsOut.freqShift(1);
        fids=fids.*exp(-1i*t*f*2*pi);
        specs=fftshift(fft(fids,[],1),1);
        MRSCont.processed.metab{kk}.MMExpSub.fids = fids;
        MRSCont.processed.metab{kk}.MMExpSub.specs = specs;
    end

    ModelProcedure = jsonToStruct(MRSCont.opts.fit.ModelProcedure.metab{2});
    if isstruct(ModelProcedure.Steps)
        ModelProcedureCell = cell(size(ModelProcedure.Steps));
        for ss = 1 : size(ModelProcedure.Steps,1)
            ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
        end
        ModelProcedure.Steps = ModelProcedureCell;
    end
    if ~isfield(ModelProcedure,'basisset') || ~isfield(ModelProcedure.basisset, 'file') || ... 
        isempty(ModelProcedure.basisset.file)
        ModelProcedure.basisset.file = {MRSCont.fit.basisSet};
    end
    [MRSCont.fit.results.metab(2,:)] = Osprey_gLCM(MRSCont.processed.metab,ModelProcedure)';

    % Run model with meanMM function
    if MRSCont.opts.fit.MeanMM
        for kk = 1:MRSCont.nDatasets(1)
            specs = MRSCont.fit.results.mm{1,kk}.Model{1, end}.fit.baseline;
            fids(:,kk)=ifft(fftshift(specs,1),[],1);
            t=MRSCont.fit.results.mm{1,kk}.Data.t';
            f = MRSCont.fit.results.mm{1,kk}.Model{1, end}.parsOut.freqShift(1);
            fids(:,kk)=fids(:,kk).*exp(-1i*t*f*2*pi);
        end
        mean_fids = mean(fids,2);
        mm_clean_spline.fids = mean_fids;
        for kk = 1:MRSCont.nDatasets(1)
            MRSCont.processed.metab{kk}.MMExpSub.fids = mean_fids;
            MRSCont.processed.metab{kk}.MMExpSub.specs = fftshift(fft(mean_fids,[],1),1);
        end
        [MRSCont.fit.results.metab(3,:)] = Osprey_gLCM(MRSCont.processed.metab,ModelProcedure)';
    end
end    

%% Timing
time = toc(metFitTime);
[~] = printLog('done',time,MRSCont.nDatasets,1,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
MRSCont.runtime.FitMet = time;

end
