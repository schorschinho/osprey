function [MRSCont] = osp_fitUnEdited(MRSCont)
%% [MRSCont] = osp_fitUnEdited(MRSCont)
%   This function performs spectral fitting of unedited MRS data.
%
%   USAGE:
%       [MRSCont] = osp_fitUnEdited(MRSCont, basisSet);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       basisSet    = Osprey basis set container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-04-12)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-04-12: First version of the code.


% Loop over all the datasets here
metFitTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end

%% New Osprey gLCM model
if strcmpi(MRSCont.opts.fit.method, 'Osprey_gLCM')
    % Read model procedure 
    ModelProcedure = jsonToStruct(MRSCont.opts.fit.ModelProcedure.metab{1});
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
    [MRSCont.fit.results.metab] = Osprey_gLCM(MRSCont.processed.metab,ModelProcedure);
end

%% Old Osprey model or LCModel
if strcmpi(MRSCont.opts.fit.method, 'Osprey') || strcmpi(MRSCont.opts.fit.method, 'LCModel')
    if MRSCont.flags.hasMM == 1
        basisSet_mm = load(MRSCont.opts.fit.basisSetFile);
        basisSet_mm = basisSet_mm.BASIS;
        basisSet_mm = osp_recalculate_basis_specs(basisSet_mm);
        basisSet_mm = fit_sortBasisSet(basisSet_mm);
        metabList_mm = fit_createMetabListMM('unedited');
        basisSet_mm = fit_selectMetabs(basisSet_mm, metabList_mm, 1);
    end
    for kk = 1:MRSCont.nDatasets(1)   
        % ----- Osprey fit pipeline -----
        if strcmpi(MRSCont.opts.fit.method, 'Osprey')
            [~] = printLog('OspreyFit',kk,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
    
            if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > size(MRSCont.fit.results.metab.fitParams,2))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
                % Apply scaling factor to the data
                dataToFit   = MRSCont.processed.metab{kk};
                dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
                % Extract fit options
                fitOpts     = MRSCont.opts.fit;
                fitModel    = fitOpts.method;
    
                % If MRSI data load priors
                if MRSCont.flags.isMRSI
                    if isfield(MRSCont.fit, 'results')
                        fitParamsOff   = MRSCont.fit.results.metab.fitParams{kk};
                        %                 dataToFit.refShift   = fitParamsOff.refShift;
                        %                 dataToFit.refFWHM   = fitParamsOff.refFWHM;
                        fitOpts.MRSIpriors = fitParamsOff;
                    end
                end
    
                % Call the fit function
                basisSet    = MRSCont.fit.basisSet;
                fitOpts.GAP = MRSCont.opts.fit.GAP.A;
                [fitParams, resBasisSet] = fit_runFit(dataToFit, basisSet, fitModel, fitOpts);
    
                % Save back the basis set and fit parameters to MRSCont
                MRSCont.fit.basisSet                    = basisSet;
                MRSCont.fit.resBasisSet.metab{1,kk,1}         = resBasisSet;
                MRSCont.fit.results.metab.fitParams{1,kk,1}   = fitParams;
    
                 % Write NIfTI-MRS results
                % MRSCont.fit.nii_mrs.metab{kk} = osp_OspreyFitToNII(MRSCont.processed.metab{kk}, fitParams,resBasisSet,MRSCont.fit.scale{kk},fitOpts);
    
    
                %Modeling MM spectra after the main spectrum re_mm
                if MRSCont.flags.hasMM == 1              %re_mm
                    dataToFit_mm   = MRSCont.processed.mm{kk};
                    dataToFit_mm   = op_ampScale(dataToFit_mm, 1/MRSCont.fit.scale{kk});
                    %add some info from the metabolite fit
                    dataToFit_mm.lineShape  = fitParams.lineShape;
                    dataToFit_mm.refFWHM  = fitParams.refFWHM;
    %                 dataToFit_mm.refShift  = 0;
                    % Extract fit options
                    fitOpts_mm    = MRSCont.opts.fit;
                    fitOpts_mm.range = [0.2 4.2];
                    fitModel_mm    = fitOpts.method;
                    fitOpts_mm.sequence = 'unedited';
                    %Specify a reduced basis set for MM modeling
                    %basisSet_mm    = MRSCont.fit.basisSet;
                    %Reduce the size of the basis set
    
                    %Adjust basis set
                    % Clear existing basis set
                    MRSCont.fit.basisSet_mm = [];
                    % Load the specified basis set
                    basisSet_mm = load(MRSCont.opts.fit.basisSetFile);
                    basisSet_mm = basisSet_mm.BASIS;
                    % Generate the list of basis functions that are supposed to be included in
                    % the basis set
                    % To do: Interface with interactive user input
                    metabList_mm = fit_createMetabListMM('unedited');
                    % Collect MMfit flag from the options determined in the job file
                    fitMM = MRSCont.opts.fit.fitMM;
                    % Add basis spectra (if they were removed to reduce the file size)
                    if ~isfield(basisSet_mm,'specs')
                        [basisSet_mm] = osp_recalculate_basis_specs(basisSet_mm);
                    end
                    % Create the modified basis set
                    basisSet_mm = fit_selectMetabs(basisSet_mm, metabList_mm, fitMM);
    
                    % Call the fit function
                    [fitParams_mm, resBasisSet_mm] = fit_runFitMM(dataToFit_mm, basisSet_mm, fitModel_mm, fitOpts_mm);
                    % Save back the basis set and fit parameters to MRSCont
                    MRSCont.fit.basisSet_mm                    = basisSet_mm;
                    MRSCont.fit.resBasisSet.mm{1,kk,1}         = resBasisSet_mm;
                    MRSCont.fit.results.mm.fitParams{1,kk,1}   = fitParams_mm;
    
                     % Create full spectralwidth clean MM spectrum
                    disp('Create clean metabolite-nulled spectrum...');
                    [mm_clean] = fit_OspreyCleanMM(dataToFit_mm, resBasisSet_mm, fitOpts_mm,MRSCont.fit.scale{kk},fitParams_mm);
                    mm_clean.names = {'A_clean'};
                    mm_clean   = op_ampScale(mm_clean, MRSCont.fit.scale{kk});
    
    
                    % Create full spectralwidth clean MM spectrum spline model
                    disp('Create clean metabolite-nulled spline model...');
                    %Fit them with a spline?
                    fitOpts_mm_clean    = MRSCont.opts.fit;
                    fitOpts_mm_clean.range=[mm_clean.ppm(1) mm_clean.ppm(end)];
                    [~,mm_clean_spline] = fit_Osprey_SplineOnly(mm_clean, 0.1, fitOpts_mm_clean.range);
                    mm_clean_spline.names = {'A_spline'};
                    temp               = op_zeropad(mm_clean_spline, 16);
                    [refShift_mm, ~] = osp_XReferencing(temp,[0.91],[1],[0 4.2],1);% determine frequency shift
                    [mm_clean]             = op_freqshift(mm_clean,-refShift_mm);
                    [mm_clean_spline]             = op_freqshift(mm_clean_spline,-refShift_mm);
                    MRSCont.processed.mm{kk} = op_mergesubspec(MRSCont.processed.mm{kk},mm_clean);
                    MRSCont.processed.mm{kk} = op_mergesubspec(MRSCont.processed.mm{kk},mm_clean_spline);
    
    %                 % Gaussian fit with(out) baseline
    %                 disp('Fit clean metabolite-nulled spectrum...');
    %                 fitOpts_mm.Params2 = fitParams_mm;
    %                 dataToFit_mm_clean   = op_ampScale(dataToFit_mm_clean, 1/MRSCont.fit.scale{kk});
    %                 [fitParams_mm_Gauss, resBasisSet_mm_Gauss] = fit_runFitMM(dataToFit_mm_clean, basisSetMM,  'OspreyMMGaussian', fitOpts_mm);
    %                 MRSCont.fit.basisSet_mm_Gauss                   = basisSetMM;
    %                 MRSCont.fit.resBasisSet.diff1_mm_Gauss{kk}         = resBasisSet_mm_Gauss;
    %                 MRSCont.fit.results.diff1_mm.fitParams{kk}.fitParams_mm_Gauss   = fitParams_mm_Gauss;
    
    %                 % Fit diff1 with individual MM model
                    disp('Fit metabolite data with metabolite-nulled spline model...');
                    fitOpts.mm_clean_spline     = mm_clean_spline;
                    fitOpts.prelimParams = fitParams.prelimParams;
                    % Call the fit function
                    [fitParamsAExpMM, resBasisSetAExpMM]  = fit_runFit(dataToFit, basisSet, fitModel, fitOpts);
                    MRSCont.fit.resBasisSet.metab{2,kk,1}         = resBasisSetAExpMM;
                    MRSCont.fit.results.metab.fitParams{2,kk,1}   = fitParamsAExpMM;
                end
            end
    
            % ----- LCModel wrapper fit pipeline -----
        elseif strcmpi(MRSCont.opts.fit.method, 'LCModel')
             [MRSCont] = LCModelWrapper(MRSCont,kk,progressText);
        end
    
    end

    %Fit again with mean MM basis function
    if MRSCont.opts.fit.MeanMM
        if MRSCont.flags.hasMM == 1
    
            fids = zeros(mm_clean_spline.sz(1),MRSCont.nDatasets(1));
            for kk = 1:MRSCont.nDatasets(1)
                fids(:,kk)     = MRSCont.processed.mm{kk}.fids(:,3);
            end
            mean_fids = mean(fids,2);
            mm_clean_spline.fids = mean_fids;
            mm_clean_spline.specs = fftshift(fft(mean_fids,[],1),1);
            fitOpts.mm_clean_spline     = mm_clean_spline;
            MRSCont.opts.fit.mm_clean_spline= mm_clean_spline;
        end
    
        for kk = 1:MRSCont.nDatasets(1)
            [~] = printLog('OspreyFit',kk,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
             % Fit diff1 with individual MM model
            disp('Fit metabolite data with mean metabolite-nulled spline model...');
            % Call the fit function
            fitOpts.prelimParams = MRSCont.fit.results.metab.fitParams{1,kk}.prelimParams;
            fitOpts.GAP = MRSCont.opts.fit.GAP.A;
    
            dataToFit   = MRSCont.processed.metab{kk};
            dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
            dataToFit.refShift   = MRSCont.fit.results.metab.fitParams{1,kk,1}.refShift;
            dataToFit.refFWHM   = MRSCont.fit.results.metab.fitParams{1,kk,1}.refFWHM;
    
            [fitParamsAExpMM, resBasisSetAExpMM]  = fit_runFit(dataToFit, basisSet, fitModel, fitOpts);
            MRSCont.fit.resBasisSet.metab{3,kk,1}         = resBasisSetAExpMM;
            MRSCont.fit.results.metab.fitParams{3,kk,1}   = fitParamsAExpMM;
    
        end
    end
end


%% Timing
time = toc(metFitTime);
[~] = printLog('done',time,MRSCont.nDatasets,1,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
MRSCont.runtime.FitMet = time;

end
