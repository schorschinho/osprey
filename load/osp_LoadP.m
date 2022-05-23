function [MRSCont] = osp_LoadP(MRSCont)
%% [MRSCont] = osp_LoadP(MRSCont)
%   This function reads raw (coil-uncombined) data from GE P (*.7) files.
%
%   USAGE:
%       [MRSCont] = osp_LoadP(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-03-17)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-03-17: First version of the code.

% Close any remaining open figures
close all;


%% Get the data (loop over all datasets)
refLoadTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end
for kk = 1:MRSCont.nDatasets(1)
    for ll = 1: 1:MRSCont.nDatasets(2)
        [~] = printLog('OspreyLoad',kk,ll,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);

        if ~(MRSCont.flags.didLoad == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'raw') && (kk > length(MRSCont.raw))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
            % Read in the raw metabolite data. Since the GE P-file loader needs
            % to know the number of sub-spectra (e.g. from spectral editing), the
            % type of sequence needs to be differentiated here already.
            metab_ll = MRSCont.opts.MultipleSpectra.metab(ll);
            if MRSCont.flags.isUnEdited
                [raw, raw_ref]  = io_loadspec_GE(MRSCont.files{metab_ll,kk},1);
                raw.flags.UnEdited = 1;
                raw_ref.flags.UnEdited = 1;
            elseif MRSCont.flags.isMEGA
                [raw, raw_ref]  = io_loadspec_GE(MRSCont.files{metab_ll,kk},2);
                raw.flags.isMEGA = 1;
                raw_ref.flags.isMEGA = 1;
            elseif MRSCont.flags.isHERMES
                [raw, raw_ref]  = io_loadspec_GE(MRSCont.files{metab_ll,kk},4);
                raw.flags.isHERMES = 1;
                raw_ref.flags.isHERMES = 1;
            elseif MRSCont.flags.isHERCULES
                [raw, raw_ref]  = io_loadspec_GE(MRSCont.files{metab_ll,kk},4);
                raw.flags.isHERCULES = 1;
                raw_ref.flags.isHERCULES = 1;
            end
            MRSCont.raw_uncomb{metab_ll,kk}      = raw;
            MRSCont.raw_ref_uncomb{metab_ll,kk}  = raw_ref;
        end
    end
end
time = toc(refLoadTime);
[~] = printLog('done',time,ll,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);  
% Set flag
MRSCont.flags.coilsCombined     = 0;
MRSCont.runtime.Load = time;
end