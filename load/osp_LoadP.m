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
reverseStr = '';
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
end
for kk = 1:MRSCont.nDatasets
    msg = sprintf('\nLoading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    fprintf([reverseStr, '\n', msg]);
    if MRSCont.flags.isGUI        
        set(progressText,'String' ,sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets));
    end
    
    if ~(MRSCont.flags.didLoadData == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'raw') && (kk > length(MRSCont.raw))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
        % Read in the raw metabolite data. Since the GE P-file loader needs
        % to know the number of sub-spectra (e.g. from spectral editing), the
        % type of sequence needs to be differentiated here already.
        if MRSCont.flags.isUnEdited
            [raw, raw_ref]  = io_loadspec_GE(MRSCont.files{kk},1);
        elseif MRSCont.flags.isMEGA
            [raw, raw_ref]  = io_loadspec_GE(MRSCont.files{kk},2);
        elseif MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
            [raw, raw_ref]  = io_loadspec_GE(MRSCont.files{kk},4);
        end
        MRSCont.raw_uncomb{kk}      = raw;
        MRSCont.raw_ref_uncomb{kk}  = raw_ref;
    end
end
time = toc(refLoadTime);
if MRSCont.flags.isGUI        
    set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
    pause(1);
end
fprintf('... done.\n Elapsed time %f seconds\n',time);
% Set flag
MRSCont.flags.coilsCombined     = 0;
MRSCont.runtime.Load = time;
end