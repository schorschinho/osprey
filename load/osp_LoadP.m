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
    progressbar = waitbar(0,'Start','Name','Osprey Load');
    waitbar(0,progressbar,sprintf('Loaded raw data from dataset %d out of %d total datasets...\n', 0, MRSCont.nDatasets))
end
for kk = 1:MRSCont.nDatasets
    msg = sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    if ((MRSCont.flags.didLoadData == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'raw') && (kk > length(MRSCont.raw))) || ~isfield(MRSCont.ver, 'Load') || ~strcmp(MRSCont.ver.Load,MRSCont.ver.CheckLoad))
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
    if MRSCont.flags.isGUI        
        waitbar(kk/MRSCont.nDatasets,progressbar,sprintf('Loaded raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets));
    end
end
fprintf('... done.\n');
if MRSCont.flags.isGUI 
    waitbar(1,progressbar,'...done')
    close(progressbar)
end
toc(refLoadTime);

% Set flag
MRSCont.flags.coilsCombined     = 0;

end