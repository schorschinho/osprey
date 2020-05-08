function [MRSCont] = osp_LoadTwix(MRSCont)
%% [MRSCont] = osp_LoadTwix(MRSCont)
%   This function reads raw (un-combined) data from Siemens TWIX files.
%
%   USAGE:
%       [MRSCont] = osp_LoadTwix(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-20)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-02-20: First version of the code.

% Close any remaining open figures
close all;

if MRSCont.flags.hasRef
    if length(MRSCont.files_ref) ~= MRSCont.nDatasets
        error('Number of specified reference files does not match number of specified metabolite files.');
    end
end
if MRSCont.flags.hasWater
    if length(MRSCont.files_w) ~= MRSCont.nDatasets
        error('Number of specified water files does not match number of specified metabolite files.');
    end
end

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
        % Read in the raw metabolite data
        raw                         = io_loadspec_twix(MRSCont.files{kk});
        raw                         = op_leftshift(raw,raw.pointsToLeftshift);
        MRSCont.raw_uncomb{kk}      = raw;
        % Read in the raw reference data. If a reference exists, perform the
        % coil combination based on the reference, and perform an eddy current
        % correction. If not, combine metabolite data based on its own signal.
        if MRSCont.flags.hasRef
            raw_ref                     = io_loadspec_twix(MRSCont.files_ref{kk});
            raw_ref                     = op_leftshift(raw_ref,raw_ref.pointsToLeftshift);
            MRSCont.raw_ref_uncomb{kk}  = raw_ref;
        end
        if MRSCont.flags.hasWater
            raw_w                       = io_loadspec_twix(MRSCont.files_w{kk});
            raw_w                       = op_leftshift(raw_w,raw_w.pointsToLeftshift);
            MRSCont.raw_w_uncomb{kk}    = raw_w;
        end

        % Perform coil combination (SENSE-based reconstruction if PRIAM flag set)
        if ~MRSCont.flags.isPRIAM
                [MRSCont] = osp_combineCoils(MRSCont,kk);
        elseif MRSCont.flags.isPRIAM
            error('Coming soon!');
            %[MRSCont] = osp_senseRecon(MRSCont);
        end
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
end