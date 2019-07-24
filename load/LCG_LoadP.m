function [MRSCont] = LCG_LoadP(MRSCont)
%% [MRSCont] = LCG_LoadP(MRSCont)
%   This function reads raw (coil-uncombined) data from GE P (*.7) files.
%
%   USAGE:
%       [MRSCont] = LCG_LoadP(MRSCont);
%
%   INPUTS:
%       MRSCont     = LCGannet MRS data container.
%
%   OUTPUTS:
%       MRSCont     = LCGannet MRS data container.
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

% Determine number of datasets
MRSCont.nDatasets = length(MRSCont.files);

%% Get the data (loop over all datasets)
refLoadTime = tic;
reverseStr = '';
for kk = 1:MRSCont.nDatasets
    msg = sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
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
fprintf('... done.\n');
toc(refLoadTime);

% Set flag
MRSCont.flags.coilsCombined     = 0;

end