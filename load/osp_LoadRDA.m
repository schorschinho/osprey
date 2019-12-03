function [MRSCont] = osp_LoadRDA(MRSCont)
%% [MRSCont] = osp_LoadRDA(MRSCont)
%   This function reads raw (but coil-combined) data from RDA files
%   produced by Siemens sequences (*.rda).
%
%   USAGE:
%       [MRSCont] = osp_LoadRDA(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-06-19)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-06-19: First version of the code.

% Close any remaining open figures
close all;
warning('off','all');

% Determine number of datasets
MRSCont.nDatasets = length(MRSCont.files);
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
for kk = 1:MRSCont.nDatasets
    msg = sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    % Read in the raw metabolite data.
    raw = io_loadspec_rda(MRSCont.files{kk});
    MRSCont.raw{kk}      = raw;
    
    % Read in the raw reference data.
    if MRSCont.flags.hasRef
        raw_ref = io_loadspec_rda(MRSCont.files_ref{kk});
        MRSCont.raw_ref{kk}  = raw_ref;
    end
    if MRSCont.flags.hasWater
        raw_w   = io_loadspec_rda(MRSCont.files_w{kk});
        MRSCont.raw_w{kk}    = raw_w;
    end
end
fprintf('... done.\n');
toc(refLoadTime);

% Set flag
MRSCont.flags.coilsCombined     = 1;

end