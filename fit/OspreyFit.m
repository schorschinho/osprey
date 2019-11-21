function [MRSCont] = OspreyFit(MRSCont)
%% [MRSCont] = OspreyFit(MRSCont)
%   This function performs spectral fitting on MRS data loaded previously
%   using OspreyLoad.
%
%   The method of fit, fitting range, included metabolites and other 
%   settings are set in the job file.
%
%   USAGE:
%       MRSCont = OspreyFit(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-24)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-02-24: First version of the code.

% Check that OspreyLoad has been run before
if ~MRSCont.flags.didLoadData
    error('Trying to fit data, but raw data has not been loaded yet. Run OspreyLoad first.')
end

% Check that OspreyProcess has been run before
if ~MRSCont.flags.didProcess
    error('Trying to fit data, but loaded data has not been process yet. Run OspreyProcess first.')
end

%% Load fit settings, prepare data and pass it on to the fitting algorithm
% Close any remaining open figures
close all;

% Initialise the fit - this step includes:
% - Parse the correct basis set
% - Apply settings on which metabolites/MM/lipids to include in the fit
% - Check for inconsistencies between basis set and data
[MRSCont] = osp_fitInitialise(MRSCont);

% Call the fit functions (depending on sequence type)
if MRSCont.flags.isUnEdited
    [MRSCont] = osp_fitUnEdited(MRSCont);
elseif MRSCont.flags.isMEGA
    [MRSCont] = osp_fitMEGA(MRSCont);
elseif MRSCont.flags.isHERMES
    [MRSCont] = osp_fitHERMES(MRSCont);
elseif MRSCont.flags.isHERCULES
    [MRSCont] = osp_fitHERCULES(MRSCont);
else
    error('No flag set for sequence type!');
end

%% Perform water reference and short-TE water fit

% If water reference exists, fit it
if MRSCont.flags.hasRef
    refFitTime = tic;
    reverseStr = '';
    % Loop over all the datasets here
    for kk = 1:MRSCont.nDatasets
        msg = sprintf('Fitting water reference from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        [MRSCont] = osp_fitWater(MRSCont, kk, 'ref');
    end
    fprintf('... done.\n');
    toc(refFitTime);
end

% If short TE water reference exists, fit it
if MRSCont.flags.hasWater
    waterFitTime = tic;
    reverseStr = '';
    % Loop over all the datasets here
    for kk = 1:MRSCont.nDatasets
        msg = sprintf('Fitting short-TE water from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        [MRSCont] = osp_fitWater(MRSCont, kk, 'w');
    end
    fprintf('... done.\n');
    toc(waterFitTime);
end

%% Clean up and save
% Set exit flags
MRSCont.flags.didFit           = 1;

end