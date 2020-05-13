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

% Version check
MRSCont.ver.CheckFit             = '100 Fit';

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
    % For now, fit HERCULES like HERMES data
    [MRSCont] = osp_fitHERCULES(MRSCont);
else
    error('No flag set for sequence type!');
end

%% Perform water reference and short-TE water fit

% If water reference exists, fit it
if MRSCont.flags.hasRef
    refFitTime = tic;
    reverseStr = '';
    if MRSCont.flags.isGUI
        progressbar = waitbar(0,'Start','Name','Osprey Fit');
        waitbar(0,progressbar,sprintf('Fitted water reference from dataset %d out of %d total datasets...\n', 0, MRSCont.nDatasets))
    end
    % Loop over all the datasets here
    for kk = 1:MRSCont.nDatasets
        msg = sprintf('\nFitting water reference from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        if ((MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.ref.fitParams))) || ~isfield(MRSCont.ver, 'Fit') || ~strcmp(MRSCont.ver.Fit,MRSCont.ver.CheckFit))
            [MRSCont] = osp_fitWater(MRSCont, kk, 'ref');
        end
    end
    if MRSCont.flags.isGUI        
        waitbar(kk/MRSCont.nDatasets,progressbar,sprintf('Fitted water reference from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets))
        waitbar(1,progressbar,'...done')
        close(progressbar)
    end
    fprintf('... done.\n');
    toc(refFitTime);
end

% If short TE water reference exists, fit it
if MRSCont.flags.hasWater
    waterFitTime = tic;
    reverseStr = '';
    if MRSCont.flags.isGUI
        progressbar = waitbar(0,'Start','Name','Osprey Fit');
        waitbar(0,progressbar,sprintf('Fitted short-TE water from dataset %d out of %d total datasets...\n', 0, MRSCont.nDatasets))
    end    
    % Loop over all the datasets here
    for kk = 1:MRSCont.nDatasets
        msg = sprintf('\nFitting short-TE water from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        if ((MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.w.fitParams))) || ~isfield(MRSCont.ver, 'Fit') || ~strcmp(MRSCont.ver.Fit,MRSCont.ver.CheckFit))
            [MRSCont] = osp_fitWater(MRSCont, kk, 'w');
        end
    end
    if MRSCont.flags.isGUI        
        waitbar(kk/MRSCont.nDatasets,progressbar,sprintf('Fitted short-TE water from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets))
        waitbar(1,progressbar,'...done')
        close(progressbar)
    end
    fprintf('... done.\n');
    toc(waterFitTime);
end

%% Clean up and save
% Set exit flags and version
MRSCont.flags.didFit           = 1;
MRSCont.ver.Fit            = '100 Fit';

% Delete redundant resBasiset entries
% FitNames = fieldnames(MRSCont.fit.results);
% NoFit = length(fieldnames(MRSCont.fit.results));
% for sf = 1 : NoFit
%     if iscell(MRSCont.fit.resBasisSet.(FitNames{sf}))
%         MRSCont.fit.resBasisSet.(FitNames{sf}) = MRSCont.fit.resBasisSet.(FitNames{sf})(MRSCont.info.A.unique_ndatapoint_ind);
%     else
%         MRSCont.fit.resBasisSet.(FitNames{sf}).water = MRSCont.fit.resBasisSet.(FitNames{sf}).water(MRSCont.info.(FitNames{sf}).unique_ndatapoint_ind); 
%     end
% end

% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

if ~MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont');
end

end