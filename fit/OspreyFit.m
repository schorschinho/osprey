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

outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));
% Check that OspreyLoad has been run before
if ~MRSCont.flags.didLoadData
    msg = 'Trying to fit data, but raw data has not been loaded yet. Run OspreyLoad first.';
    fprintf(msg);
    error(msg);
end

% Check that OspreyProcess has been run before
if ~MRSCont.flags.didProcess
    msg = 'Trying to fit data, but loaded data has not been process yet. Run OspreyProcess first.';
    fprintf(msg);
    error(msg);
end

%% Load fit settings, prepare data and pass it on to the fitting algorithm

% Version, toolbox check and updating log file
[~,MRSCont.ver.CheckOsp ] = osp_Toolbox_Check ('OspreyFit',MRSCont.flags.isGUI);
MRSCont.runtime.Fit = 0;

% Initialise the fit - this step includes:
% - Parse the correct basis set
% - Apply settings on which metabolites/MM/lipids to include in the fit
% - Check for inconsistencies between basis set and data
[MRSCont] = osp_fitInitialise(MRSCont);
MRSCont.opts.fit.outputFolder = outputFolder;
% Call the fit functions (depending on sequence type)
if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
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
        msg = 'No flag set for sequence type!';
        fprintf(msg);
        error(msg);
    end
else
    [MRSCont] = osp_fitMultiVoxel(MRSCont);
end

%% Perform water reference and short-TE water fit
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
end
% If water reference exists, fit it
if MRSCont.flags.hasRef
    refFitTime = tic;
    reverseStr = '';
    % Loop over all the datasets here
    for kk = 1:MRSCont.nDatasets
        msg = sprintf('\nFitting water reference from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, '\n', msg]);
        if MRSCont.flags.isGUI        
            set(progressText,'String' ,sprintf('Fitting water reference from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets));
            drawnow
        end
        if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.ref.fitParams))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
            [MRSCont] = osp_fitWater(MRSCont, kk, 'ref');
        end
    end
    time = toc(refFitTime);
    if MRSCont.flags.isGUI        
        set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
        pause(1);
    end
    fprintf('... done.\n Elapsed time %f seconds\n',time);
    MRSCont.runtime.FitRef = time;
    MRSCont.runtime.Fit = MRSCont.runtime.Fit + time;
end

% If short TE water reference exists, fit it
if MRSCont.flags.hasWater
    waterFitTime = tic;
    reverseStr = '';   
    % Loop over all the datasets here
    for kk = 1:MRSCont.nDatasets
        msg = sprintf('\nFitting short-TE water from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf([reverseStr, '\n', msg]);
        if MRSCont.flags.isGUI        
            set(progressText,'String' ,sprintf('Fitting short-TE water from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets));
            drawnow
        end
        if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.w.fitParams))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
            [MRSCont] = osp_fitWater(MRSCont, kk, 'w');
        end
    end
    time = toc(waterFitTime);
    if MRSCont.flags.isGUI        
        set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
        pause(1);
    end
    fprintf('... done.\n Elapsed time %f seconds\n',time);
    MRSCont.runtime.FitWater = time;
    MRSCont.runtime.Fit = MRSCont.runtime.Fit + time;
end
MRSCont.runtime.Fit = MRSCont.runtime.Fit + MRSCont.runtime.FitMet;
fprintf('Full fit time %f seconds\n',MRSCont.runtime.Fit);

%% If DualVoxel or MRSI we want to extract y-axis scaling
if MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI
    osp_scale_yaxis(MRSCont,'OspreyLoad');   
end
%% Store  and print some QM parameters
if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    [MRSCont]=osp_fit_Quality(MRSCont);

    % Store data quality measures in csv file
    if MRSCont.flags.isUnEdited
        relResA = MRSCont.QM.relAmpl.A';
        MRSCont.QM.tables.relResA = relResA;
    elseif MRSCont.flags.isMEGA
        if strcmp( MRSCont.opts.fit.style, 'Separate')
            relResA = MRSCont.QM.relAmpl.A';
            relResdiff1 = MRSCont.QM.relAmpl.diff1';
            MRSCont.QM.tables.relResA = relResA;
            MRSCont.QM.tables.relResdiff1 = relResdiff1;
        else
            relRessum = MRSCont.QM.relAmpl.sum';
            relResdiff1 = MRSCont.QM.relAmpl.diff1';
            MRSCont.QM.tables.relRessum = relRessum;
            MRSCont.QM.tables.relResdiff1 = relResdiff1;        
        end        
    elseif MRSCont.flags.isHERMES
            relRessum = MRSCont.QM.relAmpl.sum';
            relResdiff1 = MRSCont.QM.relAmpl.diff1';
            relResdiff2 = MRSCont.QM.relAmpl.diff2';
            MRSCont.QM.tables.relRessum = relRessum;
            MRSCont.QM.tables.relResdiff1 = relResdiff1; 
            MRSCont.QM.tables.relResdiff2 = relResdiff2;
    elseif MRSCont.flags.isHERCULES
        % For now, process HERCULES like HERMES data
            relRessum = MRSCont.QM.relAmpl.sum';
            relResdiff1 = MRSCont.QM.relAmpl.diff1';
            relResdiff2 = MRSCont.QM.relAmpl.diff2';
            MRSCont.QM.tables.relRessum = relRessum;
            MRSCont.QM.tables.relResdiff1 = relResdiff1; 
            MRSCont.QM.tables.relResdiff2 = relResdiff2;    
    else
        msg = 'No flag set for sequence type!';
        fprintf(msg);
        error(msg);
    end

    writetable(MRSCont.QM.tables,[outputFolder '/QM_processed_spectra.csv']);
end

%% Clean up and save
% Set exit flags and version
MRSCont.flags.didFit           = 1;

diary off
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

% Optional:  Create all pdf figures
if MRSCont.opts.savePDF
    osp_plotAllPDF(MRSCont, 'OspreyFit')
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end