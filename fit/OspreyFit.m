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

outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));

if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end

% ----- Load fit settings and fit the metabolite data -----
% Checking for version, toolbox, and previously run modules
osp_CheckRunPreviousModule(MRSCont, 'OspreyFit');
[~,MRSCont.ver.CheckOsp ] = osp_Toolbox_Check('OspreyFit',MRSCont.flags.isGUI);
% Start timer
MRSCont.runtime.Fit = 0;

% Initialise the fit - this step includes:
% - Parse the correct basis set
% - Apply settings on which metabolites/MM/lipids to include in the fit
% - Check for inconsistencies between basis set and data
[MRSCont] = osp_fitInitialise(MRSCont);

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
    elseif MRSCont.flags.isDWMRS
        [MRSCont] = osp_fitDWMRS(MRSCont);
    else
        msg = 'No flag set for sequence type!';
        fprintf(msg);
        error(msg);
    end
else
    [MRSCont] = osp_fitMultiVoxel(MRSCont);
end


% ----- Perform water reference and short-TE water fit -----
% The water signal is automatically integrated when the LCModel fit option is
% being used. In Osprey, we explicitly model the water data with a
% dedicated simulated water basis function.
if strcmpi(MRSCont.opts.fit.method, 'Osprey')

    % If water reference exists, fit it
    if MRSCont.flags.hasRef
        refFitTime = tic;
        % Loop over all the datasets here
        for kk = 1:MRSCont.nDatasets
            [~] = printLog('OspreyFitRef', kk, MRSCont.nDatasets, progressText, MRSCont.flags.isGUI, MRSCont.flags.isMRSI);
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
        % Loop over all the datasets here
        for kk = 1:MRSCont.nDatasets
            [~] = printLog('OspreyFitWater',kk,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
            if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.w.fitParams))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
                [MRSCont] = osp_fitWater(MRSCont, kk, 'w');
            end
        end
        time = toc(waterFitTime);
        fprintf('... done.\n Elapsed time %f seconds\n',time);
        MRSCont.runtime.FitWater = time;
        MRSCont.runtime.Fit = MRSCont.runtime.Fit + time;
    end

    MRSCont.runtime.Fit = MRSCont.runtime.Fit + MRSCont.runtime.FitMet;
    [~] = printLog('Fulldone',MRSCont.runtime.Fit,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);

end


%% If DualVoxel or MRSI we want to extract y-axis scaling
if MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI
    MRSCont = osp_scale_yaxis(MRSCont,'OspreyLoad');
    MRSCont.fit.resBasisSet = MRSCont.fit.resBasisSet{2,2};
end

%% Store and print some QM parameters
if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI && ~MRSCont.flags.isDWMRS
    [MRSCont] = osp_fit_Quality(MRSCont);

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
%Delete redundant resBasiset entries
if strcmpi(MRSCont.opts.fit.method, 'Osprey')
    if ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
        FitNames = fieldnames(MRSCont.fit.results);
        NoFit = length(fieldnames(MRSCont.fit.results));
        for sf = 1 : NoFit
            if iscell(MRSCont.fit.resBasisSet.(FitNames{sf}))
                MRSCont.fit.resBasisSet.(FitNames{sf}) = MRSCont.fit.resBasisSet.(FitNames{sf})(MRSCont.info.A.unique_ndatapoint_spectralwidth_ind);
                for combs = 1 : length(MRSCont.info.A.unique_ndatapoint_spectralwidth_ind)
                    resBasisSetNew.(FitNames{sf}).([MRSCont.info.A.unique_ndatapoint_spectralwidth{combs}]) = MRSCont.fit.resBasisSet.(FitNames{sf}){combs};
                end
            else
                MRSCont.fit.resBasisSet.(FitNames{sf}).water = MRSCont.fit.resBasisSet.(FitNames{sf}).water(MRSCont.info.(FitNames{sf}).unique_ndatapoint_spectralwidth_ind);
                for combs = 1 : length(MRSCont.info.(FitNames{sf}).unique_ndatapoint_spectralwidth_ind)
                    resBasisSetNew.(FitNames{sf}).water.([MRSCont.info.(FitNames{sf}).unique_ndatapoint_spectralwidth{combs}]) = MRSCont.fit.resBasisSet.(FitNames{sf}).water{combs};
                end
            end
        end
        MRSCont.fit.resBasisSet = resBasisSetNew;
    end
end

% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

% Optional:  Create all pdf figures
if MRSCont.opts.savePDF
    osp_plotAllPDF(MRSCont, 'OspreyFit');
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end
