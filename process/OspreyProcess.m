function [MRSCont] = OspreyProcess(MRSCont)
%% [MRSCont] = OspreyProcess(MRSCont)
%   This function pre-processes MRS data from all major vendors.
%   Data is read from the provided input filenames. It is shaped,
%   preprocessed, aligned, etc. according to the type of sequence
%   (un-edited data, MEGA-edited (ON/OFF), HERMES/HERCULES (A/B/C/D),
%   etc.).
%
%   USAGE:
%       MRSCont = OspreyProcess(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-19)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-02-19: First version of the code.

outputFolder = MRSCont.outputFolder;
fileID = fopen(fullfile(outputFolder, 'LogFile.txt'),'a+');
% Check that OspreyLoad has been run before
if ~MRSCont.flags.didLoadData
    msg = 'Trying to process data, but raw data has not been loaded yet. Run OspreyLoad first.';
    fprintf(fileID,msg);
    error(msg);
end


% Version, toolbox check and updating log file
MRSCont.ver.CheckPro        = '1.0.0 Pro';
fprintf(fileID,['Timestamp %s ' MRSCont.ver.Osp '  ' MRSCont.ver.CheckPro '\n'], datestr(now,'mmmm dd, yyyy HH:MM:SS'));

[~] = osp_Toolbox_Check ('OspreyProcess',MRSCont.flags.isGUI);

% Post-process raw data depending on sequence type
if ~MRSCont.flags.isPRIAM
    if MRSCont.flags.isUnEdited
        [MRSCont] = osp_processUnEdited(MRSCont);
    elseif MRSCont.flags.isMEGA
        [MRSCont] = osp_processMEGA(MRSCont);
    elseif MRSCont.flags.isHERMES
        [MRSCont] = osp_processHERMES(MRSCont);
    elseif MRSCont.flags.isHERCULES
        % For now, process HERCULES like HERMES data
        [MRSCont] = osp_processHERCULES(MRSCont);
    else
        msg = 'No flag set for sequence type!';
        fprintf(fileID,msg);
        error(msg);
    end
else
    fclose(fileID);
    [MRSCont] = osp_processMultiVoxel(MRSCont);
end



% Gather some more information from the processed data;
SubSpecNames = fieldnames(MRSCont.processed);
NoSubSpec = length(fieldnames(MRSCont.processed));
for ss = 1 : NoSubSpec
    for kk = 1 : MRSCont.nDatasets
            temp_sz(1,kk)= MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1);
            temp_sw(1,kk)= MRSCont.processed.(SubSpecNames{ss}){1,kk}.spectralwidth;
    end
    [MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint,MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint_ind,MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint_indsort]  = unique(temp_sz,'Stable');
    [MRSCont.info.(SubSpecNames{ss}).max_ndatapoint,MRSCont.info.(SubSpecNames{ss}).max_ndatapoint_ind] = max(temp_sz);
    temp_sw_store = temp_sw;
    for np = 1 : length(MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint)
        [max_ind] = find(temp_sz==MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint(np));
        temp_sw(max_ind) = nan;
        [MRSCont.info.(SubSpecNames{ss}).unique_spectralwidth{np},MRSCont.info.(SubSpecNames{ss}).unique_spectralwidth_ind{np},MRSCont.info.(SubSpecNames{ss}).unique_spectralwidth_indsort{np}]  = unique(temp_sw,'Stable');
        nanind = isnan(MRSCont.info.(SubSpecNames{ss}).unique_spectralwidth{np});
        MRSCont.info.(SubSpecNames{ss}).unique_spectralwidth{np}(isnan(MRSCont.info.(SubSpecNames{ss}).unique_spectralwidth{np}(1:end))) = [];
        MRSCont.info.(SubSpecNames{ss}).unique_spectralwidth_ind{np}(nanind ==1) = [];
        temp_sw = temp_sw_store;
    end
end
%% Clean up and save
% Set exit flags and reorder fields
MRSCont.flags.didProcess    = 1;
fclose(fileID);
[MRSCont]                   = osp_orderProcessFields(MRSCont);
MRSCont.ver.Pro             = '1.0.0 Pro';

% Save the output structure to the output folder
% Determine output folder
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

% Store data quality measures in csv file
if MRSCont.flags.isUnEdited
    names = {'NAA_SNR','NAA_FWHM','residual_water_ampl','freqShift'};
    subspec = {'A'};
    if MRSCont.flags.hasRef
        names = {'NAA_SNR','NAA_FWHM','water_FWHM','residual_water_ampl','freqShift'};
    end
elseif MRSCont.flags.isMEGA
    names = {'NAA_SNR','NAA_FWHM','residual_water_ampl','freqShift'};
    subspec = {'sum'};
    if MRSCont.flags.hasRef
        names = {'NAA_SNR','NAA_FWHM','water_FWHM','residual_water_ampl','freqShift'};
    end
elseif MRSCont.flags.isHERMES
    names = {'NAA_SNR','NAA_FWHM','residual_water_ampl','freqShift'};
    subspec = {'sum'};
    if MRSCont.flags.hasRef
        names = {'NAA_SNR','NAA_FWHM','water_FWHM','residual_water_ampl','freqShift'};
    end    
elseif MRSCont.flags.isHERCULES
    % For now, process HERCULES like HERMES data
    names = {'NAA_SNR','NAA_FWHM','residual_water_ampl','freqShift'};
    subspec = {'sum'};
    if MRSCont.flags.hasRef
        names = {'NAA_SNR','NAA_FWHM','water_FWHM','residual_water_ampl','freqShift'};
    end    
else
    msg = 'No flag set for sequence type!';
    fprintf(fileID,msg);
    error(msg);
end

if ~MRSCont.flags.isPRIAM
    if ~MRSCont.flags.hasRef
        QM = horzcat(MRSCont.QM.SNR.(subspec{1})',MRSCont.QM.FWHM.(subspec{1})',MRSCont.QM.res_water_amp.(subspec{1})',MRSCont.QM.freqShift.(subspec{1})');
    else
        QM = horzcat(MRSCont.QM.SNR.(subspec{1})',MRSCont.QM.FWHM.(subspec{1})',MRSCont.QM.FWHM.ref',MRSCont.QM.res_water_amp.(subspec{1})',MRSCont.QM.freqShift.(subspec{1})');
    end
    MRSCont.QM.tables = array2table(QM,'VariableNames',names);
    writetable(MRSCont.QM.tables,[outputFolder '/QM_processed_spectra.csv']);
end

% Optional:  Create all pdf figures
if MRSCont.opts.savePDF && ~MRSCont.flags.isPRIAM
    Names = fieldnames(MRSCont.processed);
    for kk = 1 : MRSCont.nDatasets
        for ss = 1 : length(Names)
            osp_plotModule(MRSCont, 'OspreyProcess', kk, Names{ss});
        end
    end
end

% Optional: write edited files to LCModel .RAW files
if MRSCont.opts.saveLCM && ~MRSCont.flags.isPRIAM
    [MRSCont] = osp_saveLCM(MRSCont);
end

% Optional: write edited files to jMRUI .txt files
if MRSCont.opts.savejMRUI && ~MRSCont.flags.isPRIAM
    [MRSCont] = osp_saveJMRUI(MRSCont);
end
% Optional: write edited files to vendor specific format files readable to
% LCModel and jMRUI
% SPAR/SDAT if Philips
% RDA if Siemens
if MRSCont.opts.saveVendor && ~MRSCont.flags.isPRIAM
    [MRSCont] = osp_saveVendor(MRSCont);
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end
