function [MRSCont] = osp_AverageAllDatasetsAlongExtra(MRSCont)
%% [MRSCont] = osp_AverageAllDatasetsAlongExtra(MRSCont)
%   This function averages all files of the MRSContainer along the extra
%   dimension. This can be useful to combine several acquisitons from the
%   same voxel. SNR and linewidth will be recalculated. Drift plots will
%   not be accurate anymore!

% Checking for version, toolbox, and previously run modules
osp_CheckRunPreviousModule(MRSCont, 'OspreyFit');
outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));
%% Combining data
for kk = 1:MRSCont.nDatasets(1)
    MRSCont.processed.metab{1, kk}=op_average_extra(MRSCont.processed.metab{1, kk},1);
    
    if MRSCont.flags.hasMM
        MRSCont.processed.mm{1, kk}=op_average_extra(MRSCont.processed.mm{1, kk},1);
    end
    
    if MRSCont.flags.hasMMRef
        MRSCont.processed.mm_ref{1, kk}=op_average_extra(MRSCont.processed.mm_ref{1, kk},1);
    end
    
    if MRSCont.flags.hasRef
        MRSCont.processed.ref{1, kk}=op_average_extra(MRSCont.processed.ref{1, kk},1);
    end
    
    if MRSCont.flags.hasWater
        MRSCont.processed.w{1, kk}=op_average_extra(MRSCont.processed.w{1, kk},1);
    end
end

%% Recalculate QM 
for kk = 1:MRSCont.nDatasets(1)
    SubSpec = MRSCont.processed.metab{1, kk}.names;
    % Calculate some spectral quality metrics here;
    [MRSCont.processed.metab{1, kk},SNR] = op_get_Multispectra_SNR(MRSCont.processed.metab{1, kk});
    FWHM = op_get_Multispectra_LW(MRSCont.processed.metab{1, kk});
    for ss = 1 : length(SubSpec)
        MRSCont.QM.SNR.metab(1,kk,ss)    = SNR{ss};
        MRSCont.QM.FWHM.metab(1,kk,ss)   = FWHM(ss); % in Hz
    end
    if MRSCont.flags.hasMM
        SubSpec = MRSCont.processed.mm{1, kk}.names;
        [MRSCont.processed.mm{1, kk},SNR] = op_get_Multispectra_SNR(MRSCont.processed.mm{1, kk},1);
        FWHM = op_get_Multispectra_LW(MRSCont.processed.mm{1, kk});
        for ss = 1 : length(SubSpec)
            MRSCont.QM.SNR.mm(1,kk,ss)    = SNR{ss};
            MRSCont.QM.FWHM.mm(1,kk,ss)   = FWHM(ss);
        end
    end
    if MRSCont.flags.hasRef
        MRSCont.QM.SNR.ref(1,kk)    = op_getSNR(MRSCont.processed.ref{kk});
        MRSCont.QM.FWHM.ref(1,kk)   = op_getLW(MRSCont.processed.ref{kk},4.2,5.2);
        MRSCont.processed.ref{kk}.QC_names = {'water'};
    end
    if MRSCont.flags.hasWater
        MRSCont.QM.SNR.w(1,kk)    = op_getSNR(MRSCont.processed.w{kk});
        MRSCont.QM.FWHM.w(1,kk)   = op_getLW(MRSCont.processed.w{kk},4.2,5.2);
        MRSCont.processed.w{kk}.QC_names = {'water'};
    end
    if MRSCont.flags.hasMMRef
        MRSCont.QM.SNR.mm_ref(1,kk)    = op_getSNR(MRSCont.processed.mm_ref{kk});
        MRSCont.QM.FWHM.mm_ref(1,kk)   = op_getLW(MRSCont.processed.mm_ref{kk},4.2,5.2);
        MRSCont.processed.mm_ref{kk}.QC_names = {'water'};
    end

end

%% Clean up and save
diary off
% Store data quality measures in csv file
if MRSCont.flags.isUnEdited
    subspec = 1;
    name = 'A';
elseif MRSCont.flags.isMEGA
    subspec = 1;
    name = 'A';
elseif MRSCont.flags.isHERMES
    subspec = 7;
    name = 'sum';
elseif MRSCont.flags.isHERCULES
    subspec = 7;
    name = 'sum';
else
    msg = 'No flag set for sequence type!';
    fprintf(fileID,msg);
    error(msg);
end
names = {'NAA_SNR','NAA_FWHM','residual_water_ampl','freqShift'};
if MRSCont.flags.hasRef
    names = {'NAA_SNR','NAA_FWHM','water_FWHM','residual_water_ampl','freqShift'};
end

if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    if ~MRSCont.flags.hasRef
        QM = horzcat(MRSCont.QM.SNR.metab(1,:,subspec)',MRSCont.QM.FWHM.metab(1,:,subspec)',MRSCont.QM.res_water_amp.metab(1,:,subspec)',MRSCont.QM.freqShift.metab(1,:,subspec)');
    else
        QM = horzcat(MRSCont.QM.SNR.metab(1,:,subspec)',MRSCont.QM.FWHM.metab(1,:,subspec)',MRSCont.QM.FWHM.ref(1,:)',MRSCont.QM.res_water_amp.metab(1,:,subspec)',MRSCont.QM.freqShift.metab(1,:,subspec)');
    end
    MRSCont.QM.tables = array2table(QM,'VariableNames',names);

    MRSCont.QM.tables = addprop(MRSCont.QM.tables, {'VariableLongNames'}, {'variable'}); % add long name to table properties
    % Loop over field names to populate descriptive fields of table for JSON export
    for JJ = 1:length(names)
        switch names{JJ}
            case 'NAA_SNR'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'NAA_SNR'} = 'Signal to noise ratio of NAA';
                MRSCont.QM.tables.Properties.VariableDescriptions{'NAA_SNR'} = ['The maximum amplitude of the NAA peak divided by twice the standard deviation of the noise calculated from subspectrum ' name];
                MRSCont.QM.tables.Properties.VariableUnits{'NAA_SNR'} = 'arbitrary';
            case 'NAA_FWHM'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'NAA_FWHM'} = 'Full width at half maximum of NAA';
                MRSCont.QM.tables.Properties.VariableDescriptions{'NAA_FWHM'} = ['The width of the NAA peak at half the maximum amplitude calculated as the average of the FWHM of the data and the FWHM of a lorentzian fit calculated from subspectrum ' name];
                MRSCont.QM.tables.Properties.VariableUnits{'NAA_FWHM'} = 'Hz';
            case 'water_FWHM'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'water_FWHM'} = 'Full width at half maximum of reference water peak';
                MRSCont.QM.tables.Properties.VariableDescriptions{'water_FWHM'} = 'The width of the water peak at half the maximum amplitude calculated as the average of the FWHM of the data and the FWHM of a lorentzian fit';
                MRSCont.QM.tables.Properties.VariableUnits{'water_FWHM'} = 'Hz';
            case 'residual_water_ampl'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'residual_water_ampl'} = 'Residual water amplitude';
                MRSCont.QM.tables.Properties.VariableDescriptions{'residual_water_ampl'} = ['The amplitude of the water signal removed by the HSVD filter calculated from subspectrum ' name];
                MRSCont.QM.tables.Properties.VariableUnits{'residual_water_ampl'} = 'arbitrary';
            case 'freqShift'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'freqShift'} = 'Frequency shift';
                MRSCont.QM.tables.Properties.VariableDescriptions{'freqShift'} = ['Frequency shift calculated from the cross-correlation between spectrum ' name 'and the reference peaks (creatine and choline)'];
                MRSCont.QM.tables.Properties.VariableUnits{'freqShift'} = 'Hz';
        end
    end

    % Write .tsv file and .json sidecar
    osp_WriteBIDsTable(MRSCont.QM.tables, [outputFolder filesep 'QM_processed_spectra'])
end

% Save the output structure to the output folder
% Determine output folder
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end