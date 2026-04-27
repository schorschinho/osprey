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
diary(fullfile(outputFolder, 'LogFile.txt'));

% Checking for version, toolbox, and previously run modules
osp_CheckRunPreviousModule(MRSCont, 'OspreyProcess');
[~,MRSCont.ver.CheckOsp ] = osp_Toolbox_Check('OspreyProcess',MRSCont.flags.isGUI);


% Post-process raw data depending on sequence type
if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
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
        fprintf(msg);
        error(msg);
    end
else
    refProcessTime = tic;
    [MRSCont] = osp_processMultiVoxel(MRSCont);
    switch MRSCont.opts.MRSI.NuisanceRemoval.water.type
        case 'L2-basis'
            for kk = 1:MRSCont.nDatasets(1)
                [MRSCont.processed.A{kk}] = op_CSIRemoveLipids(MRSCont.processed.A{kk}, MRSCont.opts.MRSI.NuisanceRemoval.water.basisArguments);
                if MRSCont.flags.isMEGA 
                    [MRSCont.processed.diff1{kk}] = op_CSIRemoveLipids(MRSCont.processed.diff1{kk}, MRSCont.opts.MRSI.NuisanceRemoval.water.basisArguments);
                end
            end
    end
    switch MRSCont.opts.MRSI.NuisanceRemoval.lipid.type
        case 'L2-basis'
            for kk = 1:MRSCont.nDatasets(1)
                [MRSCont.processed.A{kk}] = op_CSIRemoveLipids(MRSCont.processed.A{kk}, MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments);
                if MRSCont.flags.isMEGA 
                    [MRSCont.processed.diff1{kk}] = op_CSIRemoveLipids(MRSCont.processed.diff1{kk}, MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments);
                end
            end
        case 'L2-mask'
            for kk = 1:MRSCont.nDatasets(1)
                [MRSCont.processed.A{kk}] = op_CSIRemoveLipids(MRSCont.processed.A{kk}, MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments,squeeze(MRSCont.seg.tissue.lip(kk,:,:,:)));
                if MRSCont.flags.isMEGA 
                    [MRSCont.processed.diff1{kk}] = op_CSIRemoveLipids(MRSCont.processed.diff1{kk}, MRSCont.opts.MRSI.NuisanceRemoval.lipid.basisArguments,squeeze(MRSCont.seg.tissue.lip(kk,:,:,:)));
                end
            end
    end
    for kk = 1 :MRSCont.nDatasets
        if MRSCont.flags.isUnEdited
            MRSCont.processed.A{kk}.refFWHM = cell2mat(MRSCont.processed.A{kk}.refFWHM);
            MRSCont.processed.A{kk}.refFWHM = cell2mat(MRSCont.processed.A{kk}.refShift);
        elseif MRSCont.flags.isMEGA           
            MRSCont.processed.A{kk}.refFWHM = cell2mat(MRSCont.processed.A{kk}.refFWHM);
            MRSCont.processed.A{kk}.refFWHM = cell2mat(MRSCont.processed.A{kk}.refShift);
        elseif MRSCont.flags.isHERMES
            MRSCont.processed.A{kk}.refFWHM = cell2mat(MRSCont.processed.A{kk}.refFWHM);
            MRSCont.processed.A{kk}.refFWHM = cell2mat(MRSCont.processed.A{kk}.refShift);
        elseif MRSCont.flags.isHERCULES
            MRSCont.processed.A{kk}.refFWHM = cell2mat(MRSCont.processed.A{kk}.refFWHM);
            MRSCont.processed.A{kk}.refFWHM = cell2mat(MRSCont.processed.A{kk}.refShift);
        end       
    end
    time = toc(refProcessTime);
    fprintf('\n... done.\n Elapsed time %f seconds\n',time);
    MRSCont.runtime.Proc = time;
end


% Gather some more information from the processed data;
SubSpecNames = fieldnames(MRSCont.processed);
NoSubSpec = length(fieldnames(MRSCont.processed));
for ss = 1 : NoSubSpec
    for kk = 1 : MRSCont.nDatasets
            temp_sz(1,kk)= MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1);
            temp_sz_sw{1,kk} = ['np_sw_' num2str(MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1)) '_' num2str(MRSCont.processed.(SubSpecNames{ss}){1,kk}.spectralwidth)];   
    end
    [MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint_spectralwidth,MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint_spectralwidth_ind,~]  = unique(temp_sz_sw,'Stable');
    [MRSCont.info.(SubSpecNames{ss}).max_ndatapoint,MRSCont.info.(SubSpecNames{ss}).max_ndatapoint_ind] = max(temp_sz);
end
%% If DualVoxel or MRSI we want to extract y-axis scaling
% Creates y-axis range to align the process plots between datasets

if MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI
    MRSCont.plot.processed.match = 1; % Scaling between datasets is turned off by default
else
    MRSCont.plot.processed.match = 0; % Scaling between datasets is turned off by default
end
MRSCont = osp_scale_yaxis(MRSCont,'OspreyProcess');
%% Clean up and save
% Set exit flags and reorder fields
MRSCont.flags.didProcess    = 1;
diary off
[MRSCont]                   = osp_orderProcessFields(MRSCont);

% Store data quality measures in csv file
if MRSCont.flags.isUnEdited
    names = {'NAA_SNR','NAA_FWHM','residual_water_ampl','freqShift'};
    subspec = {'A'};
    if MRSCont.flags.hasRef
        names = {'NAA_SNR','NAA_FWHM','water_FWHM','residual_water_ampl','freqShift'};
    end
elseif MRSCont.flags.isMEGA
    names = {'NAA_SNR','NAA_FWHM','residual_water_ampl','freqShift'};
    subspec = {'A'};
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

if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    if ~MRSCont.flags.hasRef
        QM = horzcat(MRSCont.QM.SNR.(subspec{1})',MRSCont.QM.FWHM.(subspec{1})',MRSCont.QM.res_water_amp.(subspec{1})',MRSCont.QM.freqShift.(subspec{1})');
    else
        QM = horzcat(MRSCont.QM.SNR.(subspec{1})',MRSCont.QM.FWHM.(subspec{1})',MRSCont.QM.FWHM.ref',MRSCont.QM.res_water_amp.(subspec{1})',MRSCont.QM.freqShift.(subspec{1})');
    end
    MRSCont.QM.tables = array2table(QM,'VariableNames',names);
    writetable(MRSCont.QM.tables,[outputFolder '/QM_processed_spectra.csv']);
end

% Write processed nii mrsi results
if MRSCont.flags.isMRSI
    outputFolderNii = fullfile(outputFolder,'nii-export','processed_raw');
    if ~exist(outputFolderNii,'dir')
        mkdir(outputFolderNii);
    end
    if MRSCont.flags.hasWater
        outputFolderNii = fullfile(outputFolder,'nii-export','processed_raw_w');
        if ~exist(outputFolderNii,'dir')
            mkdir(outputFolderNii);
        end
    end
    for kk = 1 : MRSCont.nDatasets 
        if MRSCont.raw{kk}.nZvoxels > 1 && ~MRSCont.opts.MRSI.pseudo3D 
            reorder = flip(1:MRSCont.raw{kk}.nZvoxels);
            shift = floor(MRSCont.processed.A{kk}.nZvoxels/2);
            for ll = 1 : MRSCont.processed.A{kk}.nZvoxels
                if MRSCont.flags.isUnEdited
                    ToExport = MRSCont.processed.A{kk};
                    ToExport.fids = squeeze(ToExport.fids(:,:,:,ll));
                    ToExport.specs = squeeze(ToExport.specs(:,:,:,ll));
                    ToExport.nZvoxels = 1;
                    ToExport.sz = size(ToExport.fids);
                    ToExport.dims.Zvoxels = 0;
    
                    VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                    ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice
                    shift = shift - 1;
                    nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw',['processed_raw_slice_' num2str(reorder(ll)) '.nii.gz']));
                end
                if MRSCont.flags.isMEGA
                    % Export off spectrum
                    ToExport = MRSCont.processed.A{kk};
                    ToExport.fids = squeeze(ToExport.fids(:,:,:,ll));
                    ToExport.specs = squeeze(ToExport.specs(:,:,:,ll));
                    ToExport.nZvoxels = 1;
                    ToExport.sz = size(ToExport.fids);
                    ToExport.dims.Zvoxels = 0;
    
                    VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                    ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice
                    nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw',['processed_raw_slice_' num2str(reorder(ll)) '_A.nii.gz']));
    
                    % Export diff spectrum
                    ToExport = MRSCont.processed.diff1{kk};
                    ToExport.fids = squeeze(ToExport.fids(:,:,:,ll));
                    ToExport.specs = squeeze(ToExport.specs(:,:,:,ll));
                    ToExport.nZvoxels = 1;
                    ToExport.sz = size(ToExport.fids);
                    ToExport.dims.Zvoxels = 0;
    
                    VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                    ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice
                    shift = shift - 1;
                    nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw',['processed_raw_slice_' num2str(reorder(ll)) '_diff1.nii.gz']));
                end
            end
        else % Single slice MRSI or pseudo 3D
                if MRSCont.flags.isUnEdited
                    % Export off spectrum
                    ToExport = MRSCont.processed.A{kk};
                    ToExport.fids=flip(ToExport.fids,2);
                    ToExport.specs=flip(ToExport.specs,2);
                    if (MRSCont.raw{kk}.nZvoxels > 1)
                        ToExport.fids = flip(ToExport.fids,length(ToExport.sz));
                        ToExport.specs = flip(ToExport.fids,length(ToExport.sz));
                    end
                    nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw',['processed_raw_A.nii.gz']));
                end
                if MRSCont.flags.isMEGA
                    % Export off spectrum
                    ToExport = MRSCont.processed.A{kk};
                    ToExport.fids=flip(ToExport.fids,1);
                    ToExport.specs=flip(ToExport.specs,1);
                    if (MRSCont.raw{kk}.nZvoxels > 1)
                        ToExport.fids = flip(ToExport.fids,length(ToExport.sz));
                        ToExport.specs = flip(ToExport.fids,length(ToExport.sz));
                    end
                    nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw',['processed_raw_A.nii.gz']));
        
                    % Export diff spectrum
                    ToExport = MRSCont.processed.diff1{kk};
                    ToExport.fids=flip(ToExport.fids,2);
                    ToExport.specs=flip(ToExport.specs,2);
                    if (MRSCont.raw{kk}.nZvoxels > 1)
                        ToExport.fids = flip(ToExport.fids,length(ToExport.sz));
                        ToExport.specs = flip(ToExport.fids,length(ToExport.sz));
                    end
                    nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw',['processed_raw_diff1.nii.gz']));
                end
        end
        if MRSCont.flags.hasWater
            if MRSCont.raw_w{kk}.nZvoxels > 1 && ~MRSCont.opts.MRSI.pseudo3D 
                shift = floor(MRSCont.processed.w{kk}.nZvoxels/2);
                for ll = 1 : MRSCont.processed.w{kk}.nZvoxels
                    ToExport = MRSCont.processed.w{kk};
                    ToExport.fids = squeeze(ToExport.fids(:,:,:,ll));
                    ToExport.specs = squeeze(ToExport.specs(:,:,:,ll));
                    ToExport.nZvoxels = 1;
                    ToExport.sz = size(ToExport.fids);
                    ToExport.dims.Zvoxels = 0;
    
                    VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                    ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice
                    
                    shift = shift - 1;
                    nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw_w',['processed_raw_w_slice_' num2str(reorder(ll)) '.nii.gz']));
                end
            else
                ToExport = MRSCont.processed.w{kk};
                ToExport.fids=flip(ToExport.fids,2);
                    ToExport.specs=flip(ToExport.specs,2);
                if (MRSCont.raw_w{kk}.nZvoxels > 1)
                    ToExport.fids = flip(ToExport.fids,length(ToExport.sz));
                    ToExport.specs = flip(ToExport.fids,length(ToExport.sz));
                end
                nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw_w',['processed_raw_w.nii.gz']));
            end
        end
    end
    if isfield(MRSCont.opts.MRSI,'MaxEcho')
        if MRSCont.opts.MRSI.MaxEcho.separate
            [right] = op_SeparateMaxEcho(MRSCont.processed.A{kk},'right',MRSCont.opts.MRSI.MaxEcho.tstart);
            MRSCont.processed.AFID{kk} = right;

            [left] = op_SeparateMaxEcho(MRSCont.processed.A{kk},'flipleft',MRSCont.opts.MRSI.MaxEcho.tstart);
            

            if MRSCont.opts.MRSI.MaxEcho.AdditionalPhasing
                if MRSCont.flags.isGUI
                    progressText = MRSCont.flags.inProgress;
                else
                    progressText = '';
                end
                XVox = MRSCont.raw{kk}.nXvoxels;
                YVox = MRSCont.raw{kk}.nYvoxels;
                ZVox = MRSCont.raw{kk}.nZvoxels;
                NVox = XVox*YVox*ZVox;
                vox = 1;
                [~] = printLog('OspreyMaxEcho',[kk,vox],[MRSCont.nDatasets, NVox],progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
                for z = 1 : ZVox
                    for x = 1 : XVox
                        for y = 1 : YVox
                            if ZVox <=1
                                raw = op_takeVoxel(MRSCont.processed.AFID{kk},[x y]);
                                rawLeft  = op_takeVoxel(left,[x y]);
                                rawRight  = op_takeVoxel(right,[x y]);
                            else
                                raw = op_takeVoxel(MRSCont.processed.AFID{kk},[x y z]);
                                rawLeft = op_takeVoxel(left,[x y z]);
                                rawRight  = op_takeVoxel(right,[x y]);
                            end
                            switch MRSCont.opts.MRSI.phase.type
                                case 'none'
                                    % Do nothing
                                case 'Cr-Cho'
                                    % Fit a double-Lorentzian to the Cr-Cho area, and phase the spectrum
                                    % with the negative phase of that fit
                                    [raw,ph]       = op_phaseCrCho(raw, 1);
                                    rawLeft = op_addphase(rawLeft,-ph);
                                    rawRight = op_addphase(rawRight,ph);
                                case 'auto_phase'
                                    [raw,ph]       = op_autophase(raw, MRSCont.opts.MRSI.phase.limits(1),MRSCont.opts.MRSI.phase.limits(2));
                                    rawLeft = op_addphase(rawLeft,-ph);
                                    rawRight = op_addphase(rawRight,ph);
                            end
                            if MRSCont.opts.MRSI.MaxEcho.AdditionalFreqAlign
                                switch MRSCont.opts.MRSI.MaxEcho.FreqAlign.type
                                    case 'CC'
                                        temp = raw;
                                        if MRSCont.opts.MRSI.FreqAlign.zerofill
                                            temp = op_zeropad(temp,4);
                                        end
                                        [refShift, ~] = osp_XReferencing(temp,MRSCont.opts.MRSI.MaxEcho.FreqAlign.frequencies,MRSCont.opts.MRSI.MaxEcho.FreqAlign.polarity,...
                                                                        MRSCont.opts.MRSI.MaxEcho.FreqAlign.lim,MRSCont.opts.MRSI.MaxEcho.FreqAlign.realpart);

                                    case 'CCwithLipRemoval'
                                        temp = raw;
                                        if MRSCont.opts.MRSI.MaxEcho.FreqAlign.zerofill
                                            temp = op_zeropad(temp,4);
                                        end
                                        noise = std(real(temp.specs(temp.ppm <= 0 & temp.ppm >= -2)));
                                        lipid = max(real(temp.specs(temp.ppm <= 1.9 & temp.ppm >= 0)));
                                        ratio = lipid/noise;
                                        if ratio > MRSCont.opts.MRSI.MaxEcho.FreqAlign.thresh
                                            temp = op_Wavlet_Filter(temp, -2, 1.85, 2, 10, 0);
                                        end
                                        temp = op_Wavlet_Filter(temp, -2, 4.2, 2, 10000, 0);
                                        [refShift, ~] = osp_XReferencing(temp,MRSCont.opts.MRSI.MaxEcho.FreqAlign.frequencies,MRSCont.opts.MRSI.MaxEcho.FreqAlign.polarity,...
                                                                        MRSCont.opts.MRSI.MaxEcho.FreqAlign.lim,MRSCont.opts.MRSI.MaxEcho.FreqAlign.realpart);
                                end
                                [raw]             = op_freqshift(raw,-refShift);            % Reference spectra by cross-correlation 
                                [rawLeft]         = op_freqshift(rawLeft,-refShift);            % Reference spectra by cross-correlation 
                                [rawRight]        = op_freqshift(rawRight,-refShift);            % Reference spectra by cross-correlation 

                            end
                            if ZVox <=1
                                MRSCont.processed.AFID{kk} = op_addVoxel(MRSCont.processed.AFID{kk},raw,[x y],1);
                                left = op_addVoxel(left,rawLeft,[x y],1);
                                right = op_addVoxel(right,rawRight,[x y],1);
                            else
                                MRSCont.processed.AFID{kk} = op_addVoxel(MRSCont.processed.AFID{kk},raw,[x y z],1);
                                left = op_addVoxel(left,rawLeft,[x y z],1);
                                right = op_addVoxel(right,rawRight,[x y z],1);
                            end
                            vox = vox + 1;
                            [~] = printLog('OspreyMaxEcho',[kk,vox],[MRSCont.nDatasets, NVox],progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
                        end
                    end
                end
                [~] = printLog('MRSIdone',0,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
            end
            MRSCont.processed.Asep{kk} = op_mergeextra(right,left,'echoside');
             % Export spectra
            ToExport = MRSCont.processed.AFID{kk};
            ToExport.fids=flip(ToExport.fids,2);
            ToExport.specs=flip(ToExport.specs,2);
            if (MRSCont.raw{kk}.nZvoxels > 1)
                ToExport.fids = flip(ToExport.fids,length(ToExport.sz));
                ToExport.specs = flip(ToExport.fids,length(ToExport.sz));
            end
            nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw',['processed_raw_A_right.nii.gz']));          

            
             % Export spectra
            ToExport = left;
            ToExport.fids=flip(ToExport.fids,2);
            ToExport.specs=flip(ToExport.specs,2);
            if (MRSCont.raw{kk}.nZvoxels > 1)
                ToExport.fids = flip(ToExport.fids,length(ToExport.sz));
                ToExport.specs = flip(ToExport.fids,length(ToExport.sz));
            end
            nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw',['processed_raw_A_left.nii.gz']));
            

        end
    end
    if  (MRSCont.opts.MRSI.cosemtics.GaussianLB  > 0) || (MRSCont.opts.MRSI.cosemtics.ZeroFillFactor > 0)
        outputFolderNii = fullfile(outputFolder,'nii-export','processed_raw_enhanced');
        if ~exist(outputFolderNii,'dir')
            mkdir(outputFolderNii);
        end
        for kk = 1 : MRSCont.nDatasets 
            if MRSCont.raw{kk}.nZvoxels > 1 && ~MRSCont.opts.MRSI.pseudo3D 
                reorder = flip(1:MRSCont.raw{kk}.nZvoxels);
                shift = floor(MRSCont.processed.A{kk}.nZvoxels/2);
                for ll = 1 : MRSCont.processed.A{kk}.nZvoxels
                    if MRSCont.flags.isUnEdited
                        ToExport = MRSCont.processed.A{kk};
                        ToExport.fids = squeeze(ToExport.fids(:,:,:,ll));
                        ToExport.specs = squeeze(ToExport.specs(:,:,:,ll));
                        ToExport.nZvoxels = 1;
                        ToExport.sz = size(ToExport.fids);
                        ToExport.dims.Zvoxels = 0;
        
                        VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                        ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice
                        shift = shift - 1;
                        nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw_enhanced',['processed_raw_slice_' num2str(reorder(ll)) '.nii.gz']));
                    end
                    if MRSCont.flags.isMEGA
                        % Export off spectrum
                        ToExport = MRSCont.processed.A{kk};
                        ToExport.fids = squeeze(ToExport.fids(:,:,:,ll));
                        ToExport.specs = squeeze(ToExport.specs(:,:,:,ll));
                        ToExport.nZvoxels = 1;
                        ToExport.sz = size(ToExport.fids);
                        ToExport.dims.Zvoxels = 0;
        
                        VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                        ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice
                        nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw_enhanced',['processed_raw_slice_' num2str(reorder(ll)) '_A.nii.gz']));
        
                        % Export diff spectrum
                        ToExport = MRSCont.processed.diff1{kk};
                        ToExport.fids = squeeze(ToExport.fids(:,:,:,ll));
                        ToExport.specs = squeeze(ToExport.specs(:,:,:,ll));
                        ToExport.nZvoxels = 1;
                        ToExport.sz = size(ToExport.fids);
                        ToExport.dims.Zvoxels = 0;
        
                        VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                        ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice
                        shift = shift - 1;
                        nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw_enhanced',['processed_raw_slice_' num2str(reorder(ll)) '_diff1.nii.gz']));
                    end
                end
            else % Single slice MRSI or pseudo 3D
                    if MRSCont.flags.isUnEdited
                        % Export off spectrum
                        ToExport = MRSCont.processed.A{kk};
                        if MRSCont.opts.MRSI.cosemtics.ZeroFillFactor > 0
                            ToExport = op_zeropad(ToExport,MRSCont.opts.MRSI.cosemtics.ZeroFillFactor,1);
                        end
                        if MRSCont.opts.MRSI.cosemtics.GaussianLB > 0
                            ToExport = op_filter(ToExport,MRSCont.opts.MRSI.cosemtics.GaussianLB);
                        end
                        ToExport.fids=flip(ToExport.fids,2);
                        ToExport.specs=flip(ToExport.specs,2);
                        if (MRSCont.raw{kk}.nZvoxels > 1)
                            ToExport.fids = flip(ToExport.fids,length(ToExport.sz));
                            ToExport.specs = flip(ToExport.fids,length(ToExport.sz));
                        end
                        nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw_enhanced',['processed_raw_A.nii.gz']));
                    end
                    if MRSCont.flags.isMEGA
                        % Export off spectrum
                        ToExport = MRSCont.processed.A{kk};
                        ToExport.fids=flip(ToExport.fids,1);
                        ToExport.specs=flip(ToExport.specs,1);
                        if (MRSCont.raw{kk}.nZvoxels > 1)
                            ToExport.fids = flip(ToExport.fids,length(ToExport.sz));
                            ToExport.specs = flip(ToExport.fids,length(ToExport.sz));
                        end
                        nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw_enhanced',['processed_raw_A.nii.gz']));
            
                        % Export diff spectrum
                        ToExport = MRSCont.processed.diff1{kk};
                        ToExport.fids=flip(ToExport.fids,2);
                        ToExport.specs=flip(ToExport.specs,2);
                        if (MRSCont.raw{kk}.nZvoxels > 1)
                            ToExport.fids = flip(ToExport.fids,length(ToExport.sz));
                            ToExport.specs = flip(ToExport.fids,length(ToExport.sz));
                        end
                        nii = io_writeniimrs(ToExport, fullfile(outputFolder,'nii-export','processed_raw_enhanced',['processed_raw_diff1.nii.gz']));
                    end
            end
        end       
    end
end

% Optional:  Create all pdf figures
if MRSCont.opts.savePDF
    osp_plotAllPDF(MRSCont, 'OspreyProcess')
end

% Optional: write edited files to LCModel .RAW files
if MRSCont.opts.saveLCM && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    [MRSCont] = osp_saveLCM(MRSCont);
end

% Optional: write edited files to jMRUI .txt files
if MRSCont.opts.savejMRUI && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    [MRSCont] = osp_saveJMRUI(MRSCont);
end
% Optional: write edited files to vendor specific format files readable to
% LCModel and jMRUI
% SPAR/SDAT if Philips
% RDA if Siemens
if MRSCont.opts.saveVendor && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    [MRSCont] = osp_saveVendor(MRSCont);
end

% Optional: write edited files to NIfTI-MRS format
if MRSCont.opts.saveNII && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    [MRSCont] = osp_saveNII(MRSCont);
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
