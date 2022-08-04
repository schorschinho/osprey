function [MRSCont] = OspreySeg(MRSCont)
%% [MRSCont] = OspreySeg(MRSCont)
%   This function checks whether the structural image that the voxels were
%   coregistered to in OspreyCoreg has already been segmented by SPM12.
%
%   If it has not been, OspreySeg will call the SPM12 "New Segment"
%   function to perform segmentation into gray matter, white matter, and
%   CSF, and return fractional tissue volumes.
%
%   USAGE:
%       MRSCont = OspreySeg(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-08-21)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-08-21: First version of the code.

outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));

warning('off','all');
% Checking for version, toolbox, and previously run modules
[~,MRSCont.ver.CheckOsp ] = osp_CheckRunPreviousModule(MRSCont, 'OspreySeg');

% Set up SPM for batch processing
spm('defaults','fmri');
spm_jobman('initcfg');

% Set up saving location
saveDestinationSegMaps = fullfile(MRSCont.outputFolder, 'SegMaps');
if ~exist(saveDestinationSegMaps,'dir')
    mkdir(saveDestinationSegMaps);
end

%Do some check on the naming convention first to avoid overwriting the
%output due to non BIDS conform data names
switch MRSCont.vendor
    case {'Siemens', 'Philips'}
        % For Siemens and Philips data, this is simply the file that
        % is directly pointed to in the job file
        niftiFile = MRSCont.files_nii{1};
        if MRSCont.nDatasets(1) > 1
            niftiFile2 = MRSCont.files_nii{2};
        end
    case 'GE'
        % For GE data, SPM has created a *.nii file in the DICOM folder
        % that has been pointed to in the job file. We need to be
        % careful not to select potentially segmented files, so we'll
        % pick the filename that starts with an s (there should only be
        % one!).
        niftiList = dir([MRSCont.files_nii{1} filesep 's*.nii']);
        niftiFile = fullfile(MRSCont.files_nii{1}, niftiList.name);
        if MRSCont.nDatasets(1) > 1
            niftiList2 = dir([MRSCont.files_nii{2} filesep 's*.nii']);
            niftiFile2 = fullfile(MRSCont.files_nii{2}, niftiList2.name);
        end
    otherwise
        msg = 'Vendor not supported. Please contact the Osprey team (gabamrs@gmail.com).';
        fprintf(msg);
        error(msg);
end
SameName = 0;
if MRSCont.nDatasets(1) > 1
    [~, T1name, ~]  = fileparts(niftiFile);
    [~, T1name2, ~]  = fileparts(niftiFile2);
    SameName = strcmp(T1name,T1name2);
end

% Let's write the SPM12 propabilistic maps in a higher folder structure 
SepOutputFolder =  split(MRSCont.outputFolder, filesep);
SepOutputFolder(strcmp(SepOutputFolder,''))=[];
ind = find(strcmpi(SepOutputFolder,'derivatives')); %Do we have a folder named derivatives?
if ~isempty(ind) && ind ~= length(SepOutputFolder)
    goUpNTimes = length(SepOutputFolder) - ind;
    saveDestinationFilesSPM = fullfile(MRSCont.outputFolder,repmat(['..' filesep],[1 goUpNTimes]));
    tempDir = dir(saveDestinationFilesSPM);
    saveDestinationFilesSPM = fullfile(tempDir(1).folder,'FilesSPM');
else
    saveDestinationFilesSPM = fullfile(MRSCont.outputFolder,repmat(['..' filesep],[1 1]));
    tempDir = dir(saveDestinationFilesSPM);
    saveDestinationFilesSPM = fullfile(MRSCont.outputFolder, 'FilesSPM');
end

if ~exist(saveDestinationFilesSPM,'dir')
    mkdir(saveDestinationFilesSPM);
end
%% Loop over all datasets
refSegTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end
for kk = 1:MRSCont.nDatasets(1)
    [~] = printLog('OspreySeg',kk,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
    if ~(MRSCont.flags.didSeg == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'seg') && (kk > length(MRSCont.seg.tissue.fGM))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)

        %%% 1. CHECK WHETHER SEGMENTATION HAS BEEN RUN BEFORE %%%
        % First, we need to find the T1-NIfTI file again:
        switch MRSCont.vendor
            case {'Siemens', 'Philips'}
                % For Siemens and Philips data, this is simply the file that
                % is directly pointed to in the job file
                niftiFile = MRSCont.files_nii{kk};
            case 'GE'
                % For GE data, SPM has created a *.nii file in the DICOM folder
                % that has been pointed to in the job file. We need to be
                % careful not to select potentially segmented files, so we'll
                % pick the filename that starts with an s (there should only be
                % one!).
                niftiList = dir([MRSCont.files_nii{kk} filesep 's*.nii']);
                niftiFile = fullfile(MRSCont.files_nii{kk}, niftiList.name);
            otherwise
                msg = 'Vendor not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                fprintf(msg);
                error(msg);
        end
        
        if SameName
            [PreFix] = osp_generate_SubjectAndSessionPrefix(niftiFile,kk);
            PreFix = [PreFix '_'];
        else
            PreFix = '';
        end

        % Get the input file name
        [T1dir, T1name, T1extini]  = fileparts(niftiFile);
        if strcmp(T1extini,'.gz')
            T1name = strrep(T1name, '.nii','');
        end

        if ~isempty(MRSCont.files_seg) %Use external segmentation
            if length(MRSCont.files_seg{kk}) > 1
                segFileGM   = MRSCont.files_seg{kk}{1};
                segFileWM   = MRSCont.files_seg{kk}{2};
                segFileCSF  = MRSCont.files_seg{kk}{3};
                [~, ~, Segextini]  = fileparts(segFileGM);
                if strcmp(Segextini,'.gz')
                    gunzip(segFileGM);
                    gunzip(segFileWM);
                    gunzip(segFileCSF);
                    segFileGM   = strrep(segFileGM,'.gz','');
                    segFileWM   = strrep(segFileWM,'.gz','');
                    segFileCSF  = strrep(segFileCSF,'.gz','');
                end
                singleTissueSegFile = 0;
            else
                segFile4D   = MRSCont.files_seg{kk}{1};
                [~, ~, Segextini]  = fileparts(segFile4D);
                 if strcmp(Segextini,'.gz')
                    gunzip(segFile4D);
                    segFile4D   = strrep(segFile4D,'.gz','');
                end
                singleTissueSegFile = 1;
            end
        else
            singleTissueSegFile = 0;
            segFile               = fullfile(saveDestinationFilesSPM, [PreFix T1name '_seg8.mat']);
            if ~exist(segFile,'file') %Keep backward compabilities
                segFile               = fullfile(T1dir, [T1name '_seg8.mat']);
                if exist(segFile,'file')
                    movefile(fullfile(T1dir, ['c1' T1name '.nii.gz']),fullfile(saveDestinationFilesSPM,['c1_' PreFix T1name '_space-scanner_spm12_pseg.nii.gz']));
                    movefile(fullfile(T1dir, ['c2' T1name '.nii.gz']),fullfile(saveDestinationFilesSPM,['c2_' PreFix T1name '_space-scanner_spm12_pseg.nii.gz']));
                    movefile(fullfile(T1dir, ['c3' T1name '.nii.gz']),fullfile(saveDestinationFilesSPM,['c3_' PreFix T1name '_space-scanner_spm12_pseg.nii.gz']));
                    movefile(fullfile(T1dir, [T1name '_seg8.mat']),fullfile(saveDestinationFilesSPM,[PreFix T1name '_seg8.mat']));
                    segFile               = fullfile(saveDestinationFilesSPM, [PreFix T1name '_seg8.mat']);
                end
            end
            % If a GM-segmented file doesn't exist, start the segmentation
            if ~exist(segFile,'file')
                %Uncompress .nii.gz if needed
                if strcmp(T1extini,'.gz')
                    gunzip(niftiFile);
                    niftiFile = strrep(niftiFile,'.gz','');
                end
                T1ext = '.nii';
                createSegJob(niftiFile);
                movefile(fullfile(T1dir, ['c1' T1name T1ext]),fullfile(saveDestinationFilesSPM,['c1_' PreFix T1name '_space-scanner_spm12_pseg' T1ext]));
                movefile(fullfile(T1dir, ['c2' T1name T1ext]),fullfile(saveDestinationFilesSPM,['c2_' PreFix T1name '_space-scanner_spm12_pseg' T1ext]));
                movefile(fullfile(T1dir, ['c3' T1name T1ext]),fullfile(saveDestinationFilesSPM,['c3_' PreFix T1name '_space-scanner_spm12_pseg' T1ext]));
                movefile(fullfile(T1dir, ['y_' T1name T1ext]),fullfile(saveDestinationFilesSPM,[PreFix T1name '_spm12-transformation_field' T1ext]));
                movefile(fullfile(T1dir, [T1name '_seg8.mat']),fullfile(saveDestinationFilesSPM,[PreFix T1name '_seg8.mat']));
            else
                if strcmp(T1extini,'.gz')
                    gunzip(niftiFile);
                    niftiFile = strrep(niftiFile,'.gz','');
                    T1ext = '.nii';
                else
                    T1ext = T1extini;
                end
                if exist(fullfile(saveDestinationFilesSPM,['c1_' PreFix T1name '_space-scanner_spm12_pseg.nii.gz']),'file')
                    gunzip(fullfile(saveDestinationFilesSPM, ['c1_' PreFix T1name '_space-scanner_spm12_pseg' T1ext '.gz']));
                    gunzip(fullfile(saveDestinationFilesSPM, ['c2_' PreFix T1name '_space-scanner_spm12_pseg' T1ext '.gz']));
                    gunzip(fullfile(saveDestinationFilesSPM, ['c3_' PreFix T1name '_space-scanner_spm12_pseg' T1ext '.gz']));
                    if exist(fullfile(saveDestinationFilesSPM,[PreFix T1name '_spm12-transformation_field' T1ext '.gz']),'file')
                        gunzip(fullfile(saveDestinationFilesSPM, [PreFix T1name '_spm12-transformation_field' T1ext '.gz']));
                    end
                end
                T1ext = '.nii';
            end
        end

        %%% 2. CREATE MASKED TISSUE MAPS %%%
        % Define file names
        if isempty(MRSCont.files_seg) %Use SPM segmentation
            segFileGM   = fullfile(saveDestinationFilesSPM, ['c1_' PreFix T1name '_space-scanner_spm12_pseg' T1ext]);
            segFileWM   = fullfile(saveDestinationFilesSPM, ['c2_' PreFix T1name '_space-scanner_spm12_pseg' T1ext]);
            segFileCSF  = fullfile(saveDestinationFilesSPM, ['c3_' PreFix T1name '_space-scanner_spm12_pseg' T1ext]);
        end
        % Load volumes
        if ~singleTissueSegFile
            GMvol  = spm_vol(segFileGM);
            WMvol  = spm_vol(segFileWM);
            CSFvol = spm_vol(segFileCSF);
        else
           SEGvol  = spm_vol(segFile4D);
           if size(SEGvol,1)>1 % 4D volume
               TissueSegFile4D = 1;
           else % It's a 3D atlas with different numbers (we treat it as AAL atlas)
               TissueSegFile4D = 0;
           end


            if ~TissueSegFile4D
                % Get tissue maps from AAL atlas
                [GM,WM,CSF,AAL]=AAL_to_3TissueSeg(SEGvol);
                [AALpath,AALfilename,AALext] = fileparts(SEGvol.fname);
                %Store them in SPM12 standard
                GMvol.fname    = fullfile(AALpath,['c1' AALfilename AALext ]);
                GMvol.descrip  = ['GMmask based on AAL atlas'];
                GMvol.dim      = SEGvol.dim;
                GMvol.dt       = SEGvol.dt;
                GMvol.mat      = SEGvol.mat;
                GMvol          = spm_write_vol(GMvol, GM);

                WMvol.fname    = fullfile(AALpath,['c2' AALfilename AALext ]);
                WMvol.descrip  = ['WMmask based on AAL atlas'];
                WMvol.dim      = SEGvol.dim;
                WMvol.dt       = SEGvol.dt;
                WMvol.mat      = SEGvol.mat;
                WMvol          = spm_write_vol(WMvol, WM);

                CSFvol.fname    = fullfile(AALpath,['c3' AALfilename AALext ]);
                CSFvol.descrip  = ['CSFmask based on AAL atlas'];
                CSFvol.dim      = SEGvol.dim;
                CSFvol.dt       = SEGvol.dt;
                CSFvol.mat      = SEGvol.mat;
                CSFvol          = spm_write_vol(CSFvol, CSF);

                AALvol.fname    = fullfile(AALpath,['aal' AALfilename AALext ]);
                AALvol.descrip  = ['AAL atlas without L-R difference'];
                AALvol.dim      = SEGvol.dim;
                AALvol.dt       = SEGvol.dt;
                AALvol.mat      = SEGvol.mat;
                AALvol          = spm_write_vol(AALvol, AAL);
            end
        end

        %Loop over voxels (for DualVoxel)
        if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
            Voxels = 1;
        else
            Voxels = 2;
        end
        for rr = 1 : Voxels
            % Get voxel mask filename
            if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
                vol_mask = MRSCont.coreg.vol_mask{kk};
            else
                vol_mask = MRSCont.coreg.vol_mask{kk}{rr};
            end
            if ~exist(vol_mask.fname,'file') && exist(strrep(vol_mask.fname,'.nii','.nii.gz'),'file')
                gunzip(strrep(vol_mask.fname,'.nii','.nii.gz'));
            end
            [maskDir, maskName, maskExt] = fileparts(vol_mask.fname);

            % Get the input file name
            [~,filename,~]   = fileparts(MRSCont.files{kk});
            if MRSCont.coreg.SameName
                [PreFixMask] = osp_generate_SubjectAndSessionPrefix(MRSCont.files{kk},kk);
            else
                PreFixMask = '';
            end
            
            % <source_entities>[_space-<space>][_res-<label>][_label-<label>][_desc-<label>]_probseg.nii.gz
            % e.g.
            % sub-01_acq-press_space-individual_desc-dlpfc_label-GM_probseg.nii.gz
            saveName = [PreFixMask '_' filename '_space-scanner']; %CWDJ Check space.

            %Add voxel number for DualVoxel
            if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
                VoxelNum = '_Voxel-1';
            else
                VoxelNum = ['_Voxel-' num2str(rr)];
            end

            % GM
            vol_GMMask.fname    = fullfile(saveDestinationSegMaps, [saveName VoxelNum '_label-GM' maskExt]);
            vol_GMMask.descrip  = ['GMmasked_MRS_Voxel_Mask_' VoxelNum];
            vol_GMMask.dim      = vol_mask.dim;
            vol_GMMask.dt       = vol_mask.dt;
            vol_GMMask.mat      = vol_mask.mat;
            if ~singleTissueSegFile || ~TissueSegFile4D
                GM_voxmask_vol      = GMvol.private.dat(:,:,:) .* vol_mask.private.dat(:,:,:);
            else
                GM_voxmask_vol      = SEGvol(1).private.dat(:,:,:,1) .* vol_mask.private.dat(:,:,:);
            end
            vol_GMMask          = spm_write_vol(vol_GMMask, GM_voxmask_vol);

            % WM
            vol_WMMask.fname    = fullfile(saveDestinationSegMaps, [saveName VoxelNum '_label-WM' maskExt]);
            vol_WMMask.descrip  = ['WMmasked_MRS_Voxel_Mask_' VoxelNum];
            vol_WMMask.dim      = vol_mask.dim;
            vol_WMMask.dt       = vol_mask.dt;
            vol_WMMask.mat      = vol_mask.mat;
            if ~singleTissueSegFile || ~TissueSegFile4D
                WM_voxmask_vol      = WMvol.private.dat(:,:,:) .* vol_mask.private.dat(:,:,:);
            else
                WM_voxmask_vol      = SEGvol(2).private.dat(:,:,:,2) .* vol_mask.private.dat(:,:,:);
            end
            vol_WMMask          = spm_write_vol(vol_WMMask, WM_voxmask_vol);

            % CSF
            vol_CSFMask.fname   = fullfile(saveDestinationSegMaps, [saveName VoxelNum '_label-CSF' maskExt]);
            vol_CSFMask.descrip = ['CSFmasked_MRS_Voxel_Mask_' VoxelNum];
            vol_CSFMask.dim     = vol_mask.dim;
            vol_CSFMask.dt      = vol_mask.dt;
            vol_CSFMask.mat     = vol_mask.mat;
            if ~singleTissueSegFile || ~TissueSegFile4D
                CSF_voxmask_vol     = CSFvol.private.dat(:,:,:) .* vol_mask.private.dat(:,:,:);
            else
                CSF_voxmask_vol      = SEGvol(3).private.dat(:,:,:,3) .* vol_mask.private.dat(:,:,:);
            end
            vol_CSFMask         = spm_write_vol(vol_CSFMask, CSF_voxmask_vol);

            % Save volume structures in MRSCont
            if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
                MRSCont.seg.vol_GM{kk} = vol_GMMask;
                MRSCont.seg.vol_WM{kk} = vol_WMMask;
                MRSCont.seg.vol_CSF{kk} = vol_CSFMask;
            else
                MRSCont.seg.vol_GM{kk}{rr} = vol_GMMask;
                MRSCont.seg.vol_WM{kk}{rr} = vol_WMMask;
                MRSCont.seg.vol_CSF{kk}{rr} = vol_CSFMask;
            end
            % For MRSI data
            if MRSCont.flags.isMRSI
                GM_MRSI=spm_vol(GMvol);
                GM_MRSI=spm_read_vols(GM_MRSI);

                WM_MRSI=spm_vol(WMvol);
                WM_MRSI=spm_read_vols(WM_MRSI);

                CSF_MRSI=spm_vol(CSFvol);
                CSF_MRSI=spm_read_vols(CSF_MRSI);

                brain = (GM_MRSI > 0.5 | WM_MRSI > 0.5 | CSF_MRSI > 0.5);



                if ~exist(MRSCont.coreg.vol_image{kk}.fname,'file')
                    gunzip([MRSCont.coreg.vol_image{kk}.fname, '.gz']);
                end
                if ~exist(MRSCont.coreg.vol_mask{kk}.fname,'file')
                    gunzip([MRSCont.coreg.vol_mask{kk}.fname, '.gz']);
                end

                Vmask=spm_vol(MRSCont.coreg.vol_mask{kk}.fname);
                Vmask=spm_read_vols(Vmask);

                brain = brain .* Vmask;

                non_zero = zeros(size(brain,3),1);
                for i = 1 : size(brain,3)
                    if sum(sum(brain(:,:,i))) > 0
                        non_zero(i) = 1;
                    end

                end
                non_zero_slice = find(non_zero);
                brain_vox = brain(:,:,non_zero_slice(1):non_zero_slice(end));
                brain_vox = imresize3(double(brain_vox),[MRSCont.raw{kk}.nXvoxels,MRSCont.raw{kk}.nYvoxels,MRSCont.raw{kk}.nZvoxels]);
                brain_vox(brain_vox<(max(max(brain_vox))/200)) = 0;
                brain_vox(brain_vox > 0) = 1;
                brain_vox_rot = zeros(size(brain_vox,2),size(brain_vox,1),size(brain_vox,3));
                for i = 1 : size(brain_vox,3)
                	brain_vox_rot(:,:,i) = rot90(brain_vox(:,:,i));
                end
                MRSCont.mask{kk} = brain_vox_rot;
                if exist([MRSCont.coreg.vol_mask{kk}.fname, '.gz'],'file')
                    delete(MRSCont.coreg.vol_mask{kk}.fname);
                end
                if exist([MRSCont.coreg.vol_image{kk}.fname, '.gz'],'file')
                    delete(MRSCont.coreg.vol_image{kk}.fname);
                end
            end
            %%% 3. DETERMINE FRACTIONAL TISSUE VOLUMES %%%
            % Sum image intensities over the entire masked tissue specific volume
            GMsum  = sum(sum(sum(vol_GMMask.private.dat(:,:,:))));
            WMsum  = sum(sum(sum(vol_WMMask.private.dat(:,:,:))));
            CSFsum = sum(sum(sum(vol_CSFMask.private.dat(:,:,:))));

            % Save three plane image to container
            if MRSCont.flags.addImages
                [MRSCont.seg.img_montage{kk},MRSCont.seg.size_vox_t(kk)] = osp_extract_three_plane_image_seg(niftiFile, vol_mask,vol_GMMask,vol_WMMask,vol_CSFMask,MRSCont.coreg.voxel_ctr{kk},MRSCont.coreg.T1_max{kk});
            end

            % Apply the deformation field for overlay
            [MaskDir, MaskName, MaskExt]  = fileparts(vol_mask.fname);
            MaskNameSPM152 = strrep(MaskName,'space-scanner','space-spm152');
            DeformField = 0;
            if exist(fullfile(saveDestinationFilesSPM, [PreFix T1name '_spm12-transformation_field' T1ext]))
                DeformField = 1;
            end
            if ~exist(fullfile(MaskDir,[MaskNameSPM152 MaskExt '.gz'])) && exist(fullfile(saveDestinationFilesSPM, [PreFix T1name '_spm12-transformation_field' T1ext]))
                ApplyInverseDeformationField(fullfile(saveDestinationFilesSPM, [PreFix T1name '_spm12-transformation_field' T1ext]),vol_mask.fname)
                movefile(fullfile(MaskDir,['w' MaskName MaskExt]), fullfile(MaskDir,[MaskNameSPM152 MaskExt]));
            end

            %Compress nifit and delete uncompressed files
            gzip(vol_GMMask.fname);
            delete(vol_GMMask.fname);
            gzip(vol_WMMask.fname);
            delete(vol_WMMask.fname);
            gzip(vol_CSFMask.fname);
            delete(vol_CSFMask.fname);
            if exist(fullfile(saveDestinationFilesSPM,[PreFix T1name '_spm12-transformation_field' T1ext]),'file')
                gzip(fullfile(saveDestinationFilesSPM, [PreFix T1name '_spm12-transformation_field' T1ext]));
                delete(fullfile(saveDestinationFilesSPM, [PreFix T1name '_spm12-transformation_field' T1ext]));
            end
            if exist(fullfile(MaskDir,[MaskNameSPM152,MaskExt]),'file')
                gzip(fullfile(MaskDir,[MaskNameSPM152,MaskExt]));
                delete(fullfile(MaskDir,[MaskNameSPM152,MaskExt]));
            end
            if isempty(MRSCont.files_seg) %Standard SPM12 segmentation
                gzip(GMvol.fname);
                delete(GMvol.fname);
                gzip(WMvol.fname);
                delete(WMvol.fname);
                gzip(CSFvol.fname);
                delete(CSFvol.fname);
            elseif strcmp(Segextini,'.gz') %External segmentation
                if ~singleTissueSegFile
                    gzip(GMvol.fname);
                    delete(GMvol.fname);
                    gzip(WMvol.fname);
                    delete(WMvol.fname);
                    gzip(CSFvol.fname);
                    delete(CSFvol.fname);
                else
                    gzip(SEGvol(1).fname);
                    delete(SEGvol(1).fname);
                end

            end
            delete(vol_mask.fname);
            if strcmp(T1extini,'.gz')
                gzip(MRSCont.coreg.vol_image{kk}.fname)
                delete(MRSCont.coreg.vol_image{kk}.fname);
            end

            % Normalize
            fGM  = GMsum / (GMsum + WMsum + CSFsum);
            fWM  = WMsum / (GMsum + WMsum + CSFsum);
            fCSF = CSFsum / (GMsum + WMsum + CSFsum);

            % Save normalized fractional tissue volumes to MRSCont
            MRSCont.seg.tissue.fGM(kk,rr)  = fGM;
            MRSCont.seg.tissue.fWM(kk,rr)  = fWM;
            MRSCont.seg.tissue.fCSF(kk,rr) = fCSF;

        end
    end
end

% Calculate the overlap image
if DeformField
    for kk = 1:MRSCont.nDatasets(1)
        vol_mask = MRSCont.coreg.vol_mask{kk};
        [MaskDir, MaskName, MaskExt]  = fileparts(vol_mask.fname);
        MaskNameSPM152 = strrep(MaskName,'space-scanner','space-spm152');
        if ~exist(fullfile(MaskDir,[MaskNameSPM152, MaskExt]),'file') && exist(strrep(fullfile(MaskDir,[MaskNameSPM152, MaskExt]),'.nii','.nii.gz'),'file')
            gunzip(strrep(fullfile(MaskDir,[MaskNameSPM152, MaskExt]),'.nii','.nii.gz'));
        end
        if kk == 1
            Files{1} = fullfile(MaskDir,[MaskNameSPM152, MaskExt, ',1']);
        else
            Files{end+1} = fullfile(MaskDir,[MaskNameSPM152, MaskExt, ',1']);
        end
    end
    if isfield(MRSCont, 'exclude')
        Files(MRSCont.exclude)=[];
    end
    CalculateMaskOverlap(Files',[osp_RemovePreFix(MaskName) '_VoxelOverlap'],MaskDir);
    
    for kk = 1:MRSCont.nDatasets(1)
        vol_mask = MRSCont.coreg.vol_mask{kk};
        [MaskDir, MaskName, MaskExt]  = fileparts(vol_mask.fname);
        MaskNameSPM152 = strrep(MaskName,'space-scanner','space-spm152');
        gzip(fullfile(MaskDir,[MaskNameSPM152, MaskExt]));
        delete(fullfile(MaskDir,[MaskNameSPM152, MaskExt]));
    end
    gzip(fullfile(MaskDir,[osp_RemovePreFix(MaskName) '_VoxelOverlap.nii']));
    delete(fullfile(MaskDir,[osp_RemovePreFix(MaskName) '_VoxelOverlap.nii']));
    MRSCont.seg.overlapfile = fullfile(MaskDir,[osp_RemovePreFix(MaskName) '_VoxelOverlap.nii']);
end
time = toc(refSegTime);
[~] = printLog('done',time,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
MRSCont.runtime.Seg = time;
%% Create table and tsv file
tissueTypes = {'fGM','fWM','fCSF'};

%Loop over voxels (for DualVoxel)
for rr = 1 : Voxels
    tissue = horzcat(MRSCont.seg.tissue.fGM(:,rr),MRSCont.seg.tissue.fWM(:,rr),MRSCont.seg.tissue.fCSF(:,rr));
    MRSCont.seg.(['tables_Voxel_' num2str(rr)]) = array2table(tissue,'VariableNames',tissueTypes);
    MRSCont.seg.(['tables_Voxel_' num2str(rr)]) = addprop(MRSCont.seg.(['tables_Voxel_' num2str(rr)]), {'VariableLongNames'}, {'variable'}); % add long name to table properties

    % Populate descriptive fields of table for JSON export
    MRSCont.seg.(['tables_Voxel_' num2str(rr)]).Properties.CustomProperties.VariableLongNames{'fGM'} = 'Voxel fraction of grey matter';
    MRSCont.seg.(['tables_Voxel_' num2str(rr)]).Properties.VariableDescriptions{'fGM'} = 'Normalized fractional volume of grey matter: fGM  = GMsum / (GMsum + WMsum + CSFsum)';
    MRSCont.seg.(['tables_Voxel_' num2str(rr)]).Properties.VariableUnits{'fGM'} = 'arbitrary';

    MRSCont.seg.(['tables_Voxel_' num2str(rr)]).Properties.CustomProperties.VariableLongNames{'fWM'} = 'Voxel fraction of white matter';
    MRSCont.seg.(['tables_Voxel_' num2str(rr)]).Properties.VariableDescriptions{'fWM'} = 'Normalized fractional volume of white matter: fWM  = WMsum / (GMsum + WMsum + CSFsum)';
    MRSCont.seg.(['tables_Voxel_' num2str(rr)]).Properties.VariableUnits{'fWM'} = 'arbitrary';

    MRSCont.seg.(['tables_Voxel_' num2str(rr)]).Properties.CustomProperties.VariableLongNames{'fCSF'} = 'Voxel fraction of Cerebrospinal Fluid';
    MRSCont.seg.(['tables_Voxel_' num2str(rr)]).Properties.VariableDescriptions{'fCSF'} = 'Normalized fractional volume of Cerebrospinal Fluid: fCSF  = CSFsum / (GMsum + WMsum + CSFsum)';
    MRSCont.seg.(['tables_Voxel_' num2str(rr)]).Properties.VariableUnits{'fCSF'} = 'arbitrary';

    % Write the table to a file with json sidecar
    osp_WriteBIDsTable(MRSCont.seg.(['tables_Voxel_' num2str(rr)]), [saveDestinationSegMaps  filesep 'TissueFractions_Voxel_' num2str(rr)])
end

%% Clean up and save
% Set exit flags and version
MRSCont.flags.didSeg           = 1;
diary off

% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

% Optional:  Create all pdf figures
if MRSCont.opts.savePDF
    osp_plotAllPDF(MRSCont, 'OspreySeg')
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end




function createSegJob(T1file)

% Created with SPM12 batch manager (standard options)
spmhome = fileparts(which('spm'));
tpm = cellstr(spm_select('ExtFPList',fullfile(spmhome,'tpm'),'TPM.nii'));

matlabbatch{1}.spm.spatial.preproc.channel.vols = {[T1file ',1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = tpm(1);
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = tpm(2);
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = tpm(3);
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = tpm(4);
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = tpm(5);
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = tpm(6);
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];

spm_jobman('run',matlabbatch);

end


function [GM,WM,CSF,AAL]=AAL_to_3TissueSeg(vol_aseg_dseg)
% This function is only working for a specific AAL atlas. You have to adapt
% the lable numbers to make it work for your specific atlas.

GM = zeros(vol_aseg_dseg.dim);
WM = zeros(vol_aseg_dseg.dim);
CSF = zeros(vol_aseg_dseg.dim);
AAL = zeros(vol_aseg_dseg.dim);
aseg_dseg = spm_read_vols(vol_aseg_dseg);

% Sort AAL labels into gray matter, white matter, and CSF
% Create white matter mask
WM(aseg_dseg==2 | aseg_dseg == 41) = 1; %Cerebral-White-Matter 1
WM(aseg_dseg == 16) = 1;% Brain-Stem 2

% Create gray matter mask
GM(aseg_dseg==3 | aseg_dseg == 42) = 1;% Cerebral-Cortex 3
GM(aseg_dseg==8 | aseg_dseg == 47) = 1;% Cerebellum-Cortex 4
GM(aseg_dseg==10 | aseg_dseg == 49) = 1;% Thalamus-Proper* 5
GM(aseg_dseg==11 | aseg_dseg == 50) = 1;% Caudate 6
GM(aseg_dseg==12 | aseg_dseg == 51) = 1;% Putamen 7
GM(aseg_dseg==13 | aseg_dseg == 52) = 1;% Pallidum 8
GM(aseg_dseg==17 | aseg_dseg == 53) = 1;% Hippocampus 9
GM(aseg_dseg==18 | aseg_dseg == 54) = 1;% Amygdala 10
GM(aseg_dseg==26 | aseg_dseg == 58) = 1;% Accumbens-area 11
GM(aseg_dseg==27 | aseg_dseg == 59) = 1;% Substantia-Nigra 12
GM(aseg_dseg==28 | aseg_dseg == 60) = 1;% VentralDC 13
GM(aseg_dseg==172) = 1;% Vermis 14

% Create csf mask
CSF(aseg_dseg==4 | aseg_dseg == 43) = 1;%Lateral Ventricles 15
CSF(aseg_dseg==14) = 1;% 3rd Ventricle 16
CSF(aseg_dseg==15) = 1;% 4th Ventricle 17

% Create a single aal volume
AAL(aseg_dseg==2 | aseg_dseg == 41) = 1; %Cerebral-White-Matter
AAL(aseg_dseg==16) = 2;% Brain-Stem

AAL(aseg_dseg==3 | aseg_dseg == 42) = 3;% Cerebral-Cortex
AAL(aseg_dseg==8 | aseg_dseg == 47) = 4;% Cerebellum-Cortex
AAL(aseg_dseg==10 | aseg_dseg == 49) = 5;% Thalamus-Proper*
AAL(aseg_dseg==11 | aseg_dseg == 50) = 6;% Caudate
AAL(aseg_dseg==12 | aseg_dseg == 51) = 7;% Putamen
AAL(aseg_dseg==13 | aseg_dseg == 52) = 8;% Pallidum
AAL(aseg_dseg==17 | aseg_dseg == 53) = 9;% Hippocampus
AAL(aseg_dseg==18 | aseg_dseg == 54) = 10;% Amygdala
AAL(aseg_dseg==26 | aseg_dseg == 58) = 11;% Accumbens-area
AAL(aseg_dseg==27 | aseg_dseg == 59) = 12;% Substantia-Nigra
AAL(aseg_dseg==28 | aseg_dseg == 60) = 13;% VentralDC
AAL(aseg_dseg==172) = 14;% Vermis

AAL(aseg_dseg==4 | aseg_dseg == 43) = 15;%Lateral Ventricles
AAL(aseg_dseg==14) = 16;% 3rd Ventricle
AAL(aseg_dseg==15) = 17;% 4th Ventricle

end

function ApplyInverseDeformationField(TransformationField,maskFile)
    % Set the BB from spm152 template here
%     template = which('osprey/libraries/MRIcroGL/templates/spm152.nii');
%     template = spm_vol(template);
%     [BB,vx] = spm_get_bbox(template);
    BB = [-75.7625,-110.7625,-71.7625;...
            76.1549,77.2906,86.0546];
    vx = [0.7375,0.7375,0.7375];
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {TransformationField};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[maskFile ',1']};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = BB;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vx;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('run',matlabbatch);
end

function CalculateMaskOverlap(Files,outputName,outputDir)
    exp = '(';
    for kk = 1 : length(Files)
        exp = [exp 'i' num2str(kk)];
        if kk ~= length(Files)
            exp = [exp '+'];
        else
            exp = [exp ')/' num2str(length(Files))];
        end
    end
    matlabbatch{1}.spm.util.imcalc.input = Files;
    matlabbatch{1}.spm.util.imcalc.output = outputName;
    matlabbatch{1}.spm.util.imcalc.outdir = {outputDir};
    matlabbatch{1}.spm.util.imcalc.expression = exp;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 0;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch);
end
