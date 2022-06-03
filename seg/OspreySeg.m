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
saveDestination = fullfile(MRSCont.outputFolder, 'SegMaps');
if ~exist(saveDestination,'dir')
    mkdir(saveDestination);
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
            segFile               = fullfile(T1dir, [T1name '_seg8.mat']);
            % If a GM-segmented file doesn't exist, start the segmentation
            if ~exist(segFile,'file')
                %Uncompress .nii.gz if needed
                if strcmp(T1extini,'.gz')
                    gunzip(niftiFile);
                    niftiFile = strrep(niftiFile,'.gz','');
                end
                T1ext = '.nii';
                createSegJob(niftiFile);
            else
                if strcmp(T1extini,'.gz')
                    gunzip(niftiFile);
                    niftiFile = strrep(niftiFile,'.gz','');
                    T1ext = '.nii';
                else
                    T1ext = T1extini;
                end
                if exist(fullfile(T1dir, ['c1' T1name '.nii.gz']),'file')
                    gunzip(fullfile(T1dir, ['c1' T1name T1ext '.gz']));
                    gunzip(fullfile(T1dir, ['c2' T1name T1ext '.gz']));
                    gunzip(fullfile(T1dir, ['c3' T1name T1ext '.gz']));
                end
                T1ext = '.nii';
            end
        end

        %%% 2. CREATE MASKED TISSUE MAPS %%%
        % Define file names
        if isempty(MRSCont.files_seg) %Use SPM segmentation
            segFileGM   = fullfile(T1dir, ['c1' T1name T1ext]);
            segFileWM   = fullfile(T1dir, ['c2' T1name T1ext]);
            segFileCSF  = fullfile(T1dir, ['c3' T1name T1ext]);
        end
        % Load volumes
        if ~singleTissueSegFile
            GMvol  = spm_vol(segFileGM);
            WMvol  = spm_vol(segFileWM);
            CSFvol = spm_vol(segFileCSF);
        else
           SEGvol  = spm_vol(segFile4D);
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

            % <source_entities>[_space-<space>][_res-<label>][_label-<label>][_desc-<label>]_probseg.nii.gz
            % e.g.
            % sub-01_acq-press_space-individual_desc-dlpfc_label-GM_probseg.nii.gz
            saveName = [osp_RemoveSuffix(filename),'_space-scanner']; %CWDJ Check space.

            %Add voxel number for DualVoxel
            if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
                VoxelNum = '_Voxel-1';
            else
                VoxelNum = ['_Voxel-' num2str(rr)];
            end

            % GM
            vol_GMMask.fname    = fullfile(saveDestination, [saveName VoxelNum '_label-GM' maskExt]);
            vol_GMMask.descrip  = ['GMmasked_MRS_Voxel_Mask_' VoxelNum];
            vol_GMMask.dim      = vol_mask.dim;
            vol_GMMask.dt       = vol_mask.dt;
            vol_GMMask.mat      = vol_mask.mat;
            if ~singleTissueSegFile
                GM_voxmask_vol      = GMvol.private.dat(:,:,:) .* vol_mask.private.dat(:,:,:);
            else
                GM_voxmask_vol      = SEGvol(1).private.dat(:,:,:,1) .* vol_mask.private.dat(:,:,:);
            end
            vol_GMMask          = spm_write_vol(vol_GMMask, GM_voxmask_vol);

            % WM
            vol_WMMask.fname    = fullfile(saveDestination, [saveName VoxelNum '_label-WM' maskExt]);
            vol_WMMask.descrip  = ['WMmasked_MRS_Voxel_Mask_' VoxelNum];
            vol_WMMask.dim      = vol_mask.dim;
            vol_WMMask.dt       = vol_mask.dt;
            vol_WMMask.mat      = vol_mask.mat;
            if ~singleTissueSegFile
                WM_voxmask_vol      = WMvol.private.dat(:,:,:) .* vol_mask.private.dat(:,:,:);
            else
                WM_voxmask_vol      = SEGvol(2).private.dat(:,:,:,2) .* vol_mask.private.dat(:,:,:);
            end
            vol_WMMask          = spm_write_vol(vol_WMMask, WM_voxmask_vol);

            % CSF
            vol_CSFMask.fname   = fullfile(saveDestination, [saveName VoxelNum '_label-CSF' maskExt]);
            vol_CSFMask.descrip = ['CSFmasked_MRS_Voxel_Mask_' VoxelNum];
            vol_CSFMask.dim     = vol_mask.dim;
            vol_CSFMask.dt      = vol_mask.dt;
            vol_CSFMask.mat     = vol_mask.mat;
            if ~singleTissueSegFile
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

            %Compress nifit and delete uncompressed files
            gzip(vol_GMMask.fname);
            delete(vol_GMMask.fname);
            gzip(vol_WMMask.fname);
            delete(vol_WMMask.fname);
            gzip(vol_CSFMask.fname);
            delete(vol_CSFMask.fname);
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
                if ~exist([MRSCont.coreg.vol_image{kk}.fname,'.gz'],'file')
                    gzip(MRSCont.coreg.vol_image{kk}.fname)
                    delete(MRSCont.coreg.vol_image{kk}.fname);
                end
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


% If a pair of secondary images exists, start the process of
% segmenting the second T1 here.
% We will keep the segmented tissue maps in native space, so we don't have
% to reslice them. We can then use the resliced voxel mask from the
% OspreyCoreg to
for kk = 1:MRSCont.nDatasets
    if (MRSCont.flags.hasSecondT1 && MRSCont.flags.hasPET)

        %%% 4. CHECK WHETHER SEGMENTATION HAS BEEN RUN BEFORE %%%
        % Get the file name of the second T1 from the SPM volume stored in
        % the MRS container
        [T1_2nd_dir, T1_2nd_name, T1_2nd_ext]   = fileparts(MRSCont.coreg.vol_image_2nd{kk}.fname);
        segFileGM2nd                            = fullfile(T1_2nd_dir, ['c1' T1_2nd_name T1_2nd_ext]);
        % If a GM-segmented file doesn't exist, start the segmentation
        if ~exist(segFileGM2nd,'file')
            createSegJob(fullfile(T1_2nd_dir, [T1_2nd_name T1_2nd_ext]));
        end


        %%% 5. CREATE MASKED TISSUE MAPS %%%
        % Define file names
        segFileWM2nd   = fullfile(T1_2nd_dir, ['c2' T1_2nd_name T1_2nd_ext]);
        segFileCSF2nd   = fullfile(T1_2nd_dir, ['c3' T1_2nd_name T1_2nd_ext]);
        % Load volumes
        GMvol2nd  = spm_vol(segFileGM2nd);
        WMvol2nd  = spm_vol(segFileWM2nd);
        CSFvol2nd  = spm_vol(segFileCSF2nd);
        % Get the filename of the resliced voxel mask from the SPM volume
        % in the MRSCont container
        vol_mask_2nd = MRSCont.coreg.vol_mask_2nd{kk};
        [maskDir2nd, maskName2nd, maskExt2nd] = fileparts(vol_mask_2nd.fname);

        % Create and save masked tissue maps
        % Get the input file name
        [path,filename,~]   = fileparts(MRSCont.files{kk});
        % For batch analysis, get the last two sub-folders (e.g. site and
        % subject)
        path_split          = regexp(path,filesep,'split');
        if length(path_split) > 2
            saveName = [path_split{end-1} '_' path_split{end} '_' filename];
        end
        % GM
        vol_GMMask.fname    = fullfile(saveDestination, [saveName '_GM_2ndT1' maskExt2nd]);
        vol_GMMask.descrip  = 'GMmasked_MRS_Voxel_Mask';
        vol_GMMask.dim      = vol_mask_2nd.dim;
        vol_GMMask.dt       = vol_mask_2nd.dt;
        vol_GMMask.mat      = vol_mask_2nd.mat;
        GM_voxmask_vol      = GMvol2nd.private.dat(:,:,:) .* vol_mask_2nd.private.dat(:,:,:);
        vol_GMMask          = spm_write_vol(vol_GMMask, GM_voxmask_vol);

        % WM
        vol_WMMask.fname    = fullfile(saveDestination, [saveName '_WM_2ndT1' maskExt2nd]);
        vol_WMMask.descrip  = 'WMmasked_MRS_Voxel_Mask';
        vol_WMMask.dim      = vol_mask_2nd.dim;
        vol_WMMask.dt       = vol_mask_2nd.dt;
        vol_WMMask.mat      = vol_mask_2nd.mat;
        WM_voxmask_vol      = WMvol2nd.private.dat(:,:,:) .* vol_mask_2nd.private.dat(:,:,:);
        vol_WMMask          = spm_write_vol(vol_WMMask, WM_voxmask_vol);

        % CSF
        vol_CSFMask.fname    = fullfile(saveDestination, [saveName '_CSF_2ndT1' maskExt2nd]);
        vol_CSFMask.descrip  = 'CSFmasked_MRS_Voxel_Mask';
        vol_CSFMask.dim      = vol_mask_2nd.dim;
        vol_CSFMask.dt       = vol_mask_2nd.dt;
        vol_CSFMask.mat      = vol_mask_2nd.mat;
        CSF_voxmask_vol     = CSFvol2nd.private.dat(:,:,:) .* vol_mask_2nd.private.dat(:,:,:);
        vol_CSFMask          = spm_write_vol(vol_CSFMask, CSF_voxmask_vol);


        %%% 6. DETERMINE FRACTIONAL TISSUE VOLUMES %%%
        % Sum image intensities over the entire masked tissue specific volume
        GMsum  = sum(sum(sum(vol_GMMask.private.dat(:,:,:))));
        WMsum  = sum(sum(sum(vol_WMMask.private.dat(:,:,:))));
        CSFsum = sum(sum(sum(vol_CSFMask.private.dat(:,:,:))));

        % Normalize
        fGM  = GMsum / (GMsum + WMsum + CSFsum);
        fWM  = WMsum / (GMsum + WMsum + CSFsum);
        fCSF = CSFsum / (GMsum + WMsum + CSFsum);

        % Save normalized fractional tissue volumes to MRSCont
        MRSCont.seg.tissue.secondT1.fGM(kk)  = fGM;
        MRSCont.seg.tissue.secondT1.fWM(kk)  = fWM;
        MRSCont.seg.tissue.secondT1.fCSF(kk) = fCSF;


        %%% 7. DETERMINE PET IMAGE INTENSITY %%%
        % Get the histogram of PET image intensities inside the voxel
        % Get the input file name and the PET image
        vol_image_pet = MRSCont.coreg.vol_image_pet{kk};
        [path,filename,~]   = fileparts(MRSCont.files{kk});
        % For batch analysis, get the last two sub-folders (e.g. site and
        % subject)
        path_split          = regexp(path,filesep,'split');
        if length(path_split) > 2
            saveName = [path_split{end-1} '_' path_split{end} '_' filename];
        end

        % Loop over tissue types (for PET, do WM and GM only)
        tissueTypes = {'GM', 'WM'};
        tissueMasks = {vol_GMMask, vol_WMMask};
        for rr = 1:length(tissueTypes)
            [tissueMaskDir, tissueMaskName, tissueMaskExt] = fileparts(tissueMasks{rr}.fname);
            vol_PETMask.fname    = fullfile(saveDestination, [saveName '_PET_' tissueTypes{rr} tissueMaskExt]);
            vol_PETMask.descrip  = ['PETmasked_', tissueTypes{rr}, '_MRS_Voxel_Mask'];
            vol_PETMask.dim      = tissueMasks{rr}.dim;
            vol_PETMask.dt       = tissueMasks{rr}.dt;
            vol_PETMask.mat      = tissueMasks{rr}.mat;
            PET_voxmask_vol      = vol_image_pet.private.dat(:,:,:) .* tissueMasks{rr}.private.dat(:,:,:);
            vol_PETMask          = spm_write_vol(vol_PETMask, PET_voxmask_vol);

            % Get everything greater than zero
            PETIntensityInVoxel     = nonzeros(vol_PETMask.private.dat(:,:,:));

            % 1st metric: Raw pixel intensity sum over all voxels
            rawPETIntensitySum      = sum(PETIntensityInVoxel);

            % 2nd metric: Create a histogram and determine the maximum of
            % the pixel intensity distribution
            idx = PETIntensityInVoxel > 0.2; % cut off on lower end of greyscale range because of many dark pixels
            PETIntensityInVoxel     = PETIntensityInVoxel(idx);
            nBins                   = round(sqrt(length(PETIntensityInVoxel)));  % Number of bins
            [yIntensity, xIntensity] = hist(PETIntensityInVoxel,nBins);
            % Initial guesses for a Gaussian fit to the histogram
            [maxIntensityInit, maxIndexInit] = max(yIntensity);
            % Get starting values for linear baseline and offset
            grad_points = (yIntensity(end) - yIntensity(1)) ./ abs(xIntensity(end) - xIntensity(1));
            LinearInit  = grad_points ./ abs(xIntensity(1) - xIntensity(2));
            constInit   = (yIntensity(end)+yIntensity(1))/2;
            % Set initials and run model
            GaussModelInit  = [maxIntensityInit -0.2 xIntensity(maxIndexInit) -LinearInit constInit];
            nlinopts        = statset('nlinfit');
            nlinopts        = statset(nlinopts,'MaxIter',400,'TolX',1e-6,'TolFun',1e-6,'FunValCheck','off');
            [GaussModelParams, residPlot] = nlinfit(xIntensity, yIntensity, @GaussModel, GaussModelInit, nlinopts); % re-run for residuals for output figure

            % Save metrics
            MRSCont.seg.pet.rawPETIntensitySum.(tissueTypes{rr})(kk)                = rawPETIntensitySum;
            MRSCont.seg.pet.histogram.xIntensity.(tissueTypes{rr}){kk}              = xIntensity;
            MRSCont.seg.pet.histogram.yIntensity.(tissueTypes{rr}){kk}              = yIntensity;
            MRSCont.seg.pet.histogram.fitParams.(tissueTypes{rr}){kk}               = GaussModelParams;
            MRSCont.seg.pet.histogram.mostFrequentIntensity.(tissueTypes{rr})(kk)  = GaussModelParams(3);
            MRSCont.seg.pet.histogram.distFWHM.(tissueTypes{rr})(kk)                = abs(GaussModelParams(2));
        end

    end
end


fprintf('... done.\n');
time = toc(refSegTime);
[~] = printLog('done',time,1,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);
MRSCont.runtime.Seg = time;
fclose(fileID); %close log file


%% Create table and csv file
if MRSCont.flags.hasSecondT1
    tissueTypes         = {'fGM','fWM','fCSF', 'fGM_secondT1', 'fWM_secondT1', 'fCSF_secondT1'};
    tissue              = horzcat(MRSCont.seg.tissue.fGM',MRSCont.seg.tissue.fWM',MRSCont.seg.tissue.fCSF', MRSCont.seg.tissue.secondT1.fGM', MRSCont.seg.tissue.secondT1.fWM', MRSCont.seg.tissue.secondT1.fCSF');
else
    tissueTypes         = {'fGM','fWM','fCSF'};
    tissue              = horzcat(MRSCont.seg.tissue.fGM',MRSCont.seg.tissue.fWM',MRSCont.seg.tissue.fCSF');
end
MRSCont.seg.tables  = array2table(tissue,'VariableNames',tissueTypes);
writetable(MRSCont.seg.tables,[saveDestination  filesep 'TissueFractions.csv']);

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
    osp_WriteBIDsTable(MRSCont.seg.(['tables_Voxel_' num2str(rr)]), [saveDestination  filesep 'TissueFractions_Voxel_' num2str(rr)])
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
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

spm_jobman('run',matlabbatch);

end
