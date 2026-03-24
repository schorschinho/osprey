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
osp_CheckRunPreviousModule(MRSCont, 'OspreySeg');
[~,MRSCont.ver.CheckOsp ] = osp_Toolbox_Check ('OspreySeg',MRSCont.flags.isGUI);

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
for kk = 1:MRSCont.nDatasets  
     [~] = printLog('OspreySeg',kk,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
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
        [T1dir, T1name, T1ext]  = fileparts(niftiFile);
        if strcmp(T1ext,'.gz')
            T1name = strrep(T1name, '.nii','');
        end
        

        segFile               = fullfile(T1dir, [T1name '_seg8.mat']);
        % If a GM-segmented file doesn't exist, start the segmentation
        if ~exist(segFile,'file')
            %Uncompress .nii.gz if needed
            if strcmp(T1ext,'.gz')
                gunzip(niftiFile);
                niftiFile = strrep(niftiFile,'.gz','');
                T1ext = '.nii';
            end           
            createSegJob(niftiFile,MRSCont.flags.isMRSI);
        else
            if strcmp(T1ext,'.gz')
                gunzip(niftiFile);
                niftiFile = strrep(niftiFile,'.gz','');
                T1ext = '.nii';
            end  
            if exist(fullfile(T1dir, ['c1' T1name '.nii.gz']),'file')
                gunzip(fullfile(T1dir, ['c1' T1name T1ext '.gz']));
                gunzip(fullfile(T1dir, ['c2' T1name T1ext '.gz']));
                gunzip(fullfile(T1dir, ['c3' T1name T1ext '.gz']));  
                if MRSCont.flags.isMRSI
                    gunzip(fullfile(T1dir, ['c4' T1name T1ext '.gz'])); 
                    gunzip(fullfile(T1dir, ['iy_' T1name T1ext '.gz'])); 
                    gunzip(fullfile(T1dir, ['y_' T1name T1ext '.gz'])); 
                end
            end
            T1ext = strrep(T1ext,'.gz','');
        end


        %%% 2. CREATE MASKED TISSUE MAPS %%%
        % Define file names
        segFileGM   = fullfile(T1dir, ['c1' T1name T1ext]);
        segFileWM   = fullfile(T1dir, ['c2' T1name T1ext]);
        segFileCSF  = fullfile(T1dir, ['c3' T1name T1ext]);
        if MRSCont.flags.isMRSI
            segFileLIP  = fullfile(T1dir, ['c4' T1name T1ext]);
            segFileMRSI_iy  = fullfile(T1dir, ['iy_' T1name T1ext]);
            segFileMRSI_y  = fullfile(T1dir, ['y_' T1name T1ext]);
        end
        % Load volumes
        GMvol  = spm_vol(segFileGM);
        WMvol  = spm_vol(segFileWM);
        CSFvol = spm_vol(segFileCSF);
        if MRSCont.flags.isMRSI
            LIPvol = spm_vol(segFileLIP);
            T1vol = spm_vol(fullfile(T1dir, [T1name T1ext]));
        end
        
        %Loop over voxels (for DualVoxel)
        if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 0))
            Voxels = 1;
        else
            if ~MRSCont.flags.isMRSI
                Voxels = 2;
            else
                Voxels = MRSCont.raw{kk}.nZvoxels;
            end
        end
        if MRSCont.flags.isMRSI
            vx =1; 
            total_voxels = MRSCont.raw{kk}.nXvoxels * MRSCont.raw{kk}.nYvoxels *MRSCont.raw{kk}.nZvoxels;
            if ~exist(MRSCont.coreg.index_mask{kk}.fname,'file')
                gunzip([MRSCont.coreg.index_mask{kk}.fname,'.gz']);
            end
            index_mask = nii_tool('load', MRSCont.coreg.index_mask{kk}.fname);
            index_mask = double(index_mask.img);

            index_mask = round(index_mask);
            if ~isempty(MRSCont.coreg.gap_mask{kk})
                if ~exist(MRSCont.coreg.gap_mask{kk}.fname,'file')
                    gunzip([MRSCont.coreg.gap_mask{kk}.fname,'.gz']);
                end
                gap_mask = nii_tool('load', MRSCont.coreg.gap_mask{kk}.fname);
                gap_mask = double(gap_mask.img);
          
            end
            % Create and save masked tissue maps
            % Get the input file name
            [path,filename,~]   = fileparts(MRSCont.files{kk});
            % For batch analysis, get the last two sub-folders (e.g. site and
            % subject)
            path_split          = regexp(path,filesep,'split');
            if length(path_split) > 2
                saveName = [path_split{end-1} '_' path_split{end} '_' filename];
            end

            % SPM load
            GMvols  = spm_vol(segFileGM);
            WMvol  = spm_vol(segFileWM);
            CSFvol = spm_vol(segFileCSF);                     
            LIPvol = spm_vol(segFileLIP);
            GMvol  = GMvol.private.dat(:,:,:);
            WMvol  = WMvol.private.dat(:,:,:);
            CSFvol = CSFvol.private.dat(:,:,:);                
            LIPvol = LIPvol.private.dat(:,:,:);
            MRSCont.seg.tissue.fGM(kk,:,:,:)  = zeros(MRSCont.raw{kk}.nXvoxels,MRSCont.raw{kk}.nYvoxels,MRSCont.raw{kk}.nZvoxels);
            MRSCont.seg.tissue.fWM(kk,:,:,:)  = zeros(MRSCont.raw{kk}.nXvoxels,MRSCont.raw{kk}.nYvoxels,MRSCont.raw{kk}.nZvoxels);
            MRSCont.seg.tissue.fCSF(kk,:,:,:)  = zeros(MRSCont.raw{kk}.nXvoxels,MRSCont.raw{kk}.nYvoxels,MRSCont.raw{kk}.nZvoxels);
            MRSCont.seg.tissue.fLIP(kk,:,:,:)  = zeros(MRSCont.raw{kk}.nXvoxels,MRSCont.raw{kk}.nYvoxels,MRSCont.raw{kk}.nZvoxels);

            if isfield(MRSCont.opts.MRSI,'outerMask')
                if ~isfield(MRSCont.opts.MRSI.outerMask,'mask')  
                    MRSCont.opts.MRSI.outerMask.mask = zeros(MRSCont.raw{kk}.nXvoxels,MRSCont.raw{kk}.nYvoxels,MRSCont.raw{kk}.nZvoxels);
                    MRSCont.opts.MRSI.outerMask.mask(MRSCont.opts.MRSI.outerMask.x(1):MRSCont.opts.MRSI.outerMask.x(2),...
                                                     MRSCont.opts.MRSI.outerMask.y(1):MRSCont.opts.MRSI.outerMask.y(2),...
                                                     MRSCont.opts.MRSI.outerMask.z(1):MRSCont.opts.MRSI.outerMask.z(2)) = 1;
                    % MRSCont.opts.MRSI.outerMask.mask = squeeze(MRSCont.opts.MRSI.outerMask.mask);
                end
            end

            % Pre-compute gap mask if needed
            if MRSCont.opts.MRSI.pseudo3D && ~isempty(MRSCont.coreg.gap_mask{kk})
                gap_mask_logical = (gap_mask == 1);
            end

        end
        for rr = 1 : Voxels % This also loops over the different slice for MRSI
            if ~MRSCont.flags.isMRSI
            % Get voxel mask filename
                if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
                    vol_mask = MRSCont.coreg.vol_mask{kk};
                else
                    vol_mask = MRSCont.coreg.vol_mask{kk}{rr};
                end
                if ~exist(vol_mask.fname,'file') && exist(strrep(vol_mask.fname,'.nii','.nii.gz'),'file')
                    gunzip(strrep(vol_mask.fname,'.nii','.nii.gz'));
                end
                [~, ~, maskExt] = fileparts(vol_mask.fname);
    
                % Create and save masked tissue maps
                % Get the input file name
                [path,filename,~]   = fileparts(MRSCont.files{kk});
                % For batch analysis, get the last two sub-folders (e.g. site and
                % subject)
                path_split          = regexp(path,filesep,'split');
                if length(path_split) > 2
                    saveName = [path_split{end-1} '_' path_split{end} '_' filename];
                end
                
                %Add voxel number for DualVoxel
                if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
                    VoxelNum = '_Voxel_1';
                else
                    VoxelNum = ['_Voxel_' num2str(rr)];
                end
                
                % GM
                vol_GMMask.fname    = fullfile(saveDestination, [saveName VoxelNum '_GM' maskExt]);
                vol_GMMask.descrip  = ['GMmasked_MRS_Voxel_Mask_' VoxelNum];
                vol_GMMask.dim      = vol_mask.dim;
                vol_GMMask.dt       = vol_mask.dt;
                vol_GMMask.mat      = vol_mask.mat;
                GM_voxmask_vol      = GMvol.private.dat(:,:,:) .* vol_mask.private.dat(:,:,:);
                vol_GMMask          = spm_write_vol(vol_GMMask, GM_voxmask_vol);
    
                % WM
                vol_WMMask.fname    = fullfile(saveDestination, [saveName VoxelNum '_WM' maskExt]);
                vol_WMMask.descrip  = ['WMmasked_MRS_Voxel_Mask_' VoxelNum];
                vol_WMMask.dim      = vol_mask.dim;
                vol_WMMask.dt       = vol_mask.dt;
                vol_WMMask.mat      = vol_mask.mat;
                WM_voxmask_vol      = WMvol.private.dat(:,:,:) .* vol_mask.private.dat(:,:,:);
                vol_WMMask          = spm_write_vol(vol_WMMask, WM_voxmask_vol);
    
                % CSF
                vol_CSFMask.fname   = fullfile(saveDestination, [saveName VoxelNum '_CSF' maskExt]);
                vol_CSFMask.descrip = ['CSFmasked_MRS_Voxel_Mask_' VoxelNum];
                vol_CSFMask.dim     = vol_mask.dim;
                vol_CSFMask.dt      = vol_mask.dt;
                vol_CSFMask.mat     = vol_mask.mat;
                CSF_voxmask_vol     = CSFvol.private.dat(:,:,:) .* vol_mask.private.dat(:,:,:);
                vol_CSFMask         = spm_write_vol(vol_CSFMask, CSF_voxmask_vol);
            end

            % For MRSI data
            if MRSCont.flags.isMRSI
            

                for y = 1 : MRSCont.raw{kk}.nYvoxels
                    for x = 1 : MRSCont.raw{kk}.nXvoxels
                        if vx == 1    
                            msg = sprintf('Calculating voxelfraction from voxel %3i out of %3i total voxels...\n', vx, total_voxels);
                            fprintf(msg);
                         else
                             msg = sprintf('Calculating voxelfraction from voxel %3i out of %3i total voxels...\n', vx, total_voxels);
                                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                                fprintf([reverseStr, msg]);
                        end
                        if MRSCont.opts.MRSI.outerMask.mask(x,y,rr)
                            index = x * 1e6 + y * 1e3 + rr;
                            index_mask_temp =zeros(size(index_mask));
                            index_mask_temp(index_mask==index) =1;
                            
                            % Apply gap mask if needed
                            if MRSCont.opts.MRSI.pseudo3D  && ~isempty(MRSCont.coreg.gap_mask{kk})
                                index_mask_temp = index_mask_temp & ~gap_mask_logical;
                            end

                            mask_indices_temp = find(index_mask_temp);
                          
                           if ~isempty(mask_indices_temp)
                                % Extract values only at mask locations
                                GM_vals = GMvol(mask_indices_temp);
                                WM_vals = WMvol(mask_indices_temp);
                                CSF_vals = CSFvol(mask_indices_temp);
                                LIP_vals = LIPvol(mask_indices_temp);
                                
                                % Sum values above threshold
                                GMsum = sum(GM_vals(GM_vals > 0.9));
                                WMsum = sum(WM_vals(WM_vals > 0.9));
                                CSFsum = sum(CSF_vals(CSF_vals > 0.9));
                                LIPsum = sum(LIP_vals(LIP_vals > 0.9));
                                
                                % Calculate fractions
                                total_tissue = GMsum + WMsum + CSFsum;
                                if total_tissue == 0
                                    fGM = 0;
                                    fWM = 0;
                                    fCSF = 0;
                                else
                                    fGM = GMsum / total_tissue;
                                    fWM = WMsum / total_tissue;
                                    fCSF = CSFsum / total_tissue;
                                end
                                
                                if LIPsum == 0
                                    fLIP = 0;
                                else
                                    fLIP = LIPsum / length(mask_indices_temp);
                                end
                                
                                if fLIP < 0.01
                                    fLIP = 0;
                                end

                                MRSCont.seg.tissue.fGM(kk,x,y,rr) = fGM;
                                MRSCont.seg.tissue.fWM(kk,x,y,rr) = fWM;
                                MRSCont.seg.tissue.fCSF(kk,x,y,rr) = fCSF;
                                MRSCont.seg.tissue.fLIP(kk,x,y,rr) = fLIP;
                            end
                        end
                        vx = vx + 1;                 
                    end
                end
            end
        end
    
                if ~isfield(MRSCont.opts.MRSI,'threshBrain')
                    MRSCont.opts.MRSI.threshBrain = 0.3;
                end
                if ~isfield(MRSCont.opts.MRSI,'threshLipid')
                    MRSCont.opts.MRSI.threshLipid = 0.1;
                end

                temp = squeeze(MRSCont.seg.tissue.fGM(kk,:,:,:) + MRSCont.seg.tissue.fWM(kk,:,:,:));
                temp(temp>MRSCont.opts.MRSI.threshBrain)= 1;
                temp(temp<1)= 0;
                temp = spm_dilate(temp);
                temp = spm_erode(temp);
                MRSCont.seg.tissue.brain(kk,:,:,:) = temp;
            
                temp = MRSCont.seg.tissue.fLIP(kk,:,:,:);
                temp(temp>MRSCont.opts.MRSI.threshLipid)= 1;
                temp(temp<1)= 0;
                MRSCont.seg.tissue.lip(kk,:,:,:) = temp;
            
                in = MRSCont.raw{kk};
                if (in.nZvoxels > 1) && ~MRSCont.opts.MRSI.pseudo3D 
                    reorder = flip(1:in.nZvoxels);
                    for ll = 1 : in.nZvoxels
                            ToExport = in;
                            if isfield(ToExport.geometry,'slice_distance')
                                VoxelShift = [0 , 0, ToExport.geometry.slice_distance/ToExport.geometry.size.cc*shift];
                                ToExport = osp_shift_nii_volume(ToExport,VoxelShift); % Update slice 
                            end
                            out.hdr = ToExport.nii_mrs.hdr;
                            % out.hdr.dim(1) = 3;
                            % out.hdr.pixdim(1) = -1;
                            out.hdr.dim(1) = 3;
                            out.hdr.dim(2) = 1;
                            out.hdr.pixdim(5) = 1;
                            out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes  
                            out.img = squeeze(MRSCont.seg.tissue.fGM(kk,:,:,ll));
                            nii_tool('save', out, fullfile(saveDestination,['fGM_slice_' num2str(reorder(ll)) '.nii.gz']));  
                            out.img = squeeze(MRSCont.seg.tissue.fWM(kk,:,:,ll));
                            nii_tool('save', out, fullfile(saveDestination,['fWM_slice_' num2str(reorder(ll)) '.nii.gz']));  
                            out.img = squeeze(MRSCont.seg.tissue.fCSF(kk,:,:,ll));
                            nii_tool('save', out, fullfile(saveDestination,['fCSF_slice_' num2str(reorder(ll)) '.nii.gz']));  
                            out.img = squeeze(MRSCont.seg.tissue.brain(kk,:,:,ll));
                            nii_tool('save', out, fullfile(saveDestination,['brain_slice_' num2str(reorder(ll)) '.nii.gz']));  
                            out.img = squeeze(MRSCont.seg.tissue.fLIP(kk,:,:,ll));
                            nii_tool('save', out, fullfile(saveDestination,['fLIP_slice_' num2str(reorder(ll)) '.nii.gz']));
                            out.img = squeeze(MRSCont.seg.tissue.lip(kk,:,:,ll));
                            nii_tool('save', out, fullfile(saveDestination,['lip_slice_' num2str(reorder(ll)) '.nii.gz']));
            
                            shift = shift - 1;
                    end
                else
                    ToExport = in;
                    out.hdr = ToExport.nii_mrs.hdr;
                    % out.hdr.dim(1) = 3;
                    % out.hdr.pixdim(1) = -1;
                    out.hdr.dim(1) = 3;
                    out.hdr.dim(2) = 1;
                    out.hdr.pixdim(5) = 1;
                    out.hdr.pixdim(1) = 1; % For the other dataset this was -1 QForm changes  
                    out.img = squeeze(MRSCont.seg.tissue.fGM(kk,:,:,:));
                    nii_tool('save', out, fullfile(saveDestination,'fGM.nii.gz'));  
                    out.img = squeeze(MRSCont.seg.tissue.fWM(kk,:,:,:));
                    nii_tool('save', out, fullfile(saveDestination,'fWM.nii.gz'));  
                    out.img = squeeze(MRSCont.seg.tissue.fCSF(kk,:,:,:));
                    nii_tool('save', out, fullfile(saveDestination,'fCSF.nii.gz'));  
                    out.img = squeeze(MRSCont.seg.tissue.brain(kk,:,:,:));
                    nii_tool('save', out, fullfile(saveDestination,'brain.nii.gz'));
                    out.img = squeeze(MRSCont.seg.tissue.fLIP(kk,:,:,:));
                    nii_tool('save', out, fullfile(saveDestination,'fLIP.nii.gz'));
                    out.img = squeeze(MRSCont.seg.tissue.lip(kk,:,:,:));
                    nii_tool('save', out, fullfile(saveDestination,'lip.nii.gz'));
                end
            
                % Let's create the brain surface as gii file
                spm_imcalc({segFileGM; segFileWM; segFileCSF}, fullfile(T1dir, [T1name '_brain' T1ext]), '(i1+i2+i3)>0.1');
                spm_surf(fullfile(T1dir, [T1name '_brain' T1ext]), 2, 0.5);
                MRSCont.gii_filename_brain{kk} = fullfile(T1dir, [T1name '_brain.gii']);

    
            end
            %%% 3. DETERMINE FRACTIONAL TISSUE VOLUMES %%%
            % Sum image intensities over the entire masked tissue specific volume
            if ~MRSCont.flags.isMRSI
                GMsum  = sum(sum(sum(vol_GMMask.private.dat(:,:,:))));
                WMsum  = sum(sum(sum(vol_WMMask.private.dat(:,:,:))));
                CSFsum = sum(sum(sum(vol_CSFMask.private.dat(:,:,:))));
            end
            
            % Save three plane image to container
            if MRSCont.flags.addImages                
                [MRSCont.seg.img_montage{kk},MRSCont.seg.size_vox_t(kk)] = osp_extract_three_plane_image_seg(niftiFile, vol_mask,vol_GMMask,vol_WMMask,vol_CSFMask,MRSCont.coreg.voxel_ctr{kk},MRSCont.coreg.T1_max{kk});
            end
            
            %Compress nifit and delete uncompressed files
            if ~MRSCont.flags.isMRSI
                gzip(vol_GMMask.fname);
                delete(vol_GMMask.fname);
                gzip(vol_WMMask.fname);
                delete(vol_WMMask.fname);
                gzip(vol_CSFMask.fname);
                delete(vol_CSFMask.fname);
            
                gzip(GMvol.fname);
                delete(GMvol.fname);
                gzip(WMvol.fname);
                delete(WMvol.fname);
                gzip(CSFvol.fname);
                delete(CSFvol.fname);
                delete(vol_mask.fname);
                
                try
                    gzip(MRSCont.coreg.vol_image{kk}.fname)
                    delete(MRSCont.coreg.vol_image{kk}.fname);
                catch
                end
            else
                gzip(segFileGM);
                delete(segFileGM);
                gzip(segFileWM);
                delete(segFileWM);
                gzip(segFileCSF);
                delete(segFileCSF);
                gzip(segFileLIP);
                delete(segFileLIP);
                gzip(fullfile(T1dir, [T1name '_brain' T1ext]));
                delete(fullfile(T1dir, [T1name '_brain' T1ext]));
                gzip(segFileMRSI_iy);
                delete(segFileMRSI_iy);
                gzip(segFileMRSI_y);
                delete(segFileMRSI_y);                            
            end



            if ~MRSCont.flags.isMRSI
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

time = toc(refSegTime);
[~] = printLog('done',time,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
MRSCont.runtime.Seg = time;
%% Create table and csv file or NIfTI files
tissueTypes = {'fGM','fWM','fCSF'};
%Loop over voxels (for DualVoxel)

if ~MRSCont.flags.isMRSI
    for rr = 1 : Voxels
        tissue = horzcat(MRSCont.seg.tissue.fGM(:,rr),MRSCont.seg.tissue.fWM(:,rr),MRSCont.seg.tissue.fCSF(:,rr));
        MRSCont.seg.(['tables_Voxel_' num2str(rr)]) = array2table(tissue,'VariableNames',tissueTypes);
        writetable(MRSCont.seg.(['tables_Voxel_' num2str(rr)]),[saveDestination  filesep 'TissueFractions_Voxel_' num2str(rr) '.csv']);
    end
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


function createSegJob(T1file,isMRSI)

% Created with SPM12 batch manager (standard options)
spmhome = fileparts(which('spm'));
tpm = cellstr(spm_select('ExtFPList',fullfile(spmhome,'tpm'),'TPM.nii'));
% SPM can not handle hidden files
[~,names,~] = cellfun(@fileparts, tpm, 'UniformOutput', false); %#ok<*STRCLFH>
hidden = logical(ones(1,length(tpm)));
for jj = 1:length(tpm) 
    if ~strcmp(names{jj}(1),'.')
        hidden(jj) = 0;
    end
end
tpm = tpm(~hidden);%delete hidden files 

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
if ~isMRSI
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
else
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
end
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
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];

spm_jobman('run',matlabbatch);

end



