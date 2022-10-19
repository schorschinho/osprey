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
        if exist(fullfile(T1dir,['r' T1name '_pad.nii']))
            [T1dir, T1name, T1ext]  = fileparts(fullfile(T1dir,['r' T1name '_pad.nii']));
            niftiFile = fullfile(T1dir,[T1name T1ext]);
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
            createSegJob(niftiFile);
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
            end
            T1ext = strrep(T1ext,'.gz','');
        end


        %%% 2. CREATE MASKED TISSUE MAPS %%%
        % Define file names
        segFileGM   = fullfile(T1dir, ['c1' T1name T1ext]);
        segFileWM   = fullfile(T1dir, ['c2' T1name T1ext]);
        segFileCSF  = fullfile(T1dir, ['c3' T1name T1ext]);
        % Load volumes
        GMvol  = spm_vol(segFileGM);
        WMvol  = spm_vol(segFileWM);
        CSFvol = spm_vol(segFileCSF);
        
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
        for rr = 1 : Voxels
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
                [maskDir, maskName, maskExt] = fileparts(vol_mask.fname);
    
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
                vol_mask = MRSCont.coreg.vol_mask{kk}{rr};
                maskFile = vol_mask.fname;
                index_maskFileOut            = vol_mask.fname;  
                index_maskFileOut = strrep(index_maskFileOut,'VoxelMask_','IndexVoxelMask_');
                vol_slice = spm_vol(index_maskFileOut);
                [index_mask,~]    = spm_read_vols(vol_slice);
                index_mask = round(index_mask);
                for y = 1 : MRSCont.raw{kk}.nYvoxels
                    for x = 1 : MRSCont.raw{kk}.nXvoxels
                        index = str2num([sprintf('%02d',x) sprintf('%02d',y) ]);  
                        index_mask_temp =zeros(size(index_mask));
                        index_mask_temp(index_mask==index) =1;
                        GM_voxmask_vol = GMvol.private.dat(:,:,:) .* index_mask_temp;
                        WM_voxmask_vol = WMvol.private.dat(:,:,:) .* index_mask_temp;
                        CSF_voxmask_vol = CSFvol.private.dat(:,:,:) .* index_mask_temp;
                        GMsum  = sum(sum(sum(GM_voxmask_vol(:,:,:))));
                        WMsum  = sum(sum(sum(WM_voxmask_vol(:,:,:))));
                        CSFsum = sum(sum(sum(CSF_voxmask_vol(:,:,:))));
                        if sum([GMsum WMsum CSFsum])==0
                            fGM  = 0;
                            fWM  = 0;
                            fCSF = 0;
                        else
                            fGM  = GMsum / (GMsum + WMsum + CSFsum);
                            fWM  = WMsum / (GMsum + WMsum + CSFsum);
                            fCSF = CSFsum / (GMsum + WMsum + CSFsum);
                        end
                        MRSCont.seg.tissue.fGM(kk,x,y,rr)  = fGM;
                        MRSCont.seg.tissue.fWM(kk,x,y,rr)  = fWM;
                        MRSCont.seg.tissue.fCSF(kk,x,y,rr) = fCSF; 
                        if (fGM + fWM) > .3
                            MRSCont.seg.tissue.brain(kk,x,y,rr) = 1;  
                        else
                            MRSCont.seg.tissue.brain(kk,x,y,rr) = 0;  
                        end
                    end
                end
    
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
    end 
end
time = toc(refSegTime);
[~] = printLog('done',time,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
MRSCont.runtime.Seg = time;
%% Create table and csv file
tissueTypes = {'fGM','fWM','fCSF'};
%Loop over voxels (for DualVoxel)

for rr = 1 : Voxels
    tissue = horzcat(MRSCont.seg.tissue.fGM(:,rr),MRSCont.seg.tissue.fWM(:,rr),MRSCont.seg.tissue.fCSF(:,rr));
    MRSCont.seg.(['tables_Voxel_' num2str(rr)]) = array2table(tissue,'VariableNames',tissueTypes);
    writetable(MRSCont.seg.(['tables_Voxel_' num2str(rr)]),[saveDestination  filesep 'TissueFractions_Voxel_' num2str(rr) '.csv']);
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



