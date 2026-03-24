function [MRSCont] = OspreyAddImages(MRSCont)
%% [MRSCont] = OspreyAddImages(MRSCont)
%   This function allows adds the three plane images into the container to
%   make it transferable
%
%   USAGE:
%       MRSCont = OspreyAddIamges(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2021-05-06)
%       hzoelln2@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2021-07-23: First version of the code.


%% Pack images from coregistration
osp_CheckRunPreviousModule(MRSCont, 'OspreyCoreg')
for kk = 1 : MRSCont.nDatasets
    % Load T1 image, mask volume, T1 max value, and voxel center

        if ~exist(MRSCont.coreg.vol_image{kk}.fname,'file')
            gunzip([MRSCont.coreg.vol_image{kk}.fname, '.gz']);
        end
        if ~exist(MRSCont.coreg.vol_mask{kk}.fname,'file')
            gunzip([MRSCont.coreg.vol_mask{kk}.fname, '.gz']);
        end

        %%% 3. SET UP THREE PLANE IMAGE %%%
        if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
            [MRSCont.coreg.three_plane_img{kk}] = osp_extract_three_plane_image(MRSCont.coreg.vol_image{kk}.fname, MRSCont.coreg.vol_mask{kk}.fname,MRSCont.coreg.voxel_ctr{kk},MRSCont.coreg.T1_max{kk});
        else
            [MRSCont.coreg.three_plane_img{kk}] = osp_extract_three_plane_image(MRSCont.coreg.vol_image{kk}.fname, MRSCont.coreg.vol_mask{kk}{VoxelIndex}.fname,MRSCont.coreg.voxel_ctr{kk}(:,:,VoxelIndex),MRSCont.coreg.T1_max{kk});
        end

        if ~MRSCont.flags.didSeg
            if exist([MRSCont.coreg.vol_mask{kk}.fname, '.gz'],'file')
                delete(MRSCont.coreg.vol_mask{kk}.fname);
            end
            if exist([MRSCont.coreg.vol_image{kk}.fname, '.gz'],'file')
                delete(MRSCont.coreg.vol_image{kk}.fname);
            end
        end

        %%% 4. SET UP FIGURE LAYOUT %%%
        % Generate a new figure and keep the handle memorized

        if ~MRSCont.flags.didSeg
            if exist([MRSCont.coreg.vol_mask{kk}.fname, '.gz'],'file')
                delete(MRSCont.coreg.vol_mask{kk}.fname);
            end
            if exist([MRSCont.coreg.vol_image{kk}.fname, '.gz'],'file')
                delete(MRSCont.coreg.vol_image{kk}.fname);
            end
        end

end

%% Pack images from Segmentation
osp_CheckRunPreviousModule(MRSCont, 'OspreySeg')
for kk = 1 : MRSCont.nDatasets
    [path_voxel,filename_voxel,~]   = fileparts(MRSCont.files{kk});

    % For batch analysis, get the last two sub-folders (e.g. site and
    % subject)
    path_split          = regexp(path_voxel,filesep,'split');
    if length(path_split) > 2
        saveName = [path_split{end-1} '_' path_split{end} '_' filename_voxel];
    end

    segDestination = fullfile(MRSCont.outputFolder, 'SegMaps');
    
    %Loop over voxels (for DualVoxel)
    if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
        Voxels = 1;
    else
        Voxels = 2;
    end
    for rr = 1 : Voxels

        if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
            VoxelNum = '_Voxel_1';
        else
            VoxelNum = ['_Voxel_' num2str(rr)];
        end

        GM  = fullfile(segDestination, [saveName VoxelNum '_GM.nii']);
        WM  = fullfile(segDestination, [saveName VoxelNum '_WM.nii']);
        CSF = fullfile(segDestination, [saveName VoxelNum '_CSF.nii']);

        if exist([GM '.gz'], 'file')
            gunzip([GM '.gz']);
            gunzip([WM '.gz']);
            gunzip([CSF '.gz']);
        else
            if length(path_split) > 2
                saveName = ['jobServer_' path_split{end} '_' filename_voxel];
            end

            segDestination = fullfile(MRSCont.outputFolder, 'SegMaps');

            if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
                VoxelNum = '_Voxel_1';
            else
                VoxelNum = ['_Voxel_' num2str(rr)];
            end

            GM  = fullfile(segDestination, [saveName VoxelNum '_GM.nii']);
            WM  = fullfile(segDestination, [saveName VoxelNum '_WM.nii']);
            CSF = fullfile(segDestination, [saveName VoxelNum '_CSF.nii']);
            if exist([GM '.gz'], 'file')
                gunzip([GM '.gz']);
                gunzip([WM '.gz']);
                gunzip([CSF '.gz']);
            end
        end

        if ~exist(GM, 'file') % This is just for combability of older Osprey versions
            GM  = fullfile(segDestination, [saveName '_GM.nii']);
            WM  = fullfile(segDestination, [saveName '_WM.nii']);
            CSF = fullfile(segDestination, [saveName '_CSF.nii']);
        end

        vol_GM_mask  = spm_vol(GM);
        vol_WM_mask  = spm_vol(WM);
        vol_CSF_mask = spm_vol(CSF);

         [~, ~, T1ext]  = fileparts(MRSCont.coreg.vol_image{kk}.fname);
        if strcmp(T1ext, '.gz') && ~exist(MRSCont.coreg.vol_image{kk}.fname,'file') 
            gunzip(MRSCont.coreg.vol_image{kk}.fname)
        end

        %%% 3. SET UP THREE PLANE IMAGE %%%
        if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
            [MRSCont.seg.img_montage{kk},MRSCont.seg.size_vox_t(kk)] = osp_extract_three_plane_image_seg(MRSCont.coreg.vol_image{kk}.fname, MRSCont.coreg.vol_mask{kk}.fname,vol_GM_mask,vol_WM_mask,vol_CSF_mask,MRSCont.coreg.voxel_ctr{kk},MRSCont.coreg.T1_max{kk});
        else
            [MRSCont.seg.img_montage{kk},MRSCont.seg.size_vox_t(kk)] = osp_extract_three_plane_image_seg(MRSCont.coreg.vol_image{kk}.fname, MRSCont.coreg.vol_mask{kk}{rr}.fname,vol_GM_mask,vol_WM_mask,vol_CSF_mask,MRSCont.coreg.voxel_ctr{kk}(:,:,VoxelIndex),MRSCont.coreg.T1_max{kk});
        end
    end

    if exist([GM, '.gz'],'file')
        delete(GM);
        delete(WM);
        delete(CSF);
    end
    if exist([MRSCont.coreg.vol_image{kk}.fname, '.gz'],'file')
        delete(MRSCont.coreg.vol_image{kk}.fname);
    end
    if exist([MRSCont.coreg.vol_mask{kk}.fname, '.gz'],'file')
        delete(MRSCont.coreg.vol_mask{kk}.fname);
    end
end
%% Update flags
    MRSCont.flags.addImages = 1;
end