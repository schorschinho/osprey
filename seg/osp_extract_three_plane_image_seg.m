% osp_extract_image_montage.m
% Helge Zoellner, Johns Hopkins University 2021.
%
% USAGE:
% [three_plane_image] = osp_extract_image_montage(ImageFile, MaskFile);
%
% DESCRIPTION:
% Creates a SPM volume containing a voxel mask with the same dimensions
% as the SPM volume containing a structural image. The voxel mask will be
% created with the geometry information stored in the FID-A structure "in".
% This routine works for the GE P data type.

% INPUTS:
% ImageFile   = Path to Image File (NifTI).
% MaskFile    = Path to Mask Volume (NifTI).
% VoxelIndex  = Voxel index for PRIAM
%
% OUTPUTS:
% three_plane_img  = three plane image for plots and GUI.


function [img_montage,size_vox_t] = osp_extract_three_plane_image_seg(ImageFile, MaskFile,GM,WM,CSF,voxel_ctr,T1_max);
%%% 1. LOAD IMAGE AND MASK VOLUME %%%

    Vimage=spm_vol(ImageFile);
    Vmask=spm_vol(MaskFile);    
    
    vol_GM_mask  = spm_vol(GM);
    vol_WM_mask  = spm_vol(WM);
    vol_CSF_mask = spm_vol(CSF);

    %%% 2. SET UP THREE PLANE IMAGE %%%
    % Generate three plane image for the output
    % Transform structural image and co-registered voxel mask from voxel to
    % world space for output (MM: 180221)

    img_t     = flipud(voxel2world_space(Vimage, voxel_ctr));
    vox_t     = flipud(voxel2world_space(Vmask, voxel_ctr));
    vox_t_GM  = flipud(voxel2world_space(vol_GM_mask, voxel_ctr));
    vox_t_WM  = flipud(voxel2world_space(vol_WM_mask, voxel_ctr));
    vox_t_CSF = flipud(voxel2world_space(vol_CSF_mask, voxel_ctr));
    img_t = img_t/T1_max;
    img_montage = [img_t+0.225*vox_t, img_t+0.3*vox_t_GM, img_t+0.225*vox_t_WM, img_t+0.4*vox_t_CSF];
    size_vox_t = size(vox_t,2);
end