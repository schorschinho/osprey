% osp_extract_three_plane_image_overlay.m
% Helge Zoellner, Johns Hopkins University 2021.
%
% USAGE:
% [three_plane_image] = osp_extract_three_plane_image(ImageFile, MaskFile);
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


function [three_plane_img, three_plane_overlay,size_vox_t] = osp_extract_three_plane_image_overlay(ImageFile, MaskFile,voxel_ctr,T1_max);
%%% 1. LOAD IMAGE AND MASK VOLUME %%%

    Vimage=spm_vol(ImageFile);
    Vmask=spm_vol(MaskFile);   

    [mask,~] = spm_read_vols(Vmask);
    if isempty(voxel_ctr)
        center = regionprops3(logical(mask),'Centroid');
        center = round(center{1,:});
        center = [center , 1];
        centerWorldSpace = Vmask.mat*center';
    else
        centerWorldSpace=voxel_ctr;
    end
    

    %%% 2. SET UP THREE PLANE IMAGE %%%
    % Generate three plane image for the output
    % Transform structural image and co-registered voxel mask from voxel to
    % world space for output (MM: 180221)

    [img_t,img_c,img_s] = voxel2world_space(Vimage,centerWorldSpace(1:3));
    [mask_t,mask_c,mask_s] = voxel2world_space(Vmask,centerWorldSpace(1:3));

    img_t = flipud(img_t/T1_max);
    img_c = flipud(img_c/T1_max);
    img_s = flipud(img_s/T1_max);


    size_max = max([max(size(img_t)) max(size(img_c)) max(size(img_s))]);
    three_plane_img = zeros([size_max 3*size_max]);
    three_plane_img(:,1:size_max)              = image_center(img_t, size_max);
    three_plane_img(:,size_max+(1:size_max))   = image_center(img_s, size_max);
    three_plane_img(:,size_max*2+(1:size_max)) = image_center(img_c, size_max);


    mask_t = flipud(mask_t);
    mask_c = flipud(mask_c);
    mask_s = flipud(mask_s);

%     mask_t(mask_t==0) = nan;
%     mask_c(mask_c==0) = nan;
%     mask_s(mask_s==0) = nan;




    size_max = max([max(size(mask_t)) max(size(mask_c)) max(size(mask_s))]);
    three_plane_overlay = zeros([size_max 3*size_max]);
    three_plane_overlay(:,1:size_max)              = image_center(mask_t, size_max);
    three_plane_overlay(:,size_max+(1:size_max))   = image_center(mask_s, size_max);
    three_plane_overlay(:,size_max*2+(1:size_max)) = image_center(mask_c, size_max);
    size_vox_t = size(mask_t,2);
end