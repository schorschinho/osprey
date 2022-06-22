% coreg_nifti.m
% Chris Davies-Jenkins, Johns Hopkins University 2022.
% 
% USAGE:
% [vol_mask, T1_max, voxel_ctr] = coreg_nifti(in, vol_image, maskFile)
% 
% DESCRIPTION:
% Creates a SPM volume containing a voxel mask with the same dimensions 
% as the SPM volume containing a structural image. The voxel mask is
% created by transforming the NIfTI MRS voxel into T1 resolution, removing
% MRS deminsions, then writing this to maskFile using nii_tool. The 
% resulting file is then loaded using SPM to remove NaN background and 
% edit the header information.
% 
% CREDITS:
% Thanks to Xiangrui Li, Ph.D. for his helpful suggestions using nii_tool.
%
% INPUTS:
% in        = Input data structure.
% vol_image	= SPM volume of the structural image that the voxel defined in
%           the 'geometry' field of the input data structure 'in' is
%           supposed to be co-registered to.
% maskFile  = Filename under which the SPM volume of the co-registered
%           voxel mask is supposed to be saved.
%
% OUTPUTS:
% vol_mask  = SPM volume of the coregistered voxel mask.
% T1_max    = maximum intensity of the image volume.

function [vol_mask, T1_max, voxel_ctr] = coreg_nifti(in, vol_image, maskFile)

image_fname = vol_image.fname;

NiiStruct = nii_tool('load',image_fname);%load structural nifti
NiiVox = nii_tool('load',in.nii_mrs.hdr.file_name);%load voxel nifti

% fh = nii_viewer(NiiStruct,NiiVox);%plot voxel over structural

% Assume Nii voxel and structural are in same space:
NiiVox.hdr.sform_code = NiiStruct.hdr.sform_code;
NiiVox.hdr.qform_code = NiiStruct.hdr.qform_code;

NiiVox.img = 1; % Overwrites image, so mask 
NiiVox.hdr.dim(4:end) = 1; % remove additional MRS dimensions from header

 %transform voxel to image resolution and save under maskFile for now
nii_xform(NiiVox, NiiStruct.hdr, maskFile, 'linear', 0);

% Load maskFile into spm to adapt some fields:
vol_mask = spm_vol(maskFile);
vol_mask.dt = vol_image.dt;

%if ~isstruct(DualVoxel) 
    vol_mask.descrip = 'MRS_voxel_mask';
%else % For PRIAM data create two voxel masks
%    vol_mask.descrip = [ 'MRS_voxel_mask_' num2str(rr)];
%end

vol_mask = spm_write_vol(vol_mask,vol_mask.private.dat(:,:,:)); % write mask to vol

[T1,~] = spm_read_vols(vol_image);
T1_max = max(T1(:));

voxel_ctr = [NiiVox.hdr.qoffset_x, NiiVox.hdr.qoffset_y, NiiVox.hdr.qoffset_z];

end