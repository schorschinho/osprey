% coreg_sdat.m
% Georg Oeltzschner, Johns Hopkins University 2019.
% 
% USAGE:
% [vol_mask] = coreg_sdat(in, vol_image, maskFile);
% 
% DESCRIPTION:
% Creates a SPM volume containing a voxel mask with the same dimensions 
% as the SPM volume containing a structural image. The voxel mask will be
% created with the geometry information stored in the FID-A structure "in".
% This routine works for the Philips SDAT format.
% 
% CREDITS:
% The routine for correct determination of the phase and readout
% directions of the MRS voxel is adapted from 
% vox2ras_rsolveAA.m
% (Dr. Rudolph Pienaar, Massachusetts General Hospital, Boston)
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

function [vol_mask, T1_max] = coreg_sdat(in, vol_image, maskFile)

% Deactivate MATLAB warnings and load geometry parameters
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:qhullmx:InternalWarning');
geom = in.geometry;

%%% 1. PREPARE THE STRUCTURAL IMAGE
% Create SPM volume and read in the NIfTI file with the structural image.
[T1,XYZ]    = spm_read_vols(vol_image);
T1_max      = max(T1(:));
%Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
%tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
[~,voxdim] = spm_get_bbox(vol_image,'fv');
voxdim = abs(voxdim)';
halfpixshift = -voxdim(1:3)/2;
halfpixshift(3) = -halfpixshift(3);
XYZ = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);


%%% 2. GENERATE THE COORDINATES OF THE VOXEL CORNERS
% Get information from SPAR - change later to be read in
ap_size = geom.size.ap;
lr_size = geom.size.lr;
cc_size = geom.size.cc;
ap_off  = geom.pos.ap;
lr_off  = geom.pos.lr;
cc_off  = geom.pos.cc;
ap_ang  = geom.rot.ap;
lr_ang  = geom.rot.lr;
cc_ang  = geom.rot.cc;

% We need to flip ap and lr axes to match NIFTI convention
ap_off = -ap_off;
lr_off = -lr_off;
ap_ang = -ap_ang;
lr_ang = -lr_ang;

% Define voxel coordinates before rotation and transition
vox_ctr = ...
    [lr_size/2 -ap_size/2  cc_size/2;
    -lr_size/2 -ap_size/2  cc_size/2;
    -lr_size/2  ap_size/2  cc_size/2;
     lr_size/2  ap_size/2  cc_size/2;
    -lr_size/2  ap_size/2 -cc_size/2;
     lr_size/2  ap_size/2 -cc_size/2;
     lr_size/2 -ap_size/2 -cc_size/2;
    -lr_size/2 -ap_size/2 -cc_size/2];

% Make rotations on voxel
rad = pi/180;
initrot = zeros(3,3);

xrot      = initrot;
xrot(1,1) = 1;
xrot(2,2) = cos(lr_ang *rad);
xrot(2,3) = -sin(lr_ang*rad);
xrot(3,2) = sin(lr_ang*rad);
xrot(3,3) = cos(lr_ang*rad);

yrot      = initrot;
yrot(1,1) = cos(ap_ang*rad);
yrot(1,3) = sin(ap_ang*rad);
yrot(2,2) = 1;
yrot(3,1) = -sin(ap_ang*rad);
yrot(3,3) = cos(ap_ang*rad);

zrot      = initrot;
zrot(1,1) = cos(cc_ang*rad);
zrot(1,2) = -sin(cc_ang*rad);
zrot(2,1) = sin(cc_ang*rad);
zrot(2,2) = cos(cc_ang*rad);
zrot(3,3) = 1;

% Apply rotation as prescribed
vox_rot = xrot * yrot * zrot * vox_ctr.';

% Shift rotated voxel by the center offset to its final position
vox_ctr_coor = [lr_off ap_off cc_off];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner = vox_rot+vox_ctr_coor;


%%% 3. CREATE AND SAVE THE VOXEL MASK
% Create a mask with all voxels that are inside the voxel
mask = zeros(1,size(XYZ,2));
sphere_radius = sqrt((lr_size/2)^2+(ap_size/2)^2+(cc_size/2)^2);
distance2voxctr = sqrt(sum((XYZ-repmat([lr_off ap_off cc_off].',[1 size(XYZ,2)])).^2,1));
sphere_mask(distance2voxctr <= sphere_radius) = 1;

mask(sphere_mask == 1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);

tri = delaunayn([vox_corner.'; [lr_off ap_off cc_off]]);
tn = tsearchn([vox_corner.'; [lr_off ap_off cc_off]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask==1) = isinside;

% Take over the voxel dimensions from the structural
mask = reshape(mask, vol_image.dim);

% Fill in the SPM volume header information
vol_mask.fname   = maskFile;
vol_mask.dim     = vol_image.dim;
vol_mask.dt      = vol_image.dt;
vol_mask.mat     = vol_image.mat;
vol_mask.pinfo   = vol_image.pinfo;
vol_mask.n       = vol_image.n;
vol_mask.descrip = 'MRS_voxel_mask';
vol_mask.private = vol_image.private;

% Write the SPM volume to disk
vol_mask = spm_write_vol(vol_mask,mask);

% Reactivate MATLAB warnings
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:qhullmx:InternalWarning');

end