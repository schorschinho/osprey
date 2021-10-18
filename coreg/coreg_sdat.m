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
% DualVoxel = struct with Dualvoxel information
%
% OUTPUTS:
% vol_mask  = SPM volume of the coregistered voxel mask.
% T1_max    = maximum intensity of the image volume.

function [vol_mask, T1_max, voxel_ctr, vol_mask_mrsi] = coreg_sdat(in, vol_image, maskFile,DualVoxel)

if nargin < 4
    DualVoxel = 0;
end

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

% For PRIAM MRS we need to add the second voxel dimensions and locations
% Depending on the PRIAM offset direction parameter, calculate the voxel
% corner and voxel center coordinates
if isstruct(DualVoxel)
    % Calculate unit vectors pointing along the voxel edges in real space
    lr_vec = vox_corner(:,1) - vox_corner(:,2);
    ap_vec = vox_corner(:,4) - vox_corner(:,1);
    cc_vec = vox_corner(:,1) - vox_corner(:,7);
    lr_vec = lr_vec ./ (norm(lr_vec));
    ap_vec = ap_vec ./ (norm(ap_vec));
    cc_vec = cc_vec ./ (norm(cc_vec));
    
    secondVoxelOffset    = DualVoxel.priam_offset;
    secondVoxelDirection = DualVoxel.priam_direction;
    switch secondVoxelDirection
        case 'R'
            shift = secondVoxelOffset*lr_vec;
        case 'L'
            shift = -secondVoxelOffset*lr_vec;
        case 'A'
            shift = secondVoxelOffset*ap_vec;
        case 'P'
            shift = -secondVoxelOffset*ap_vec;
        case 'F'
            shift = -secondVoxelOffset*cc_vec;
        case 'H'
            shift = -secondVoxelOffset*cc_vec;
    end
    vox_corner(:,:,2) = vox_corner + shift;
    vox_ctr_coor(:,:,2) = vox_ctr_coor + shift;
end



%%% 3. CREATE AND SAVE THE VOXEL MASK
% Create a mask with all voxels that are inside the voxel
for rr = 1:size(vox_ctr_coor,3)
    mask = zeros(1,size(XYZ,2));
    sphere_radius = sqrt((lr_size/2)^2+(ap_size/2)^2+(cc_size/2)^2);
    distance2voxctr = sqrt(sum((XYZ-repmat([vox_ctr_coor(1,1,rr) vox_ctr_coor(2,1,rr) vox_ctr_coor(3,1,rr)].',[1 size(XYZ,2)])).^2,1));
    sphere_mask(distance2voxctr <= sphere_radius) = 1;

    mask(sphere_mask == 1) = 1;
    XYZ_sphere = XYZ(:,sphere_mask == 1);

    tri = delaunayn([vox_corner(:,:,rr).'; [vox_ctr_coor(1,1,rr) vox_ctr_coor(2,1,rr) vox_ctr_coor(3,1,rr)]]);
    tn = tsearchn([vox_corner(:,:,rr).'; [vox_ctr_coor(1,1,rr) vox_ctr_coor(2,1,rr) vox_ctr_coor(3,1,rr)]], tri, XYZ_sphere.');
    isinside = ~isnan(tn);
    mask(sphere_mask==1) = isinside;

    % Take over the voxel dimensions from the structural
    mask = reshape(mask, vol_image.dim);
    
    if ~isstruct(DualVoxel) % For svs data create one voxel mask
        maskFileOut = maskFile;
        maskFileOut            = strrep(maskFile,'VoxelMask.nii',['Voxel_1_VoxelMask.nii']);     
        
    else  % For PRIAM data create two voxel masks        
        maskFileOut            = strrep(maskFile,'VoxelMask.nii',['Voxel_' num2str(rr) '_VoxelMask.nii']);       
    end

    % Fill in the SPM volume header information
    vol_mask.fname   = maskFileOut;
    vol_mask.dim     = vol_image.dim;
    vol_mask.dt      = vol_image.dt;
    vol_mask.mat     = vol_image.mat;
    if ~isstruct(DualVoxel) 
        vol_mask.descrip = 'MRS_voxel_mask';
    else % For PRIAM data create two voxel masks
        vol_mask.descrip = [ 'MRS_voxel_mask_' num2str(rr)];
    end
        

    % Write the SPM volume to disk
    if ~isstruct(DualVoxel) 
        vol_mask = spm_write_vol(vol_mask,mask);
    else % For PRIAM data store two voxel masks
        vol_mask_out{rr} = spm_write_vol(vol_mask,mask);
    end

    % Store voxel centre for output figure
    if ~isstruct(DualVoxel)
        voxel_ctr = [lr_off ap_off cc_off];
    else % For PRIAM data store two voxel center coordnates
        voxel_ctr(:,:,rr) = [lr_off ap_off cc_off]' + (rr -1) *shift;
    end
end

if exist('vol_mask_out','var')
    vol_mask = vol_mask_out;
end
if exist('vol_mask_out_mrsi','var')
    vol_mask_mrsi = vol_mask_out_mrsi;
else
    vol_mask_mrsi = [];
end
% Reactivate MATLAB warnings
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:qhullmx:InternalWarning');

end