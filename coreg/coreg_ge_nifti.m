% coreg_p.m
% dr. Peter Van Schuerbbeek, UZ Brussel (VUB), 2020.
%
% USAGE:
% [vol_mask] = coreg_nifti(in, vol_image, maskFile);
%
% DESCRIPTION:
% Creates a SPM volume containing a voxel mask with the same dimensions
% as the SPM volume containing a structural image. The voxel mask will be
% created with the geometry information stored in the FID-A structure "in".
% This routine works for the nifti data type.

% INPUTS:
% in        = Input data structure.
% vol_image = SPM volume with the T1 data
% maskFile  = Filename under which the SPM volume of the co-registered
%           voxel mask is supposed to be saved.
%
% OUTPUTS:
% vol_mask  = SPM volume of the coregistered voxel mask.
% T1_max    = maximum intensity of the image volume.

function [vol_mask, T1_max, voxel_ctr] = coreg_ge_nifti(in, vol_image, maskFile)

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

%%% 2. PREPARE THE MRS VOXEL COORDINATES %%%
% Convert from RAS to LPS
VoxOffs  = geom.pos .* [-1 1 1];
tlhc_LPS = geom.rot.tlhc .* [-1 1 1];
trhc_LPS = geom.rot.trhc .* [-1 1 1];
brhc_LPS = geom.rot.brhc .* [-1 1 1];

e1_SVS_n = trhc_LPS - tlhc_LPS;
e1_SVS_n = e1_SVS_n ./ norm(e1_SVS_n);
e2_SVS_n = brhc_LPS - trhc_LPS;
e2_SVS_n = e2_SVS_n ./ norm(e2_SVS_n);
e3_SVS_n = -cross(e1_SVS_n, e2_SVS_n);

[~,orientation_SVS] = max(abs(e3_SVS_n));

if orientation_SVS == 3     % axial
    e1_SVS_n2 = e1_SVS_n;
    e2_SVS_n2 = e2_SVS_n;
    e3_SVS_n2 = e3_SVS_n;
elseif orientation_SVS == 2 % coronal
    e1_SVS_n2 = e1_SVS_n;
    e2_SVS_n2 = e3_SVS_n;
    e3_SVS_n2 = e2_SVS_n;
elseif orientation_SVS == 1 % sagittal
    e1_SVS_n2 = e3_SVS_n;
    e2_SVS_n2 = e1_SVS_n;
    e3_SVS_n2 = e2_SVS_n;
end

e1_SVS = geom.size.dim1 * e1_SVS_n2;
e2_SVS = geom.size.dim2 * e2_SVS_n2;
e3_SVS = geom.size.dim3 * e3_SVS_n2;

% LPS gives center of voxel
LPS_SVS_edge = VoxOffs - 0.5 * e1_SVS ...
                        - 0.5 * e2_SVS ...
                        - 0.5 * e3_SVS;
                                
dXYZ = sqrt((XYZ(1,:)-LPS_SVS_edge(1)).^2+(XYZ(2,:)-LPS_SVS_edge(2)).^2+(XYZ(3,:)-LPS_SVS_edge(3)).^2);
[~,refvox] = min(dXYZ);
[refvox_x,refvox_y,refvox_z] = ind2sub(vol_image.dim,refvox(1));

e1_edge = LPS_SVS_edge + e1_SVS;
e2_edge = LPS_SVS_edge + e2_SVS;
e3_edge = LPS_SVS_edge + e3_SVS;
                
de1XYZ = sqrt((XYZ(1,:)-e1_edge(1)).^2+(XYZ(2,:)-e1_edge(2)).^2+(XYZ(3,:)-e1_edge(3)).^2);
[~,e1vox] = min(de1XYZ);
[e1vox_x,e1vox_y,e1vox_z] = ind2sub(vol_image.dim,e1vox(1));

de2XYZ = sqrt((XYZ(1,:)-e2_edge(1)).^2+(XYZ(2,:)-e2_edge(2)).^2+(XYZ(3,:)-e2_edge(3)).^2);
[~,e2vox] = min(de2XYZ);
[e2vox_x,e2vox_y,e2vox_z] = ind2sub(vol_image.dim,e2vox(1));

de3XYZ = sqrt((XYZ(1,:)-e3_edge(1)).^2+(XYZ(2,:)-e3_edge(2)).^2+(XYZ(3,:)-e3_edge(3)).^2);
[~,e3vox] = min(de3XYZ);
[e3vox_x,e3vox_y,e3vox_z] = ind2sub(vol_image.dim,e3vox(1));

%%% 3. CREATE AND SAVE THE VOXEL MASK
% Create a mask with all voxels that are inside the voxel
mask = zeros(vol_image.dim);

nx = floor(sqrt((e1vox_x-refvox_x)^2+(e1vox_y-refvox_y)^2+(e1vox_z-refvox_z)^2))*2;
ny = floor(sqrt((e2vox_x-refvox_x)^2+(e2vox_y-refvox_y)^2+(e2vox_z-refvox_z)^2))*2;
nz = floor(sqrt((e3vox_x-refvox_x)^2+(e3vox_y-refvox_y)^2+(e3vox_z-refvox_z)^2))*2;

stepx = ([e1vox_x,e1vox_y,e1vox_z]-[refvox_x,refvox_y,refvox_z])/nx;
stepy = ([e2vox_x,e2vox_y,e2vox_z]-[refvox_x,refvox_y,refvox_z])/ny;
stepz = ([e3vox_x,e3vox_y,e3vox_z]-[refvox_x,refvox_y,refvox_z])/nz;

mrs_box_ind = [1:1:nx*ny*nz];
mrs_box_sub = zeros(3,nx*ny*nz);

[mrs_box_sub_x,mrs_box_sub_y,mrs_box_sub_z] = ind2sub([nx,ny,nz],mrs_box_ind);

mrs_box_sub(1,:)=mrs_box_sub_x;
mrs_box_sub(2,:)=mrs_box_sub_y;
mrs_box_sub(3,:)=mrs_box_sub_z;

e1_stepx = repmat(stepx,[numel(mrs_box_sub(1,:)),1])';
e2_stepy = repmat(stepy,[numel(mrs_box_sub(1,:)),1])';
e3_stepz = repmat(stepz,[numel(mrs_box_sub(1,:)),1])';

mrs_box_sub = (mrs_box_sub(1,:)-1) .* e1_stepx + (mrs_box_sub(2,:)-1) .* e2_stepy + (mrs_box_sub(3,:)-1) .* e3_stepz;

refvox_rep = repmat([refvox_x,refvox_y,refvox_z],[numel(mrs_box_sub(1,:)),1])';

mrs_box_sub = round(mrs_box_sub+refvox_rep);

mrs_box_ind = sub2ind(vol_image.dim,mrs_box_sub(1,:),mrs_box_sub(2,:),mrs_box_sub(3,:));

mask(mrs_box_ind) = 1;

%mask = flip(mask,2);

% Fill in the SPM volume header information
vol_mask.fname   = maskFile;
vol_mask.dim     = vol_image.dim;
vol_mask.dt      = vol_image.dt;
vol_mask.mat     = vol_image.mat;
vol_mask.descrip = 'MRS_voxel_mask';

% Write the SPM volume to disk
vol_mask = spm_write_vol(vol_mask,mask);

% Store voxel centre for output figure
% VoxOffs(1:2) = -VoxOffs(1:2); Apparently, this is not necessary for GE
% niftis (HZ 20201102)
voxel_ctr = VoxOffs;

% Reactivate MATLAB warnings
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:qhullmx:InternalWarning');

end