function [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_MRSIloc_rMRSI(MRSCont, vertices)

%   This function loads the MRSI localization image and transforms voxel
%   vertex coordinates from world space to image voxel space.
%
%   The function reads the NIfTI image used for MRSI localization and
%   applies the inverse affine transformation to convert vertex coordinates
%   from RAS+ world coordinates to voxel indices for overlay display.
%
%   USAGE:
%       [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_MRSIloc(MRSCont, vertices);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container with field:
%                     .files_nii_MRSIloc - Cell array containing path to
%                                          MRSI localization NIfTI file.
%       vertices    = Nx3 matrix of vertex coordinates in world space (mm),
%                     where N is the number of vertices.
%
%   OUTPUTS:
%       Coreg_img      = 3D image volume array from the NIfTI file.
%       AffineMat      = 4x4 affine transformation matrix (voxel to world).
%       vertices_voxel = Nx3 matrix of vertex coordinates in voxel space.
%       MRSI_vol       = SPM volume structure containing NIfTI header info.
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2024-01-06)
%       hzoelln2@jhmi.edu
%
%
%   HISTORY:
%       2026-01-06: First version of the code.
%%

    MRSIloc_vol = spm_vol(MRSCont.files_nii_MRSIloc{1});
    
    gunzip(fullfile(MRSCont.outputFolder, 'quickMaps', 'raw_tNAA.nii.gz'));
    MRSI_vol = spm_vol(fullfile(MRSCont.outputFolder, 'quickMaps', 'raw_tNAA.nii'));
    delete(fullfile(MRSCont.outputFolder, 'quickMaps', 'raw_tNAA.nii'));
    
    MRSIloc_vox_size = sqrt(sum(MRSIloc_vol.mat(1:3,1:3).^2));
    MRSI_vox_size = sqrt(sum(MRSI_vol.mat(1:3,1:3).^2));
    
    FOV_x = MRSI_vol.dim(1) * MRSI_vox_size(1);
    FOV_y = MRSI_vol.dim(2) * MRSI_vox_size(2);
    
    MRSIloc_resolution = min(MRSIloc_vox_size(1:2));
    
    temp_vol = MRSI_vol;
    temp_vol.dim(1) = round(FOV_x / MRSIloc_resolution);
    temp_vol.dim(2) = round(FOV_y / MRSIloc_resolution);
    
    MRSI_rot = MRSI_vol.mat(1:3,1:3);
    for i = 1:3
        MRSI_rot(:,i) = MRSI_rot(:,i) / norm(MRSI_rot(:,i));
    end
    
    temp_vol.mat = eye(4);
    temp_vol.mat(1:3,1) = MRSI_rot(:,1) * MRSIloc_resolution;
    temp_vol.mat(1:3,2) = MRSI_rot(:,2) * MRSIloc_resolution;
    temp_vol.mat(1:3,3) = MRSI_vol.mat(1:3,3);
    temp_vol.mat(1:3,4) = MRSI_vol.mat(1:3,4);
    
    [pth, nm, ext] = fileparts(MRSI_vol.fname);
    temp_fname = fullfile(pth, ['temp_ref_' nm ext]);
    
    temp_vol.fname = temp_fname;
    temp_vol = spm_create_vol(temp_vol);
    spm_write_vol(temp_vol, zeros(temp_vol.dim));
    
    flags = struct('mask', 0, 'mean', 0, 'which', 1, 'interp', 1);
    spm_reslice({temp_fname, MRSIloc_vol.fname}, flags);
    
    delete(temp_fname);
    
    AffineMat = temp_vol.mat;
    
    [MRSIlocpath, MRSIlocname, MRSIlocext] = fileparts(MRSCont.files_nii_MRSIloc{1});
    MRSIloc_MRSIspace_vol = spm_vol(fullfile(MRSIlocpath, ['r' MRSIlocname MRSIlocext]));
    Coreg_img = spm_read_vols(MRSIloc_MRSIspace_vol);
    
    vertices_homogeneous = [vertices, ones(size(vertices, 1), 1)];
    vertices_voxel = (AffineMat \ vertices_homogeneous')';
    vertices_voxel = vertices_voxel(:, 1:3);
end