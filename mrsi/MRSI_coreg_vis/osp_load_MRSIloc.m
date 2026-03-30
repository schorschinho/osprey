function [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_MRSIloc(MRSCont, vertices)
%% [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_MRSIloc(MRSCont, vertices)
%   This function loads the MRSI localization image and transforms vertex
%   coordinates from world space to voxel space for overlay display.
%
%   The function reads the NIfTI image acquired for MRSI localization and
%   applies the inverse affine transformation to convert vertex coordinates
%   from RAS+ world coordinates (mm) to voxel indices.
%
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
%       Coreg_img      = 3D image volume array from the localization NIfTI.
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
    Coreg_img = spm_read_vols(MRSIloc_vol);
    AffineMat = MRSIloc_vol.mat;
    MRSI_vol = MRSIloc_vol;
    
    vertices_homogeneous = [vertices, ones(size(vertices, 1), 1)];
    vertices_voxel = (AffineMat \ vertices_homogeneous')';
    vertices_voxel = vertices_voxel(:, 1:3);
end