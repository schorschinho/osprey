function [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_T1w_rMRSIloc(MRSCont, vertices)
%% [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_T1w_rMRSIloc(MRSCont, vertices)
%   This function reslices a T1-weighted structural image to match the MRSI
%   localization image geometry for overlay visualization.
%
%
%   USAGE:
%       [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_T1w_rMRSIloc(MRSCont, vertices);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container with fields:
%                     .files_nii         - Cell array containing path to T1 NIfTI
%                     .files_nii_MRSIloc - Cell array containing path to MRSI
%                                          localization NIfTI file
%       vertices    = Nx3 matrix of vertex coordinates in world space (mm),
%                     where N is the number of vertices.
%
%   OUTPUTS:
%       Coreg_img      = 3D image volume array of resliced T1 in MRSI space.
%       AffineMat      = 4x4 affine transformation matrix of MRSI localization.
%       vertices_voxel = Nx3 matrix of vertex coordinates in voxel space.
%       MRSI_vol       = SPM volume structure of MRSI localization image.
%
%   AUTHOR:
%       Dr. Helge Zollner (Johns Hopkins University, 2024-01-06)
%       hzoelln2@jhmi.edu
%
%
%   HISTORY:
%       2026-01-06: First version of the code.
%%
    T1_struc_vol = spm_vol(MRSCont.files_nii{1});
    MRSIloc_vol = spm_vol(MRSCont.files_nii_MRSIloc{1});
    
    spm_reslice({MRSIloc_vol.fname, T1_struc_vol.fname}, struct('mask', 0, 'mean', 0, 'which', 1));
    
    AffineMat = MRSIloc_vol.mat;
    MRSI_vol = MRSIloc_vol;
    
    vertices_homogeneous = [vertices, ones(size(vertices, 1), 1)];
    vertices_voxel = (AffineMat \ vertices_homogeneous')';
    vertices_voxel = vertices_voxel(:, 1:3);
    
    [T1path, T1name, T1ext] = fileparts(MRSCont.files_nii{1});
    T1_struc_MRSIspace_vol = spm_vol(fullfile(T1path, ['r' T1name T1ext]));
    Coreg_img = spm_read_vols(T1_struc_MRSIspace_vol);
end
