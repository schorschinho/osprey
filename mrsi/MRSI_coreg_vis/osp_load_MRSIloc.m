function [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_MRSIloc(MRSCont, vertices)
    MRSIloc_vol = spm_vol(MRSCont.files_nii_MRSIloc{1});
    Coreg_img = spm_read_vols(MRSIloc_vol);
    AffineMat = MRSIloc_vol.mat;
    MRSI_vol = MRSIloc_vol;
    
    vertices_homogeneous = [vertices, ones(size(vertices, 1), 1)];
    vertices_voxel = (AffineMat \ vertices_homogeneous')';
    vertices_voxel = vertices_voxel(:, 1:3);
end