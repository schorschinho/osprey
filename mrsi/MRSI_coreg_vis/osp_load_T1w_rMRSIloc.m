function [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_T1w_rMRSIloc(MRSCont, vertices)
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
