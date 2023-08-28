function out = osp_shift_nii_volume(in,voxelshift)
%% nii_hdr_ext = osp_shift_nii_volume(in)
%   This function updates the afine transformation matrix accoring to the
%   slice. This is useful for multi-slice MRSI with a GAP parameter. 
%
%   USAGE:
%       out             = osp_update_nii_mrsi_slice(in)
%
%   INPUTS:
%       in              = input data in matlab structure format.
%
%   OUTPUTS:
%       out             = output data in matlab structure format.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-06)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%       This code is also based on numerous functions from the spec2nii
%       conversion tool by Dr. William Clarke (University of Oxford)
%       https://github.com/wtclarke/spec2nii
%       Clarke et al., Magn Reson Med 88:2358-2370 (2022)
%
%   HISTORY:
%       2023-03-06: First version of the code.

if nargin < 2
    voxelshift = [0 0 0];
end

    out = in;

    m44 = [out.nii_mrs.hdr.srow_x; out.nii_mrs.hdr.srow_y; out.nii_mrs.hdr.srow_z];
    v = voxelshift * m44(1:3, 1:3)';
    m44(1:3, 4) = m44(1:3, 4) + v';
    
    out.nii_mrs.hdr.srow_x = m44(1,:);
    out.nii_mrs.hdr.srow_y = m44(2,:);
    out.nii_mrs.hdr.srow_z = m44(3,:);
    
    [out.nii_mrs.hdr.quatern_b, out.nii_mrs.hdr.quatern_c, out.nii_mrs.hdr.quatern_d,...
        out.nii_mrs.hdr.qoffset_x, out.nii_mrs.hdr.qoffset_y, out.nii_mrs.hdr.qoffset_z,...
        out.nii_mrs.hdr.pixdim(2), out.nii_mrs.hdr.pixdim(3), out.nii_mrs.hdr.pixdim(4),...
        out.nii_mrs.hdr.pixdim(1)] = nifti_mat44_to_quatern(m44);


end