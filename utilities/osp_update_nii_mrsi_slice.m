function out = osp_update_nii_mrsi_slice(in,dim,pixdim,distance)
%% nii_hdr_ext = osp_update_nii_mrsi_slice(in)
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

if nargin < 5
    recalculate_xyz = 0;
end

    out = in;
    if distance ~= 0
        out.nii_mrs.hdr.qoffset_z = out.nii_mrs.hdr.qoffset_z + distance;       % Shift according to distance
    end
    
    % Recalculate the sform from the existing qform and pixdims
    out.nii_mrs.hdr.srow_x= [0,0,0,0];
    out.nii_mrs.hdr.srow_y= [0,0,0,0];
    out.nii_mrs.hdr.srow_z= [0,0,0,0];
    
    a = sqrt(1-out.nii_mrs.hdr.quatern_b^2-out.nii_mrs.hdr.quatern_c^2-out.nii_mrs.hdr.quatern_d^2);
    b = out.nii_mrs.hdr.quatern_b;
    c = out.nii_mrs.hdr.quatern_c;
    d = out.nii_mrs.hdr.quatern_d;
    
    out.nii_mrs.hdr.srow_x(1)= (a^2+b^2-c^2-d^2)*out.nii_mrs.hdr.pixdim(2);
    out.nii_mrs.hdr.srow_x(2)= 2*(b*c-a*d)*out.nii_mrs.hdr.pixdim(3);
    out.nii_mrs.hdr.srow_x(3)= 2*(b*d+a*c)*out.nii_mrs.hdr.pixdim(4);
    out.nii_mrs.hdr.srow_x(4)= out.nii_mrs.hdr.qoffset_x;
    
    out.nii_mrs.hdr.srow_y(1)= 2*(b*c+a*d)*out.nii_mrs.hdr.pixdim(2);
    out.nii_mrs.hdr.srow_y(2)= (a^2+c^2-b^2-d^2)*out.nii_mrs.hdr.pixdim(3);
    out.nii_mrs.hdr.srow_y(3)= 2*(c*d-a*b)*out.nii_mrs.hdr.pixdim(4);
    out.nii_mrs.hdr.srow_y(4)= out.nii_mrs.hdr.qoffset_y;
    
    out.nii_mrs.hdr.srow_z(1)= 2*(b*d-a*c)*out.nii_mrs.hdr.pixdim(2);
    out.nii_mrs.hdr.srow_z(2)= 2*(c*d+a*b)*out.nii_mrs.hdr.pixdim(3);
    out.nii_mrs.hdr.srow_z(3)= (a^2+d^2-b^2-c^2)*out.nii_mrs.hdr.pixdim(4);
    out.nii_mrs.hdr.srow_z(4)= out.nii_mrs.hdr.qoffset_z;

end