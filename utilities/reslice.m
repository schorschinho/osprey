function reslice(PI,PO,dim,mat,hld)
% FORMAT reslice(PI,PO,dim,mat,hld)
%   PI - input filename
%   PO - output filename
%   dim - 1x3 matrix of image dimensions
%   mat - 4x4 affine transformation matrix mapping
%         from vox to mm (for output image).
%         To define M from vox and origin, then
%             off = -vox.*origin;
%              M   = [vox(1) 0      0      off(1)
%                     0      vox(2) 0      off(2)
%                     0      0      vox(3) off(3)
%                     0      0      0      1];
%
%   hld - interpolation method.
%___________________________________________________________________________
% @(#)JohnsGems.html	1.42 John Ashburner 05/02/02
 
VI          = spm_vol(PI);
VO          = VI;
VO.fname    = deblank(PO);
VO.mat      = mat;
VO.dim(1:3) = dim;
 
VO = spm_create_vol(VO); 
for x3 = 1:VO.dim(3)
        M  = inv(spm_matrix([0 0 -x3 0 0 0 1 1 1])*inv(VO.mat)*VI.mat);
        v  = spm_slice_vol(VI,M,VO.dim(1:2),hld);
        VO = spm_write_plane(VO,v,x3);
end
end