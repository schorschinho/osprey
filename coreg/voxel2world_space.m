function [img_t,img_c,img_s] = voxel2world_space(V,voxel_ctr)

% MM (180221): Code borrowed from spm_orthviews.m

bb   = spm_get_bbox(V,'fv');
dims = round(diff(bb)'+1);
cent = voxel_ctr';

M = V.mat;
TM0 = [1 0 0 -bb(1,1)+1
       0 1 0 -bb(1,2)+1
       0 0 1 -cent(3)
       0 0 0 1];
TM = inv(TM0*M);
TD = dims([1 2]);

CM0 = [1 0 0 -bb(1,1)+1
       0 0 1 -bb(1,3)+1
       0 1 0 -cent(2)
       0 0 0 1];
CM = inv(CM0*M);
CD = dims([1 3]);


SM0 = [0 -1 0 +bb(2,2)+1
       0  0 1 -bb(1,3)+1
       1  0 0 -cent(1)
       0  0 0 1];
SM = inv(SM0*M);
SD = dims([2 3]);

img_t = spm_slice_vol(V,TM,TD,2)';
img_c = spm_slice_vol(V,CM,CD,2)';
img_s = spm_slice_vol(V,SM,SD,2)';



