%%

kk = 1;
which_slice = 2;
ppmmin = 0.5;
ppmmax = 4.2;
lb = 1;
add_No_MoCo = 1;
coord = [3 6; 4 6];
yLim = [-2 1.3] * 1e-3;

which_spec = {'diff1'};
out = osp_plotProcessMRSIMatrix(MRSCont, kk, which_spec,coord, ppmmin, ppmmax, lb,add_No_MoCo,which_slice,yLim);
mkdir(fullfile(MRSCont.outputFolder,'Figures'));
saveas(out,fullfile(MRSCont.outputFolder,'Figures','diff1'),'pdf');

%%
kk = 1;
which_slice = 2;
ppmmin = 0.5;
ppmmax = 4.2;
lb = 1;
add_No_MoCo = 1;
coord = [3 6; 4 6];

yLim = [-0.5 1.3] * 1e-2;
which_spec = {'A','B'};
out = osp_plotProcessMRSIMatrix(MRSCont, kk, which_spec,coord, ppmmin, ppmmax, lb,add_No_MoCo,which_slice,yLim);
saveas(out,fullfile(MRSCont.outputFolder,'Figures','A_B'),'pdf');
%%
kk = 1;
which_slice = 2;
ppmmin = 0.5;
ppmmax = 4.2;
lb = 1;
add_No_MoCo = 0;
coord = [3 6; 4 6];

which_spec = {'A','B'};
out = osp_plotProcessMRSIMatrix(MRSCont, kk, which_spec,coord, ppmmin, ppmmax, lb,add_No_MoCo,which_slice);
mkdir(fullfile(MRSCont.outputFolder,'Figures'));
saveas(out,fullfile(MRSCont.outputFolder,'Figures','A_B_NoMoCo'),'pdf');