load('/Volumes/Samsung/working/OspreyMRSIbeta/basissets/BASIS_Philips_Edited_se_MRSI_PRESS_GABA68_wMM.mat')
basisSetOff = BASIS;
basisSetOff.fids = basisSetOff.fids(:,:,1);
basisSetOff.specs = basisSetOff.specs(:,:,1);
basisSetOff.sz = size(basisSetOff.specs);
basisSetOff.ppm = basisSetOff.ppm;
basisDIFF = BASIS;
basisDIFF.fids = basisDIFF.fids(:,:,3);
basisDIFF.specs = basisDIFF.specs(:,:,3);
basisDIFF.sz = size(basisDIFF.specs);
fitOpts.coMM3 = '1to1GABA';
fitOpts.FWHMcoMM3 = 14;
fitOpts.CrFactor = 1;
[basisDIFF] = osp_addDiffMMPeaks(basisDIFF,basisSetOff,fitOpts);

BASIS.fids(:,37,:) = zeros(8192,1,4);
BASIS.specs(:,37,:) = zeros(8192,1,4);
BASIS.fids(:,37,3) = basisDIFF.fids(:,9);
BASIS.specs(:,37,3) = basisDIFF.specs(:,9);
BASIS.sz = size(BASIS.fids);
BASIS.nMets = BASIS.nMets + 1;
BASIS.name{37} = 'GABAp';
save('/Volumes/Samsung/working/OspreyMRSIbeta/basissets/BASIS_Philips_Edited_se_MRSI_PRESS_GABA68_wMM.mat','BASIS')
%%
% load(which('/libraries/AAL/ROI_MNI_V7_1mm_List.mat'));
atlas_path =  which('/libraries/AAL/AAL3v1_1mm.nii'); 
AALvol  = spm_vol(atlas_path);
AALimg  = AALvol.private.dat(:,:,:);

AALimg(AALimg == 121 |...
       AALimg == 123 |...
       AALimg == 125 |...
       AALimg == 127 |...
       AALimg == 129 |...
       AALimg == 131 |...
       AALimg == 133 |...
       AALimg == 135 |...
       AALimg == 137 |...
       AALimg == 139 |...
       AALimg == 141 |...
       AALimg == 143 |...
       AALimg == 145 |...
       AALimg == 147 |...
       AALimg == 149) = 121;

AALimg(AALimg == 122 |...
       AALimg == 124 |...
       AALimg == 126 |...
       AALimg == 128 |...
       AALimg == 130 |...
       AALimg == 132 |...
       AALimg == 134 |...
       AALimg == 136 |...
       AALimg == 138 |...
       AALimg == 140 |...
       AALimg == 142 |...
       AALimg == 144 |...
       AALimg == 146 |...
       AALimg == 148 |...
       AALimg == 150) = 122;

cur_id = 123;
for ii = 151 : 170
    AALimg(AALimg == ii) = cur_id;
    cur_id = cur_id +1;
end

[folder,name,ext]=fileparts(AALvol.fname);
AALvol.fname = fullfile(folder,['mod' name ext]);
AALvol          = spm_write_vol(AALvol, AALimg);

JHUvol  = spm_vol('/Volumes/Samsung/working/OspreyMRSIbeta/code/Osprey_MRSI-branch/libraries/JHU-ICBM/JHU-ICBM-labels-1mm.nii');
JHUimg  = JHUvol.private.dat(:,:,:);
JHUimg = JHUimg(1:181,1:217,1:181);

AtlasTable = readtable('/Volumes/Samsung/working/OspreyMRSIbeta/code/Osprey_MRSI-branch/libraries/JHU-ICBM/labels.csv'); 
%%
start_id = 143;
for ii= 1:48
    AALimg(AALimg == 0 & JHUimg == ii) = start_id;
    start_id = start_id + 1;
end
AALvol          = spm_write_vol(AALvol, AALimg);
%%
for ii = 139 : 186
    ROI(ii).Nom_C = ROI(ii).Nom_C{1};
    ROI(ii).Nom_L = ROI(ii).Nom_L{1};
end
%%
% max_scale = 0.1;
figure, montage(flip(rot90(MRSCont.processed.quickMapsInt.A.tLip),3),'DisplayRange',[min(MRSCont.processed.quickMapsInt.A.tLip,[],'all') max(MRSCont.processed.quickMapsInt.A.tLip,[],'all')])
title('tLip')
figure, montage(flip(rot90(MRSCont.processed.quickMapsInt.A.tNAA),3),'DisplayRange',[min(MRSCont.processed.quickMapsInt.A.tNAA,[],'all') max(MRSCont.processed.quickMapsInt.A.tNAA,[],'all')*max_scale])
title('tNAA')
figure, montage(flip(rot90(MRSCont.processed.quickMapsInt.A.tCr),3),'DisplayRange',[min(MRSCont.processed.quickMapsInt.A.tCr,[],'all') max(MRSCont.processed.quickMapsInt.A.tCr,[],'all')*max_scale])
title('tCr')
figure, montage(flip(rot90(MRSCont.processed.quickMapsInt.A.tCho),3),'DisplayRange',[min(MRSCont.processed.quickMapsInt.A.tCho,[],'all') max(MRSCont.processed.quickMapsInt.A.tCho,[],'all')*max_scale])
title('tCho')


figure, montage(flip(rot90(MRSCont.processed.quickMapsInt.w.tLip),3),'DisplayRange',[min(MRSCont.processed.quickMapsInt.w.tLip,[],'all') max(MRSCont.processed.quickMapsInt.w.tLip,[],'all')])
title('tLip from water ref')
figure, montage(flip(rot90(MRSCont.processed.quickMapsInt.w.H2O),3),'DisplayRange',[min(MRSCont.processed.quickMapsInt.w.H2O,[],'all') max(MRSCont.processed.quickMapsInt.w.H2O,[],'all')])
title('H2O')
%%
test = MRSCont.processed.AFID{1};

test.fids = squeeze(test.fids(:,12,21,3,:));
test.specs = squeeze(test.specs(:,12,21,3,:));
test.centerFreq=test.centerFreq(1);
test.dwelltime=test.dwelltime(1);
test.spectralwidth=test.spectralwidth(1);
test.txfrq=test.txfrq(1);
test.dims.extras=0;
test.dims.Xvoxels=0;
test.dims.Yvoxels=0;
test.dims.Zvoxels=0;
test.sz = size(test.fids);
% test=op_addphase(test,180,0,4.65,1);

ModelProcedure = jsonToStruct('/Volumes/Samsung/working/OspreyMRSIbeta/model-procedures/1Step_final_in_vivo_longTE.json');
if isstruct(ModelProcedure.Steps)
    ModelProcedureCell = cell(size(ModelProcedure.Steps));
    for ss = 1 : size(ModelProcedure.Steps,1)
        ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
    end
    ModelProcedure.Steps = ModelProcedureCell;
end

ModelProcedure.Steps{1, 1}.parametrizations.freqShift.lb = -1000;
ModelProcedure.Steps{1, 1}.parametrizations.freqShift.ub = 1000;
ModelProcedure.Steps{1, 1}.parametrizations.freqShift.gr.g1 = ModelProcedure.Steps{1, 1}.basisset.include;


res = Osprey_gLCM(test,ModelProcedure,0,0);
res{1}.plotFit1DStack(1,1,1)
res{1}.plotFit1DStack(1,1,1,res{1}.Options{1}.optimFreqFitRange,1)

test = op_zeropad(test,2,1);
test=op_freqshift(test,res{1, 1}.Model{1, 1}.parsOut.freqShift(1));

res = Osprey_gLCM(test,'/Volumes/Samsung/working/OspreyMRSIbeta/model-procedures/1Step_final_in_vivo_longTE.json',0,1);
res{1}.plotFit1DStack(1,1,1)
res{1}.plotFit1DStack(1,1,1,res{1}.Options{1}.optimFreqFitRange,1)
%%
test = MRSCont.processed.A{1};
vox = [1,1,1];

test.fids = squeeze(test.fids(:,vox(1),vox(2),vox(3),:));
test.specs = squeeze(test.specs(:,vox(1),vox(2),vox(3),:));
test.centerFreq=test.centerFreq(1);
test.dwelltime=test.dwelltime(1);
test.spectralwidth=test.spectralwidth(1);
test.txfrq=test.txfrq(1);
test.dims.extras=0;
test.dims.Xvoxels=0;
test.dims.Yvoxels=0;
test.dims.Zvoxels=0;
test.sz = size(test.fids);

ModelProcedure = jsonToStruct('/Volumes/Samsung/working/OspreyMRSIbeta/model-procedures/3Step_Spline_invivo_Reg_Optim_MRSI.json');
if isstruct(ModelProcedure.Steps)
    ModelProcedureCell = cell(size(ModelProcedure.Steps));
    for ss = 1 : size(ModelProcedure.Steps,1)
        ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
    end
    ModelProcedure.Steps = ModelProcedureCell;
end




res = Osprey_gLCM(test,ModelProcedure,0,0);
res{1}.plotFit1DStack(1,1,1)
res{1}.plotFit1DStack(1,2,1)
res{1}.plotFit1DStack(1,3,1)
% res{1}.plotFit1DStack(1,1,1,res{1}.Options{1}.optimFreqFitRange,1)

%%
test = MRSCont.processed.AFID{1};
vox = [17,22,3];

test.fids = squeeze(test.fids(:,vox(1),vox(2),vox(3),:));
test.specs = squeeze(test.specs(:,vox(1),vox(2),vox(3),:));
test.centerFreq=test.centerFreq(1);
test.dwelltime=test.dwelltime(1);
test.spectralwidth=test.spectralwidth(1);
test.txfrq=test.txfrq(1);
test.dims.extras=0;
test.dims.Xvoxels=0;
test.dims.Yvoxels=0;
test.dims.Zvoxels=0;
test.sz = size(test.fids);
% op_plotfid(test)

ModelProcedure = jsonToStruct('/Volumes/Samsung/working/OspreyMRSIbeta/model-procedures/1Step_final_in_vivo_longTE.json');
if isstruct(ModelProcedure.Steps)
    ModelProcedureCell = cell(size(ModelProcedure.Steps));
    for ss = 1 : size(ModelProcedure.Steps,1)
        ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
    end
    ModelProcedure.Steps = ModelProcedureCell;
end


  [~, refFWHM] = osp_XReferencing(test,[2.01,3.03,3.22],[1,1,1],[1.85,4],0);
  refFWHM = refFWHM * test.txfrq*1e-6;

  [refShift, ~] = osp_XReferencing(test,[2.01],[1],[1.85,4],0);
 
  % test.FWHM = refFWHM;


 % test=op_freqshift(test,-refShift);
ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.init = 0;
ModelProcedure.Steps{1, 1}.parametrizations.GlobFreqShift.init = 0;
% ModelProcedure.Steps{1, 1}.parametrizations.GlobFreqShift.ex = refShift;
% ModelProcedure.Steps{1, 1}.parametrizations.GlobFreqShift.sd = 0.1*abs(refShift);
ModelProcedure.Steps{1, 1}.parametrizations.GlobFreqShift.ub = 15;
ModelProcedure.Steps{1, 1}.parametrizations.GlobFreqShift.lb = -15;
% ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.init = refFWHM*2;
% ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.ub = refFWHM*2;
% ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.lb = refFWHM*2;
% ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.ex = refFWHM;
% ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.sd = 0.5 * refFWHM;



res = Osprey_gLCM(test,ModelProcedure,0,0);
res{1}.plotFit1DStack(1,1,1)
% res{1}.plotFit1DStack(1,1,1,res{1}.Options{1}.optimFreqFitRange,1)

%%

RSCont.opts.MRSI.MaxEcho.tstart = 75;

vox = [18,22,3];
test = MRSCont.processed.A{kk};
test.fids = squeeze(test.fids(:,vox(1),vox(2),vox(3),:));
test.specs = squeeze(test.specs(:,vox(1),vox(2),vox(3),:));
test.dims.extras=0;
test.dims.Xvoxels=0;
test.dims.Yvoxels=0;
test.dims.Zvoxels=0;
test.sz = size(test.fids);

% test =op_autophase(test,1.8,2.1);

[right] = op_SeparateMaxEcho(test,'right',MRSCont.opts.MRSI.MaxEcho.tstart);
[left] = op_SeparateMaxEcho(test,'flipleft',MRSCont.opts.MRSI.MaxEcho.tstart);
test = op_mergeextra(right,left,'echoside');
test.fids = squeeze(test.fids);
test.specs = squeeze(test.specs);
test.centerFreq=test.centerFreq(1);
test.dwelltime=test.dwelltime(1);
test.spectralwidth=test.spectralwidth(1);
test.txfrq=test.txfrq(1);
test.dims.extras=2;
test.sz = size(test.fids);



ModelProcedure = jsonToStruct('/Volumes/Samsung/working/OspreyMRSIbeta/model-procedures/1Step_final_in_vivo_longTE_2D_MaxEcho.json');
if isstruct(ModelProcedure.Steps)
    ModelProcedureCell = cell(size(ModelProcedure.Steps));
    for ss = 1 : size(ModelProcedure.Steps,1)
        ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
    end
    ModelProcedure.Steps = ModelProcedureCell;
end


  % [~, refFWHM] = osp_XReferencing(right,[2.01,3.03,3.22],[1,1,1],[1.85,4],0);
  % refFWHM = refFWHM * test.txfrq*1e-6;

  % [refShift, ~] = osp_XReferencing(test,[2.01],[1],[1.85,4],0);
 
  % test.FWHM = refFWHM;


 % test=op_freqshift(test,-refShift);
% ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.init = 0;
% ModelProcedure.Steps{1, 1}.parametrizations.GlobFreqShift.init = 0;
% ModelProcedure.Steps{1, 1}.parametrizations.GlobFreqShift.ex = refShift;
% ModelProcedure.Steps{1, 1}.parametrizations.GlobFreqShift.sd = 0.1*abs(refShift);
% ModelProcedure.Steps{1, 1}.parametrizations.GlobFreqShift.ub = 15;
% ModelProcedure.Steps{1, 1}.parametrizations.GlobFreqShift.lb = -15;
% ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.init = refFWHM;
% ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.ub = refFWHM*2;
% ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.lb = refFWHM*0.5;
% % ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.ex = refFWHM;
% ModelProcedure.Steps{1, 1}.parametrizations.gaussLB.sd = 0.5 * refFWHM;



res = Osprey_gLCM(test,ModelProcedure,0,0);
res{1}.plotFit1DStack(1,1,1)
res{1}.plotFit1DStack(1,1,2)
% res{1}.plotFit1DStack(1,1,1,res{1}.Options{1}.optimFreqFitRange,1)
%% Segmentation results

% Start with voxel fractions

figure
% set(gcf, 'Color', [0 0 0]);

plotMap = squeeze((MRSCont.atlas{1, 1}.fAAL(77,:,:,:))) + squeeze((MRSCont.atlas{1, 1}.fAAL(78,:,:,:)));
map_cat_fWM = flip(rot90(flip(squeeze(plotMap(:, :, :)),3)),2);

imagesc(squeeze(map_cat_fWM(:,:,4)));
colormap gray
clim([0 1]) 
axis off
hold on
%%
for i = 1 : 166
    MRSCont.atlas{1, 1}.fAAL(i,:,:,:) = MRSCont.atlas{1, 1}.fAAL(i,:,:,:)/i;
end
%%
[BASIS_LS] = create_maximumEchoBasis(BASIS,75,'left');
[BASIS_LSflip] = create_maximumEchoBasis(BASIS,75,'leftflip');
figure, plot(BASIS.ppm,real(squeeze(BASIS.specs(:,3)))), hold on
plot(BASIS_LS.ppm,real(squeeze(BASIS_LS.specs(:,3))))
plot(BASIS_LSflip.ppm,real(squeeze(BASIS_LSflip.specs(:,3))))

fit_plotBasis(BASIS_LS)
fit_plotBasis(BASIS_LSflip)
%%
[BASIS_max] = create_maximumEchoBasis(BASIS,73,'max');
figure, plot(BASIS.ppm,real(squeeze(BASIS.specs(:,3)))), hold on
plot(BASIS_max.ppm,real(squeeze(BASIS_max.specs(:,3))))

fit_plotBasis(BASIS_max)
%%
load('/Volumes/Samsung/working/OspreyMRSIbeta/basissets/BASIS_Siemens_UnEdited_sLASER_GABA135_noMM.mat')
fit_plotBasis(BASIS)
BASIS.fids(:,31)= -BASIS.fids(:,5);
BASIS.specs(:,31)= -BASIS.specs(:,5);
BASIS.name{31} = 'CrCH2';
BASIS.sz = size(BASIS.fids);
BASIS.nMets =31;
fit_plotBasis(BASIS)
save('/Volumes/Samsung/working/OspreyMRSIbeta/basissets/BASIS_Siemens_UnEdited_sLASER_GABA135_noMM.mat','BASIS');
