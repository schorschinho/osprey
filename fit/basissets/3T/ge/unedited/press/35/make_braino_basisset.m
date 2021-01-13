function make_braino_basisset

ref_basis = io_readlcmraw_basis('basis_braino.basis');
load('ref_water.mat')
load('basis_ge_press35.mat')

nBASIS.spectralwidth = ref_basis.Cho.spectralwidth;
nBASIS.dwelltime = ref_basis.Cho.dwelltime;
nBASIS.n = ref_basis.Cho.n;
nBASIS.linewidth = ref_basis.Cho.linewidth;
nBASIS.Bo = ref_basis.Cho.Bo;
nBASIS.seq = {'press'};
nBASIS.te = ref_basis.Cho.te;
nBASIS.centerFreq = 4.83;
nBASIS.ppm = ref_basis.Cho.ppm';
nBASIS.t = ref_basis.Cho.t';
nBASIS.flags = ref_basis.Cho.flags;
nBASIS.nMM = 0;
nBASIS.dims = ref_basis.Cho.dims;
nBASIS.nMets = 8;
nBASIS.sz = [nBASIS.n,nBASIS.nMets+nBASIS.nMM];

nBASIS.fids = zeros(nBASIS.sz);
nBASIS.specs = zeros(nBASIS.sz);

nBASIS.name = {'PCh','Cr','CrCH2','Glu','Ins','Lac','NAA','H2O'};

nBASIS.fids(:,1) = ref_basis.Cho.fids;
nBASIS.fids(:,2) = ref_basis.Cr.fids;
nBASIS.fids(:,3) = ref_basis.CrCH2.fids;
nBASIS.fids(:,4) = ref_basis.Glu.fids;
nBASIS.fids(:,5) = ref_basis.Ins.fids;
nBASIS.fids(:,6) = ref_basis.Lac.fids;
nBASIS.fids(:,7) = ref_basis.NAA.fids;
%nBASIS.fids(:,8) = ref_water.fids;

nBASIS.specs(:,1) = ref_basis.Cho.specs;
nBASIS.specs(:,2) = ref_basis.Cr.specs;
nBASIS.specs(:,3) = ref_basis.CrCH2.specs;
nBASIS.specs(:,4) = ref_basis.Glu.specs;
nBASIS.specs(:,5) = ref_basis.Ins.specs;
nBASIS.specs(:,6) = ref_basis.Lac.specs;
nBASIS.specs(:,7) = ref_basis.NAA.specs;
%nBASIS.specs(:,8) = ref_water.specs;

out = op_addphase(ref_water,-7,0,4.83,0);

nBASIS.fids(:,8) = out.fids*0.46;
nBASIS.specs(:,8) = flipud(out.specs)*0.46;

BASIS = fit_resampleBasis(BASIS,nBASIS);

BASIS.ppm = BASIS.ppm';

save('basis_ge_press35_braino.mat','BASIS');