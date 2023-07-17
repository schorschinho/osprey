% This script generates the original Osprey LCM (without the arbitrary lineshape function) according to the new format
%% Step 1 Referencing 
% LCM.Steps{1}.module = 'Referencing';
% LCM.Steps{1}.frequencies = [2.01, 3.00, 3.22];
% LCM.Steps{1}.polarity = [1 1 1];
% LCM.Steps{1}.ppm_limits = [1.85, 4.2];
% LCM.Steps{1}.realpart = 0;
% LCM.Steps{1}.extra = 1;

%% Step 2 Preliminary fit to generate initial guess
LCM.basisset.file = {'/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA30_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA50_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA70_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA90_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA110_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA130_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA150_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA170_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA190_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA210_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA230_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA250_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA270_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA290_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA310_wMM.mat',...
                      '/Volumes/Samsung/working/Shandong2/MSM/basissets/sLASER/BASIS_Philips_UnEdited_sLASER_GABA330_wMM.mat'};
LCM.Steps(1).module = 'Model';
LCM.Steps(1).ModelFunction = 'BasicPhysicsModel';
LCM.Steps(1).extra.flag = 1;
LCM.Steps(1).extra.linkedParameter = {'metAmpl'};
LCM.Steps(1).basisset.spec = 1;
LCM.Steps(1).basisset.include = {'Cr'};
LCM.Steps(1).fit_opts.ppm = [0.5 4];
LCM.Steps(1).fit_opts.gap = [];
LCM.Steps(1).fit_opts.optimDomain = 'FD';
LCM.Steps(1).fit_opts.optimSignalPart = 'R';
LCM.Steps(1).fit_opts.baseline.type = 'spline';
LCM.Steps(1).fit_opts.baseline.dkntmn = [0.5];
LCM.Steps(1).fit_opts.solver = 'lsqnonlin';
LCM.Steps(1).parameter = {'ph0','ph1','gaussLB','lorentzLB',...
                            'freqShift','metAmpl','baseAmpl'};
LCM.Steps(1).parametrizations.gaussLB.ub = 'Inf';
LCM.Steps(1).parametrizations.lorentzLB.ub = 'Inf';

%%
fid=fopen(['/Volumes/Samsung/working/Shandong2/MSM/code/ModelProcedures/TE_Series_1Step_Spline.json'],'w');
fprintf(fid, jsonencode(LCM, "PrettyPrint", true)); 
fclose(fid);