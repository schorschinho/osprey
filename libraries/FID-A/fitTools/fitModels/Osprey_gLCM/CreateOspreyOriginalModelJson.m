% This script generates the original Osprey LCM (without the arbitrary lineshape function) according to the new format
%% Step 1 Referencing 
% LCM.Steps{1}.module = 'Referencing';
% LCM.Steps{1}.frequencies = [2.01, 3.00, 3.22];
% LCM.Steps{1}.polarity = [1 1 1];
% LCM.Steps{1}.ppm_limits = [1.85, 4.2];
% LCM.Steps{1}.realpart = 0;
% LCM.Steps{1}.extra = 1;

%% Step 2 Preliminary fit to generate initial guess
LCM.Steps{1}.module = 'Model';
LCM.Steps{1}.extra = 1;
LCM.Steps{1}.basisset.file = '/Users/helge/Documents/GitHub/osprey/fit/basissets/3T/philips/unedited/press/35/basis_philips_press35.mat';
LCM.Steps{1}.basisset.spec = 1;
LCM.Steps{1}.basisset.include = {'Cr'};
LCM.Steps{1}.fit_opts.ppm = [0.5 4];
LCM.Steps{1}.fit_opts.gap = [];
LCM.Steps{1}.fit_opts.optimDomain = 'FD';
LCM.Steps{1}.fit_opts.optimSignalPart = 'RI';
LCM.Steps{1}.fit_opts.baseline.type = 'spline';
LCM.Steps{1}.fit_opts.baseline.dkntmn = [0.5];

%% Step 3 Final fit 

LCM.Steps{2}.module = 'Model';
LCM.Steps{2}.extra = 1;
LCM.Steps{2}.basisset.file = '/Users/helge/Documents/GitHub/osprey/fit/basissets/3T/philips/unedited/press/35/basis_philips_press35.mat';
LCM.Steps{2}.basisset.spec = 1;
LCM.Steps{2}.basisset.include = {'Cr'};
LCM.Steps{2}.fit_opts.ppm = [0.5 4];
LCM.Steps{2}.fit_opts.gap = [];
LCM.Steps{2}.fit_opts.optimDomain = 'FD';
LCM.Steps{2}.fit_opts.optimSignalPart = 'R';
LCM.Steps{2}.fit_opts.baseline.type = 'spline';
LCM.Steps{2}.fit_opts.baseline.dkntmn = [0.5];
LCM.Steps{2}.parameter = {'ph0','ph1','gaussLB','lorentzLB_mets','lorentzLB_MM',...
                            'freqShift_mets','freqShift_MM','ampl','baseline'};
LCM.Steps{2}.initials.ph0 = {'Step 1'};
LCM.Steps{2}.initials.ph1 = {'Step 1'};
LCM.Steps{2}.initials.gaussLB = {'Step 1'};
LCM.Steps{2}.initials.lorentzLB = [0];
LCM.Steps{2}.initials.freqShift = [0];
LCM.Steps{2}.initials.ampl = [0];
LCM.Steps{2}.initials.baseline = [0];

%%
fid=fopen(['/Volumes/Samsung/working/Shandong2/MSM/code/DefaultOsprey.json'],'w');
fprintf(fid, jsonencode(LCM, "PrettyPrint", true)); 
fclose(fid);
