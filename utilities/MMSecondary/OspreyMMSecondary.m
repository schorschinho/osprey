function [MRSCont] = OspreyMMSecondary(MRSCont)
%% [MRSCont] = OspreyMMSecondary(MRSCont)
%   This calculates the MM concentrations from a cleaned MM spectrum using
%   the methods discribed in Hui et al. MRM 2021. It is based on the pre
%   v.3.0.0 model and the scripts that were published with the manuscript.
%
%   USAGE:
%       MRSCont = OspreyMMSecondary(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Helge Zoellner(Johns Hopkins University, 2026-05-05)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2026-05-05: First version of the code.
%% Inital checks
% You need run OspreyLoad, OspreyProcess, OspreyFit, OspreyCoreg, OspreySeg
% first

% Checking for version, toolbox, and previously run modules
[~,MRSCont.ver.CheckOsp ] = osp_CheckRunPreviousModule(MRSCont, 'OspreyLoad');
[~,MRSCont.ver.CheckOsp ] = osp_CheckRunPreviousModule(MRSCont, 'OspreyProcess');
[~,MRSCont.ver.CheckOsp ] = osp_CheckRunPreviousModule(MRSCont, 'OspreyFit');
[~,MRSCont.ver.CheckOsp ] = osp_CheckRunPreviousModule(MRSCont, 'OspreyCoreg');
[~,MRSCont.ver.CheckOsp ] = osp_CheckRunPreviousModule(MRSCont, 'OspreySeg');
%% Get phased MM spectrum

for kk = 1:MRSCont.nDatasets(1) %Subject loop
    basisSet    = MRSCont.fit.resBasisSet.mm.(MRSCont.info.mm.unique_ndatapoint_spectralwidth{1}){1}; 
    dataToPlot  = MRSCont.processed.mm{kk}; 
    fitParams   = MRSCont.fit.results.mm.fitParams{kk};
    inputData.dataToFit                 = op_takesubspec(dataToPlot,1);
    inputData.basisSet                  = basisSet;
    inputSettings.scale                 = MRSCont.fit.scale{kk};
    inputSettings.fitRangePPM           = [0.2 4.2];
    inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;
    inputSettings.fitStyle              = MRSCont.opts.fit.style;
    inputSettings.flags.isMEGA          = MRSCont.flags.isMEGA;
    inputSettings.flags.isHERMES        = MRSCont.flags.isHERMES;
    inputSettings.flags.isHERCULES      = MRSCont.flags.isHERCULES;
    inputSettings.flags.isPRIAM         = MRSCont.flags.isPRIAM;

    [ModelOutput]=fit_OspreyParamsToModel_nophase(inputData, inputSettings, fitParams);
    fit=ModelOutput.completeFit;
    baseline=ModelOutput.baseline;
    ppm=ModelOutput.ppm;
    res=ModelOutput.residual;
    data=ModelOutput.data;
    MM_clean(kk,:) = ModelOutput.data -sum(ModelOutput.indivMets(:,1:4),2);
    MRSCont.MMSecondary.clean_MM_specs(kk,:)=MM_clean(kk,:)*MRSCont.fit.scale{kk};

end
MRSCont.MMSecondary.ppm = ppm;
%% Get water area for weighting
for kk = 1:MRSCont.nDatasets(1) %Subject loop
    fitRangePPM2 = MRSCont.opts.fit.rangeWater; 
    basisSet    = MRSCont.fit.resBasisSet.ref.(MRSCont.info.ref.unique_ndatapoint_spectralwidth{1}){1}; 
    dataToPlot  = MRSCont.processed.ref{kk}; 
    fitParams   = MRSCont.fit.results.ref.fitParams{kk}; 
    inputData.dataToFit                 = dataToPlot;
    inputData.basisSet                  = basisSet;
    inputSettings.scale                 = MRSCont.fit.scale{kk};
    inputSettings.fitRangePPM           = fitRangePPM2;
    inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;
    inputSettings.fitStyle              = MRSCont.opts.fit.style;
    inputSettings.flags.isMEGA          = MRSCont.flags.isMEGA;
    inputSettings.flags.isHERMES        = MRSCont.flags.isHERMES;
    inputSettings.flags.isHERCULES      = MRSCont.flags.isHERCULES;
    inputSettings.flags.isPRIAM         = MRSCont.flags.isPRIAM;

    [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
    fit      = ModelOutput.completeFit;
    ppm      =  ModelOutput.ppm;
    
    fit  = fit * MRSCont.fit.scale{kk};
    MRSCont.MMSecondary.waterArea(kk) = sum(fit);
end
%% Apply relaxation and waterArea weighting to the spectra
for kk = 1:MRSCont.nDatasets(1) %Subject loop

    % Water relaxation 3 T
    % From Lu et al. 2005 (JMRI)
    % CSF T1 = 3817 +/- 424msec - but state may underestimated and that 4300ms
    % is likely more accurate - but the reference is to an ISMRM 2001 abstract
    % MacKay (last author) 2006 ISMRM abstract has T1 CSF = 3300 ms
    % CSF T2 = 503.0 +/- 64.3 Piechnik MRM 2009; 61: 579
    % However, other values from Stanisz et al:
    % CPMG for T2, IR for T1
    % T2GM = 99 +/ 7, lit: 71+/- 10 (27)
    % T1GM = 1820 +/- 114, lit 1470 +/- 50 (29)
    % T2WM = 69 +/-3 lit 56 +/- 4 (27)
    % T1WM = 1084 +/- 45 lit 1110 +/- 45 (29) 
    T1w_WM    = 0.832;
    T2w_WM    = 0.0792;
    T1w_GM    = 1.331;
    T2w_GM    = 0.110;
    T1w_CSF   = 3.817;
    T2w_CSF   = 0.503;

    % Determine concentration of water in GM, WM and CSF
    % Gasparovic et al. 2006 (MRM) uses relative densities, ref to
    % Ernst et al. 1993 (JMR)
    % fGM = 0.78
    % fWM = 0.65
    % fCSF = 0.97
    % such that
    % concw_GM = 0.78 * 55.51 mol/kg = 43.30
    % concw_WM = 0.65 * 55.51 mol/kg = 36.08
    % concw_CSF = 0.97 * 55.51 mol/kg = 53.84
    concW_GM    = 43.30*1e3;
    concW_WM    = 36.08*1e3;
    concW_CSF   = 53.84*1e3;
    molal_concW = 55.51*1e3;

    % Segmentation results
    f_CSF=MRSCont.seg.tissue.fCSF(kk);
    f_GM=MRSCont.seg.tissue.fGM(kk);
    f_WM=MRSCont.seg.tissue.fWM(kk);
    
    % Gasparovic et al. method
    % Calculate molal fractions from volume fractions (equivalent to eqs. 5-7 in Gasparovic et al., 2006)

    molal_fGM  = (f_GM*concW_GM) ./ (f_GM*concW_GM + f_WM*concW_WM + f_CSF*concW_CSF);
    molal_fWM  = (f_WM*concW_WM) ./ (f_GM*concW_GM + f_WM*concW_WM + f_CSF*concW_CSF);
    molal_fCSF = (f_CSF*concW_CSF) ./ (f_GM*concW_GM + f_WM*concW_WM + f_CSF*concW_CSF);
    
    % Get sequence parameter
    waterTR=MRSCont.processed.ref{kk}.tr * 1e-3;
    waterTE=MRSCont.processed.ref{kk}.te * 1e-3;

    % Calculate relaxation-weighting
    Relaxation_weighting= molal_concW * (molal_fGM  * (1 - exp(-waterTR/T1w_GM)) * exp(-waterTE/T2w_GM)  + ...
                molal_fWM  * (1 - exp(-waterTR/T1w_WM)) * exp(-waterTE/T2w_WM)  + ...
                molal_fCSF * (1 - exp(-waterTR/T1w_CSF)) * exp(-waterTE/T2w_CSF));

    
    MRSCont.MMSecondary.clean_MM_specs_scaled(kk,:)=MRSCont.MMSecondary.clean_MM_specs(kk,:)*Relaxation_weighting/MRSCont.MMSecondary.waterArea(kk)*2;
    MRSCont.MMSecondary.Relaxation_weighting(kk) = Relaxation_weighting;
end
%% Calculate concentration estimates
% Run fit on averaged spectrum first for initals. Then run per subejct fit

MRSCont.MMSecondary.Mean_spec=mean(real(MRSCont.MMSecondary.clean_MM_specs_scaled),1);
input = [0.9  -80.0000    0.55  -80.0000    0.75  -80.0000    0.8  -50.0000    1.3  -50.0000...
    0.85  -50.0000    0.4  -80.0000    0.55  -80.0000    0.2588  -80.0000    0.5576  -60.0000...
    0.15  -60.0000    0.15  -60.0000    0.5576  -60.0000    0.2588  -60.0000   -0.0006    0.0294...
   -0.0400];
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts, 'MaxIter', 1e5);

% % For debugging
% % Plot initals
% figure,
% plot(MRSCont.MMSecondary.ppm, MRSCont.MMSecondary.Mean_spec,MRSCont.MMSecondary.ppm, MM_spectrum_model(input,MRSCont.MMSecondary.ppm))

% Run initial fit on mean spectrum
lb = 1.0e+04 *[ 0   -1  0   -1 0   -1 0   -1 0   -1 0   -1 0   -1  0   -1  0 -1  0   -1 0 -1 0  -1 0 -1 0   -1   -1   -0.001 -0.001];
ub = 1.0e+04 *[0.01 -0.002 0.01 -0.002 0.01 -0.002 0.01 -0.002 0.01 -0.002 0.01 -0.002 0.01 -0.002 0.01 -0.002 0.01 -0.002 0.01 -0.002 0.01 -0.002 0.01 -0.002  0.01  -0.002 0.01 -0.002 1 0.001 0.001];
options = optimset('lsqcurvefit');
 options = optimset(options,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',10000);
 
output = lsqcurvefit(@(xdummy,ydummy) ...
                MM_spectrum_model(xdummy,ydummy),...
                input, MRSCont.MMSecondary.ppm,MRSCont.MMSecondary.Mean_spec,...
                lb,ub,options);

MRSCont.MMSecondary.initalFitPars = output;
MRSCont.MMSecondary.initalFit = MM_spectrum_model(output,MRSCont.MMSecondary.ppm);

% For debugging
% Plot mean spectrum fit
% figure,
% plot(MRSCont.MMSecondary.ppm,MM_spectrum_model(output,MRSCont.MMSecondary.ppm),MRSCont.MMSecondary.ppm,MRSCont.MMSecondary.Mean_spec);

% Run individual fits
for kk = 1:MRSCont.nDatasets(1) %Subject loop
    MRSCont.MMSecondary.FitPars(kk,:) = lsqcurvefit(@(xdummy,ydummy) ...
                                            MM_spectrum_model(xdummy,ydummy),...
                                            output, MRSCont.MMSecondary.ppm,MRSCont.MMSecondary.clean_MM_specs_scaled(kk,:),...
                                            lb,ub,options);
    MRSCont.MMSecondary.MM_specs_fits(kk,:) = MM_spectrum_model(MRSCont.MMSecondary.FitPars(kk,:),MRSCont.MMSecondary.ppm);
end

% Get areas
for kk = 1:MRSCont.nDatasets(1) %Subject loop
   MRSCont.MMSecondary.Areas(kk,1)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,1:2)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,2)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,3:4)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,3)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,5:6)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,4)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,7:8)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,5)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,9:10)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,6)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,11:12)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,7)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,13:14)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,8)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,15:16)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,9)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,17:18)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,10)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,19:20)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,11)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,21:22)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,12)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,23:24)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,13)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,25:26)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
   MRSCont.MMSecondary.Areas(kk,14)=sum(GaussModel([MRSCont.MMSecondary.FitPars(kk,27:28)  2.5 0 0 ],MRSCont.MMSecondary.ppm));
end

%% Save results
% Save the output structure to the output folder
% Determine output folder
outputFolder = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end
end