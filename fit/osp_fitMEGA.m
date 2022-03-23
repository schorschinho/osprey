function [MRSCont] = osp_fitMEGA(MRSCont)
%% [MRSCont] = osp_fitMEGA(MRSCont)
%   This function performs spectral fitting of MEGA-edited MRS data.
%
%   USAGE:
%       [MRSCont] = osp_fitMEGA(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-04-12)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-04-12: First version of the code.


% Loop over all the datasets here
metFitTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end
for kk = 1:MRSCont.nDatasets
    [~] = printLog('OspreyFit',kk,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
    
    %%% 1. DETERMINE THE FITTING STYLE %%%
    % Extract fit options
    fitOpts     = MRSCont.opts.fit;
    fitModel    = fitOpts.method;
    fitStyle    = fitOpts.style;

    %%% 2. SEPARATE FIT %%%
    % For the separate (classic) MEGA fit, model the EDIT-OFF and DIFF
    % spectra separately.
    if strcmp(fitStyle, 'Separate')
        if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.off.fitParams))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)  
            %%% 2a. FIT OFF-SPECTRUM
            % Apply scaling factor to the data
            dataToFit   = MRSCont.processed.A{kk};
            basisSetOff = MRSCont.fit.basisSet;
            basisSetOff.fids = basisSetOff.fids(:,:,1);
            basisSetOff.specs = basisSetOff.specs(:,:,1);
            dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});

            % Call the fit function
            [fitParamsOff, resBasisSetOff] = fit_runFit(dataToFit, basisSetOff, fitModel, fitOpts);

            % Save back the basis set and fit parameters to MRSCont
            MRSCont.fit.resBasisSet.off{kk}             = resBasisSetOff;
            MRSCont.fit.results.off.fitParams{kk}   = fitParamsOff;


            %%% 2b. FIT DIFF1-SPECTRUM
            % Apply scaling factor to the data
            fitOpts     = MRSCont.opts.fit;
            dataToFit   = MRSCont.processed.diff1{kk};
            basisSetDiff1 = MRSCont.fit.basisSet;
            basisSetDiff1.fids = basisSetDiff1.fids(:,:,3);
            basisSetDiff1.specs = basisSetDiff1.specs(:,:,3);
            dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
            dataToFit.refShift   = fitParamsOff.refShift;
            dataToFit.refFWHM   = fitParamsOff.refFWHM;
            
            if ~strcmp(fitOpts.coMM3, 'none')
                [basisSetDiff1] = osp_addDiffMMPeaks(basisSetDiff1,basisSetOff,fitOpts);
            end

            % Call the fit function
            [fitParamsDiff1, resBasisSetDiff1]  = fit_runFit(dataToFit, basisSetDiff1, fitModel, fitOpts);

            % Save back the basis set and fit parameters to MRSCont
            MRSCont.fit.resBasisSet.diff1{kk}           = resBasisSetDiff1;
            MRSCont.fit.results.diff1.fitParams{kk} = fitParamsDiff1;
        end
    end


    %%% 3. CONCATENATED FIT %%%
    % For the concatenated MEGA fit, model the DIFF1 and SUM spectra
    % together.
    if strcmp(fitStyle, 'Concatenated')
        if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.conc.fitParams))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)          
            %%% 3a. FIT CONCATENATED SPECTRUM
            % Select the difference and sum spectrum to put into the
            % concatenated fit
            dataToFit   = {MRSCont.processed.diff1{kk} MRSCont.processed.sum{kk}};
            basisSetConc = MRSCont.fit.basisSet;
            % Create the basis set with difference and sum basis functions
            basisSetConc.fids(:,:,1) = basisSetConc.fids(:,:,3);
            basisSetConc.specs(:,:,1) = basisSetConc.specs(:,:,3);
            basisSetConc.fids(:,:,2) = basisSetConc.fids(:,:,4);
            basisSetConc.specs(:,:,2) = basisSetConc.specs(:,:,4);
            basisSetConc.fids(:,:,3:4) = [];
            basisSetConc.specs(:,:,3:4) = [];
            
            if ~strcmp(fitOpts.coMM3, 'none')
                basisSetDiff1 = MRSCont.fit.basisSet;
                basisSetDiff1.fids = basisSetDiff1.fids(:,:,3);
                basisSetDiff1.specs = basisSetDiff1.specs(:,:,3);
                [basisSetDiff1] = osp_addDiffMMPeaks(basisSetDiff1,fitOpts);
                basisSetConc.fids(:,:,1) = basisSetDiff1.fids(:,:);
                basisSetConc.specs(:,:,1) = basisSetDiff1.specs(:,:);
            end
            
            % Apply scaling factor to the data
            for rr = 1:length(dataToFit)
                dataToFit{rr}   = op_ampScale(dataToFit{rr}, 1/MRSCont.fit.scale{kk});
            end
            % Call the multi-spectrum fit function
            [fitParamsConc, resBasisSetConc] = fit_runFitMultiplex(dataToFit, basisSetConc, fitModel, fitOpts);

            % Save back the basis set and fit parameters to MRSCont
            MRSCont.fit.resBasisSet.conc{kk}           = resBasisSetConc;
            MRSCont.fit.results.conc.fitParams{kk} = fitParamsConc;
        end
    end           
end
time = toc(metFitTime);
[~] = printLog('done',time,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
MRSCont.runtime.FitMet = time;

end