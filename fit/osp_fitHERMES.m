function [MRSCont] = osp_fitHERMES(MRSCont)
%% [MRSCont] = osp_fitHERMES(MRSCont)
%   This function performs spectral fitting of HERMES MRS data.
%
%   USAGE:
%       [MRSCont] = osp_fitHERMES(MRSCont);
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
    % For the separate (classic) HERMES fit, model the two DIFF
    % spectra and the SUM spectrum separately.
    if strcmp(fitStyle, 'Separate')
        if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.sum.fitParams))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)  

            %%% 2a. FIT SUM-SPECTRUM
            % Apply scaling factor to the data
            dataToFit   = MRSCont.processed.sum{kk};
            basisSetSum = MRSCont.fit.basisSet;
            basisSetSum.fids = basisSetSum.fids(:,:,7);
            basisSetSum.specs = basisSetSum.specs(:,:,7);
            dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
            
            % Call the fit function
            [fitParamsSum, resBasisSetSum] = fit_runFit(dataToFit, basisSetSum, fitModel, fitOpts);

            % Save back the basis set and fit parameters to MRSCont
            MRSCont.fit.resBasisSet.sum{kk}             = resBasisSetSum;
            MRSCont.fit.results.sum.fitParams{kk}   = fitParamsSum;


            %%% 2b. FIT DIFF1-SPECTRUM
            % Apply scaling factor to the data
            dataToFit   = MRSCont.processed.diff1{kk};
            basisSetDiff1 = MRSCont.fit.basisSet;
            basisSetDiff1.fids = basisSetDiff1.fids(:,:,5);
            basisSetDiff1.specs = basisSetDiff1.specs(:,:,5);
            dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
            dataToFit.refShift   = fitParamsSum.refShift;
            dataToFit.refFWHM   = fitParamsSum.refFWHM;
            
            if isfield(fitOpts, 'coMM3') && ~strcmp(fitOpts.coMM3, 'none')
                [basisSetDiff1] = osp_addDiffMMPeaks(basisSetDiff1,basisSetSum,fitOpts);
            end

            % Call the fit function
            [fitParamsDiff1, resBasisSetDiff1]  = fit_runFit(dataToFit, basisSetDiff1, fitModel, fitOpts);

            % Save back the basis set and fit parameters to MRSCont
            MRSCont.fit.resBasisSet.diff1{kk}           = resBasisSetDiff1;
            MRSCont.fit.results.diff1.fitParams{kk} = fitParamsDiff1;


            %%% 2c. FIT DIFF2-SPECTRUM
            % Apply scaling factor to the data
            dataToFit   = MRSCont.processed.diff2{kk};
            basisSetDiff2 = MRSCont.fit.basisSet;
            basisSetDiff2.fids = basisSetDiff2.fids(:,:,6);
            basisSetDiff2.specs = basisSetDiff2.specs(:,:,6);
            dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
            dataToFit.refShift   = fitParamsSum.refShift;
            dataToFit.refFWHM   = fitParamsSum.refFWHM;
            fitOpts.Diff2 = 1;

            % Call the fit function
            [fitParamsDiff2, resBasisSetDiff2]  = fit_runFit(dataToFit, basisSetDiff2, fitModel, fitOpts);

            % Save back the basis set and fit parameters to MRSCont
            MRSCont.fit.resBasisSet.diff2{kk}           = resBasisSetDiff2;
            MRSCont.fit.results.diff2.fitParams{kk} = fitParamsDiff2;
        end
    end


    %%% 3. CONCATENATED FIT %%%
    % For the concatenated MEGA fit, model the DIFF1 and SUM spectra
    % together.
    if strcmp(fitStyle, 'Concatenated')
        if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.conc.fitParams))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)   
        
            %%% 3a. FIT CONCATENATED SPECTRUM
            % Apply scaling factor to the data
            dataToFit   = {MRSCont.processed.diff1{kk}, MRSCont.processed.diff2{kk}, MRSCont.processed.sum{kk}};
            basisSetConc = MRSCont.fit.basisSet;
            basisSetConc.fids = basisSetConc.fids(:,:,5:7);
            basisSetConc.specs = basisSetConc.specs(:,:,5:7);
            for rr = 1:length(dataToFit)
                dataToFit{rr}   = op_ampScale(dataToFit{rr}, 1/MRSCont.fit.scale{kk});
            end

            % Call the fit function
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