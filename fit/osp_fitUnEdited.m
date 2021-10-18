function [MRSCont] = osp_fitUnEdited(MRSCont)
%% [MRSCont] = osp_fitUnEdited(MRSCont)
%   This function performs spectral fitting of unedited MRS data.
%
%   USAGE:
%       [MRSCont] = osp_fitUnEdited(MRSCont, basisSet);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       basisSet    = Osprey basis set container.
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
    
    if ~(MRSCont.flags.didFit == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'fit') && (kk > length(MRSCont.fit.results.off.fitParams))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
        % Apply scaling factor to the data
        dataToFit   = MRSCont.processed.A{kk};
        dataToFit   = op_ampScale(dataToFit, 1/MRSCont.fit.scale{kk});
        % Extract fit options
        fitOpts     = MRSCont.opts.fit;
        fitModel    = fitOpts.method;
        
        % If MRSI data load priors
        if MRSCont.flags.isMRSI
            if isfield(MRSCont.fit, 'results')
                fitParamsOff   = MRSCont.fit.results.off.fitParams{kk};
%                 dataToFit.refShift   = fitParamsOff.refShift;
%                 dataToFit.refFWHM   = fitParamsOff.refFWHM;
                fitOpts.MRSIpriors = fitParamsOff;
            end
        end

        % Call the fit function
        basisSet    = MRSCont.fit.basisSet;
        [fitParams, resBasisSet] = fit_runFit(dataToFit, basisSet, fitModel, fitOpts);

        % Save back the basis set and fit parameters to MRSCont
        MRSCont.fit.basisSet                    = basisSet;
        MRSCont.fit.resBasisSet.off{kk}         = resBasisSet;
        MRSCont.fit.results.off.fitParams{kk}   = fitParams;
    %Modeling MM spectra after the main spectrum re_mm
        if MRSCont.flags.hasMM == 1              %re_mm
            dataToFit_mm   = MRSCont.processed.mm{kk};
            dataToFit_mm   = op_ampScale(dataToFit_mm, 1/MRSCont.fit.scale{kk});
            %add some info from the metabolite fit
            dataToFit_mm.lineShape  = fitParams.lineShape;
            dataToFit_mm.refFWHM  = fitParams.refFWHM;
            % Extract fit options
            fitOpts_mm    = MRSCont.opts.fit;
            fitModel_mm    = fitOpts.method;
            fitOpts_mm.sequence = 'unedited';
            %Specify a reduced basis set for MM modeling
            %basisSet_mm    = MRSCont.fit.basisSet;
            %Reduce the size of the basis set
            
            %Adjust basis set
            % Clear existing basis set
            MRSCont.fit.basisSet_mm = [];
            % Load the specified basis set
            basisSet_mm = load(MRSCont.opts.fit.basisSetFile);
            basisSet_mm = basisSet_mm.BASIS;
            % Generate the list of basis functions that are supposed to be included in
            % the basis set
            % To do: Interface with interactive user input
            metabList_mm = fit_createMetabListMM('unedited');
            % Collect MMfit flag from the options determined in the job file
            fitMM = MRSCont.opts.fit.fitMM;
            % Create the modified basis set
            basisSet_mm = fit_selectMetabs(basisSet_mm, metabList_mm, fitMM);
             % Call the fit function
            [fitParams_mm, resBasisSet_mm] = fit_runFitMM(dataToFit_mm, basisSet_mm, fitModel_mm, fitOpts_mm);
             % Save back the basis set and fit parameters to MRSCont
            MRSCont.fit.basisSet_mm                    = basisSet_mm;
            MRSCont.fit.resBasisSet.mm{kk}         = resBasisSet_mm;
            MRSCont.fit.results.mm.fitParams{kk}   = fitParams_mm;
            MRSCont.fit.basisSet_mm
        end                                         % re_mm
    end

end
time = toc(metFitTime);
[~] = printLog('done',time,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 
MRSCont.runtime.FitMet = time;
end