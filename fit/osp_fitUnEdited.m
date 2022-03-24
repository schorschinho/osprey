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


    % ----- Osprey fit pipeline -----
    if strcmpi(MRSCont.opts.fit.method, 'Osprey')
        [~] = printLog('OspreyFit', kk, MRSCont.nDatasets, progressText, MRSCont.flags.isGUI, MRSCont.flags.isMRSI);

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
                MRSCont.fit.basisSet_mm;
            end                                         % re_mm
        end

        % ----- LCModel wrapper fit pipeline -----
    elseif strcmpi(MRSCont.opts.fit.method, 'LCModel')
        [~] = printLog('OspreyFit', kk, MRSCont.nDatasets, progressText, MRSCont.flags.isGUI, MRSCont.flags.isMRSI,1);


        % Put the scale to be 1 for now. Might have to adjust.
        MRSCont.fit.scale{kk} = 1;

        % Create the sub-folder where the LCModel results will be saved,
        % otherwise LCModel will throw 'FATAL ERROR MAIN 2'.
        if ~exist(fullfile(MRSCont.outputFolder, 'LCMoutput'), 'dir')
            mkdir(fullfile(MRSCont.outputFolder, 'LCMoutput'));
        end


        % Find the path to LCModel binary
        pathLCModelBinaryStruct = osp_platform('lcmodel');
        pathLCModelBinary = fullfile(MRSCont.ospFolder, 'libraries', 'LCModel',pathLCModelBinaryStruct.os,pathLCModelBinaryStruct.osver,['LCModel_' pathLCModelBinaryStruct.os,'_',pathLCModelBinaryStruct.osver]);
        switch  pathLCModelBinaryStruct.os
            case 'macos'
                bin = 'lcmodel';
            case 'unix'
                bin = 'lcmodel';
            case 'win'
                bin = 'LCModel.exe';
        end

        if ~exist([pathLCModelBinary filesep bin],'file') %Unzip binary
            unzip([pathLCModelBinary filesep bin '.zip'],[pathLCModelBinary filesep]);
        end

        if ~exist([pathLCModelBinary filesep bin],'file') %Binary is still missing
            error('ERROR: No LCModel binary found.  Aborting!!');
        else
            addpath(which('libraries/LCModel'));
        end

        callLCModel(MRSCont, MRSCont.opts.fit.lcmodel.controlfileA{kk},[pathLCModelBinary filesep bin]);

        % Save the parameters and information about the basis set
        MRSCont.fit.results.off.fitParams{kk} = readLCMFitParams(MRSCont, 'A', kk);

    end

end

time = toc(metFitTime);
[~] = printLog('done', time, MRSCont.nDatasets, progressText, MRSCont.flags.isGUI, MRSCont.flags.isMRSI);
MRSCont.runtime.FitMet = time;

end


function callLCModel(MRSCont, controlFile,pathLCModelBinary)
% Wrapper function for LCModel binary

callLCMCommand = ['"' pathLCModelBinary '" < "' controlFile '"'];
system(callLCMCommand);

end

function fitParams = readLCMFitParams(MRSCont, which, kk)

% Fall back to defaults if not provided
if nargin < 3
    which = 'A';
    if nargin < 2
        kk = 1;
        if nargin < 1
            error('ERROR: no input Osprey container specified.  Aborting!!');
        end
    end
end

% Determine output file strings
switch which
    case {'A','B','C','D'}
        lcmOutputFile = ['outputfile' which];
    case {'diff1'}
        lcmOutputFile = 'outputfileDiff1';
    case {'diff2'}
        lcmOutputFile = 'outputfileDiff2';
    case {'sum'}
        lcmOutputFile = 'outputfileSum';
end

% Read the fit results from the .table files
tab                 = mrs_readLcmodelTABLE([MRSCont.opts.fit.lcmodel.outputfileA{kk} '.table']);
fitParams.name      = tab.name;
fitParams.CRLB      = tab.SDpct;
fitParams.relConc   = tab.relative_conc;
fitParams.ph0       = tab.ph0;
fitParams.ph1       = tab.ph1;
fitParams.refShift  = tab.refShift;
fitParams.refFWHM   = tab.fwhm;
fitParams.SNR       = tab.snr;

%Remove the - in -CrCH2 because it interferes with the downstream functions
idx = find(strcmp(fitParams.name,'-CrCH2'));
if ~isempty(idx)
    fitParams.name{idx} = 'CrCH2';
end
%Remove the + in combinations because it interferes with the downstream functions
idx = find(contains(fitParams.name,'+'));
if ~isempty(idx)
    for combs = 1 : length(idx)
        fitParams.name{idx(combs)} = strrep(fitParams.name{idx(combs)},'+','_');
    end
end

% Read the spectrum, fit, and baseline from the .coord files
[ spectra, spectra_metabolites, x_ppm, info ] = mrs_readLcmodelCOORD( [MRSCont.opts.fit.lcmodel.(lcmOutputFile){kk} '.coord'] );
fitParams.ppm           = x_ppm;
fitParams.data          = spectra(:,1);
fitParams.completeFit   = spectra(:,2);
fitParams.baseline      = spectra(:,3);
fitParams.residual      = fitParams.data - fitParams.completeFit;

% The .coord files also contain the individual metabolite fits, BUT only if
% the estimate is not zero, and the individual metabolite fits include the
% baseline.
fitParams.indivMets     = zeros(info.n, length(fitParams.name));
for rr = 1:length(fitParams.name)
    % Check whether a particular metabolite has been fit
    idxMatch = find(strcmp(info.metabolites, fitParams.name{rr}));
    % If it has been fit, save the fit with the correct index after
    % subtracting the baseline
    if idxMatch
        fitParams.indivMets(:,rr) = spectra_metabolites(:,idxMatch) - fitParams.baseline;
    end
end
% Store amplitudes regardless whether the are fit or not
fitParams.ampl        = tab.concentration';

% Read the raw area of the unsuppressed water peak from the .print file
infoPrint = mrs_readLcmodelPRINT( [MRSCont.opts.fit.lcmodel.(lcmOutputFile){kk} '.print'] );
if isfield(infoPrint,'h2oarea')
    fitParams.h2oarea = infoPrint.h2oarea;
end



end
