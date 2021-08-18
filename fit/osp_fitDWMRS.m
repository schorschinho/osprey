function [MRSCont] = osp_fitDWMRS(MRSCont)
%% [MRSCont] = osp_fitDWMRS(MRSCont)
%   This function performs spectral fitting of diffusion-weighted MRS data.
%
%   USAGE:
%       [MRSCont] = osp_fitDWMRS(MRSCont);
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
    
    [~] = printLog('OspreyFit', kk, MRSCont.nDatasets, progressText, MRSCont.flags.isGUI, MRSCont.flags.isMRSI);
    
    % Loop over number of diffusion-weighted spectra and
    % create a control file for each spectrum that needs to
    % be fit
    nDW = length(MRSCont.opts.fit.lcmodel.outfileA{kk});
    for dd = 1:nDW
        
        % ----- Osprey fit pipeline -----
        if strcmpi(MRSCont.opts.fit.method, 'Osprey')
            
            error('Diffusion MRS data is currently modeled with LCModel.');
            
            % ----- LCModel wrapper fit pipeline -----
        elseif strcmpi(MRSCont.opts.fit.method, 'LCModel')
            
            % Put the scale to be 1 for now. Might have to adjust.
            MRSCont.fit.scale{kk} = 1;
            
            % Create the sub-folder where the LCModel results will be saved,
            % otherwise LCModel will throw 'FATAL ERROR MAIN 2'.
            if ~exist(fullfile(MRSCont.outputFolder, 'LCMoutput'), 'dir')
                mkdir(fullfile(MRSCont.outputFolder, 'LCMoutput'));
            end
            callLCModel(MRSCont, MRSCont.opts.fit.lcmodel.controlfileA{kk}{dd});
            
            % Save the parameters and information about the basis set
            MRSCont.fit.results.off.fitParams{kk}{dd} = readLCMFitParams(MRSCont, 'A', kk, dd);
            
        end
    end
    
end

time = toc(metFitTime);
[~] = printLog('done', time, MRSCont.nDatasets, progressText, MRSCont.flags.isGUI, MRSCont.flags.isMRSI);
MRSCont.runtime.FitMet = time;

end


function callLCModel(MRSCont, controlFile)
% Wrapper function for LCModel binary

callLCMCommand = ['"' fullfile(MRSCont.ospFolder, 'libraries', 'LCModel') filesep 'lcmodel" < "' controlFile '"'];
system(callLCMCommand);

end

function fitParams = readLCMFitParams(MRSCont, which, kk, dd)

% Fall back to defaults if not provided
if nargin < 4
    dd = 1;
    if nargin < 3
        kk = 1;
        if nargin < 2
            which = 'A';
            if nargin < 1
                error('ERROR: no input Osprey container specified.  Aborting!!');
            end
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
tab                 = mrs_readLcmodelTABLE([MRSCont.opts.fit.lcmodel.outputfileA{kk}{dd} '.table']);
fitParams.name      = tab.name;
fitParams.CRLB      = tab.SDpct;
fitParams.relConc   = tab.relative_conc;
fitParams.ph0       = tab.ph0;
fitParams.ph1       = tab.ph1;
fitParams.FWHM      = tab.fwhm;
fitParams.SNR       = tab.snr;

% Read the spectrum, fit, and baseline from the .coord files
[ spectra, spectra_metabolites, x_ppm, info ] = mrs_readLcmodelCOORD( [MRSCont.opts.fit.lcmodel.(lcmOutputFile){kk}{dd} '.coord'] );
fitParams.ppm           = x_ppm;
fitParams.data          = spectra(:,1);
fitParams.completeFit   = spectra(:,2);
fitParams.baseline      = spectra(:,3);
fitParams.residual      = fitParams.data - fitParams.completeFit;

% The .coord files also contain the individual metabolite fits, BUT only if
% the estimate is not zero, and the individual metabolite fits include the
% baseline.
fitParams.ampl          = tab.concentration;
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

% Read the raw area of the unsuppressed water peak from the .print file
infoPrint = mrs_readLcmodelPRINT( [MRSCont.opts.fit.lcmodel.(lcmOutputFile){kk}{dd} '.print'] );
if isfield(infoPrint, 'h2oarea')
    fitParams.h2oarea = infoPrint.h2oarea;
end

  
end

