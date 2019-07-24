function [MRSCont] = LCG_processUnEdited(MRSCont)
%% [MRSCont] = LCG_processUnEdited(MRSCont)
%   This function performs the following steps to process un-edited MRS
%   data (e.g. PRESS, STEAM, sLASER):
%       - Removal of bad averages (determine by likeliness metric)
%       - Alignment of individual averages using spectral registration
%       - Averaging
%       - Removal of residual water using HSVD filtering
%       - Klose Eddy current correction (if a reference scan is provided)
%       - Correct referencing of the ppm frequency axis
%
%   USAGE:
%       [MRSCont] = LCG_processUnEdited(MRSCont);
%
%   INPUTS:
%       MRSCont     = LCGannet MRS data container.
%
%   OUTPUTS:
%       MRSCont     = LCGannet MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-20)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-02-20: First version of the code.

% Close any remaining open figures
close all;
warning('off','all');

%% Loop over all datasets
refProcessTime = tic;
reverseStr = '';
for kk = 1:MRSCont.nDatasets
    msg = sprintf('Processing data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    %%% 1. GET RAW DATA %%%
    raw                         = MRSCont.raw{kk};                                          % Get the kk-th dataset
    
    %%% 2. GET REFERENCE DATA / EDDY CURRENT CORRECTION %%%
    % If there are reference scans, load them here to allow eddy-current
    % correction of the raw data.
    if MRSCont.flags.hasRef
        raw_ref                         = MRSCont.raw_ref{kk};              % Get the kk-th dataset
        if raw_ref.averages > 1 && raw_ref.flags.averaged == 0
            [raw_ref,fs_ref,phs_ref]    = op_robustSpecReg(raw_ref, 0);
            raw_ref                     = op_averaging(raw_ref);            % Average
        else
            raw_ref.flags.averaged  = 1;
            raw_ref.dims.averages   = 0;
        end
        [raw,raw_ref]                   = op_eccKlose(raw, raw_ref);        % Klose eddy current correction
        [raw_ref,~]                     = op_ppmref(raw_ref,4.6,4.8,4.68);  % Reference to water @ 4.68 ppm
        MRSCont.processed.ref{kk}       = raw_ref;                          % Save back to MRSCont container
    end
    
    %%% 3. FREQUENCY/PHASE CORRECTION AND AVERAGING %%%
    if raw.averages > 1 && raw.flags.averaged == 0
        [raw,fs,phs]                = op_robustSpecReg(raw, 0);
        raw                         = op_averaging(raw);                                        % Average
    else
        raw.flags.averaged  = 1;
        raw.dims.averages   = 0;
    end
    
    %%% 4. DETERMINE POLARITY OF SPECTRUM (EG FOR MOIST WATER SUPP) %%%
    % Automate determination whether the NAA peak has positive polarity.
    % For water suppression methods like MOIST, the residual water may
    % actually have negative polarity, but end up positive in the data, so
    % that the spectrum needs to be flipped.
    raw_NAA     = op_freqrange(raw,1.9,2.1);
    % Determine the polarity of the respective peak: if the absolute of the
    % maximum minus the absolute of the minimum is positive, the polarity 
    % of the respective peak is positive; if the absolute of the maximum 
    % minus the absolute of the minimum is negative, the polarity is negative.
    polResidNAA = abs(max(real(raw_NAA.specs))) - abs(min(real(raw_NAA.specs)));
    if polResidNAA < 0
        raw = op_ampScale(raw,-1);
    end
    
    %%% 5. REMOVE RESIDUAL WATER %%%
    [raw,K_A,amp_A]             = op_removeWater(raw,[4.6 4.8],16,0.5*length(raw.fids),0); % Remove the residual water
    raw                         = op_fddccorr(raw,100);                                     % Correct back to baseline
    
    %%% 6. REFERENCE SPECTRUM CORRECTLY TO FREQUENCY AXIS
    [raw,~]                     = op_ppmref(raw,1.9,2.1,2.008);           % Reference to NAA @ 2.008 ppm
    MRSCont.processed.A{kk}     = raw;                                      % Save back to MRSCont container
    
    %%% 7. GET SHORT-TE WATER DATA %%%
    if MRSCont.flags.hasWater
        raw_w                           = MRSCont.raw_w{kk};                % Get the kk-th dataset
        if raw_w.averages > 1 && raw_w.flags.averaged == 0
            [raw_w,fs_w,phs_w]              = op_robustSpecReg(raw_w, 0);
            raw_w                           = op_averaging(raw_w);              % Average
        else
            raw_w.flags.averaged    = 1;
            raw_w.dims.averages     = 0;
        end
        [raw_w,~]                       = op_eccKlose(raw_w, raw_w);        % Klose eddy current correction
        [raw_w,~]                       = op_ppmref(raw_w,4.6,4.8,4.68);    % Reference to water @ 4.68 ppm
        MRSCont.processed.w{kk}         = raw_w; % Save back to MRSCont container
    end
    
    %%% 8. QUALITY CONTROL PARAMETERS %%%
    % Calculate some spectral quality metrics here;
    MRSCont.QM.SNR.A{kk}    = op_getSNR(MRSCont.processed.A{kk}); % NAA amplitude over noise floor
    FWHM_Hz                 = op_getLW(MRSCont.processed.A{kk},1.8,2.2); % in Hz
    MRSCont.QM.FWHM.A{kk}   = FWHM_Hz./MRSCont.processed.A{kk}.txfrq*1e6; % convert to ppm
   
end
fprintf('... done.\n');
toc(refProcessTime);

%%% 10. SET FLAGS %%%
MRSCont.flags.badAvgsRemoved    = 1;
MRSCont.flags.avgsAligned       = 1;
MRSCont.flags.averaged          = 1;
MRSCont.flags.ECCed             = 1;
MRSCont.flags.waterRemoved      = 1;

% Close any remaining open figures
close all;

end