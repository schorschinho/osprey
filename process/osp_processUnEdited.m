function [MRSCont] = osp_processUnEdited(MRSCont)
%% [MRSCont] = osp_processUnEdited(MRSCont)
%   This function performs the following steps to process un-edited MRS
%   data (e.g. PRESS, STEAM, sLASER):
%       - Alignment of individual averages using robust spectral registration
%       - Averaging
%       - Removal of residual water using HSVD filtering
%       - Klose Eddy current correction (if a reference scan is provided)
%       - Correct referencing of the ppm frequency axis
%
%   USAGE:
%       [MRSCont] = osp_processUnEdited(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
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
            [raw_ref,fs,phs]               = op_alignAverages(raw_ref, 1, 'n');
            raw_ref                     = op_averaging(raw_ref);            % Average
            raw_ref.specReg.fs              = fs; % save align parameters
            raw_ref.specReg.phs             = phs; % save align parameters            
        else
            raw_ref.flags.averaged  = 1;
            raw_ref.dims.averages   = 0;
            raw_ref.specReg.fs = 0;
            raw_ref.specReg.phs = 0;               
        end
        [raw,raw_ref]                   = op_eccKlose(raw, raw_ref);        % Klose eddy current correction
        [raw_ref,~]                     = op_ppmref(raw_ref,4.6,4.8,4.68);  % Reference to water @ 4.68 ppm
        MRSCont.processed.ref{kk}       = raw_ref;                          % Save back to MRSCont container
    end
    
    
    %%% 3. FREQUENCY/PHASE CORRECTION AND AVERAGING %%%
    % Measure drift pre-alignment
    driftPre = op_measureDrift(raw);
    if raw.averages > 1 && raw.flags.averaged == 0
        [raw, fs, phs, weights]     = op_robustSpecReg(raw, 'unedited', 0); % Align and average
        raw.specReg.fs              = fs; % save align parameters
        raw.specReg.phs             = phs; % save align parameters
        raw.specReg.weights         = weights; % save align parameters
    else
        raw.flags.averaged  = 1;
        raw.dims.averages   = 0;
        raw.specReg.fs              = 0; % save align parameters
        raw.specReg.phs             = 0; % save align parameters
        raw.specReg.weights         = 1; % save align parameters
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
    [raw_temp,~,~]   = op_removeWater(raw,[4.6 4.8],20,0.75*length(raw.fids),0); % Remove the residual water
    if isnan(real(raw_temp.fids))
        rr = 30;
        while isnan(real(raw_temp.fids))
            [raw_temp,~,~]   = op_removeWater(raw,[4.6 4.8],rr,0.75*length(raw.fids),0); % Remove the residual water
            rr = rr-1;
        end
    end
    raw     = raw_temp;
    raw     = op_fddccorr(raw,100);                                     % Correct back to baseline
    
    
    %%% 6. REFERENCE SPECTRUM CORRECTLY TO FREQUENCY AXIS AND PHASE SIEMENS
    %%% DATA
    [raw, refShift]             = op_ppmref(raw,1.9,2.1,2.008);             % Reference to NAA @ 2.008 ppm
                                      
    % Save back to MRSCont container
    if strcmp(MRSCont.vendor,'Siemens')
        % Fit a double-Lorentzian to the Cr-Cho area, and phase the spectrum
        % with the negative phase of that fit
        [raw,~]       = op_phaseCrCho(raw, 1);
    end
    MRSCont.processed.A{kk}     = raw;   
    
    
    %%% 7. GET SHORT-TE WATER DATA %%%
    if MRSCont.flags.hasWater
        raw_w                           = MRSCont.raw_w{kk};                % Get the kk-th dataset
        if raw_w.averages > 1 && raw_w.flags.averaged == 0
            [raw_w,fs,phs]                 = op_alignAverages(raw_w, 1, 'n');
            raw_w.specReg.fs              = fs; % save align parameters
            raw_w.specReg.phs             = phs; % save align parameters
            raw_w                       = op_averaging(raw_w);              % Average
        else
            raw_w.flags.averaged    = 1;
            raw_w.dims.averages     = 0;
            raw_w.specReg.fs = 0;
            raw_w.specReg.phs = 0;   
        end
        [raw_w,~]                       = op_eccKlose(raw_w, raw_w);        % Klose eddy current correction
        [raw_w,~]                       = op_ppmref(raw_w,4.6,4.8,4.68);    % Reference to water @ 4.68 ppm
        MRSCont.processed.w{kk}         = raw_w; % Save back to MRSCont container
    end
    
    
    %%% 8. QUALITY CONTROL PARAMETERS %%%
    % Calculate some spectral quality metrics here;
    MRSCont.QM.SNR.A(kk)    = op_getSNR(MRSCont.processed.A{kk}); % NAA amplitude over noise floor
    FWHM_Hz                 = op_getLW(MRSCont.processed.A{kk},1.8,2.2); % in Hz
    MRSCont.QM.FWHM.A(kk)   = FWHM_Hz./MRSCont.processed.A{kk}.txfrq*1e6; % convert to ppm
    MRSCont.QM.drift.A{kk}= driftPre;
    MRSCont.QM.freqShift.A(kk) = refShift;
    
    if MRSCont.flags.hasRef
        MRSCont.QM.SNR.ref(kk)  = op_getSNR(MRSCont.processed.ref{kk},4.2,5.2); % water amplitude over noise floor
        FWHM_Hz                 = op_getLW(MRSCont.processed.ref{kk},4.2,5.2); % in Hz
        MRSCont.QM.FWHM.ref(kk) = FWHM_Hz./MRSCont.processed.ref{kk}.txfrq*1e6; % convert to ppm
    end
    if MRSCont.flags.hasWater
        MRSCont.QM.SNR.w(kk)    = op_getSNR(MRSCont.processed.w{kk},4.2,5.2); % water amplitude over noise floor
        FWHM_Hz                 = op_getLW(MRSCont.processed.w{kk},4.2,5.2); % in Hz
        MRSCont.QM.FWHM.w(kk)   = FWHM_Hz./MRSCont.processed.w{kk}.txfrq*1e6; % convert to ppm
    end
    
end
fprintf('... done.\n');
toc(refProcessTime);


%%% 9. SET FLAGS %%%
MRSCont.flags.avgsAligned       = 1;
MRSCont.flags.averaged          = 1;
MRSCont.flags.ECCed             = 1;
MRSCont.flags.waterRemoved      = 1;

% Close any remaining open figures
close all;

end
