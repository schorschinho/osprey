function [MRSCont] = LCG_processMEGA(MRSCont)
%% [MRSCont] = LCG_processMEGA(MRSCont)
%   This function performs the following steps to process MEGA-edited
%   (2-step) MRS data (e.g. MEGA-PRESS, MEGA-sLASER):
%       - Alignment of individual averages using spectral registration
%       - Averaging
%       - Removal of residual water using HSVD filtering
%       - Klose Eddy current correction (if a reference scan is provided)
%       - Automated zero-order phase correction
%       - Correct referencing of the ppm frequency axis
%
%   USAGE:
%       [MRSCont] = LCG_processMEGA(MRSCont);
%
%   INPUTS:
%       MRSCont     = LCGannet MRS data container.
%
%   OUTPUTS:
%       MRSCont     = LCGannet MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-22)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-02-22: First version of the code.

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
        
    % Siemens TWIX data sometimes comes out misphased - add a first-point phase
    % correction here
    if strcmpi(MRSCont.vendor, 'Siemens') && strcmpi(MRSCont.datatype, 'TWIX')
        corrph = conj(raw.fids(2,:))./abs(raw.fids(2,:));
        corrph = repmat(corrph, [raw.sz(1) 1]);
        if length(raw.sz) == 2
            corrph = reshape(corrph, [raw.sz(1) raw.sz(2)]);
        elseif length(raw.sz) == 3
            corrph = reshape(corrph, [raw.sz(1) raw.sz(2) raw.sz(3)]);
        end
        raw.fids    = raw.fids .* corrph;
        raw.specs   = fftshift(fft(raw.fids,[],raw.dims.t),raw.dims.t);
       
    end
    
    % Data may be mis-phased because of poor or incomplete
    % water suppression. Get an estimate for the right phase by phasing a
    % water-filtered OFF spectrum.
    raw_A                       = op_takesubspec(raw,1);                                    % Get first subspectrum
    raw_B                       = op_takesubspec(raw,2);                                    % Get second subspectrum

    % Perform frequency-and-phase correction. This can obviously only be done,
    % if the spectra have not been pre-averaged, i.e. in some older RDA and
    % DICOM files (which should, generally, not be used).
    if ~raw.flags.averaged
        raw                         = op_concatAverages(raw_A,raw_B);                           % Concatenate all averages together again
        raw                         = op_alignAverages(raw,0.2,'n');                              % Align averages
        raw_A                       = op_averaging(op_takeaverages(raw,1:size(raw.fids,2)/2));  %Take ONs and average
        raw_B                       = op_averaging(op_takeaverages(raw,(size(raw.fids,2)/2+1):(size(raw.fids,2)))); %Take OFFs and average
    end
        
    
    %%% 2. GET REFERENCE DATA / EDDY CURRENT CORRECTION %%%
    % If there are reference scans, perform the same operations
    if MRSCont.flags.hasRef
        raw_ref                         = MRSCont.raw_ref{kk};              % Get the kk-th dataset
        
        % Siemens TWIX data sometimes come out misphased - add a first-point phase
        % correction here
        if strcmpi(MRSCont.vendor, 'Siemens') && strcmpi(MRSCont.datatype, 'TWIX')
            corrph = conj(raw_ref.fids(2,:))./abs(raw_ref.fids(2,:));
            corrph = repmat(corrph, [raw_ref.sz(1) 1]);
            if length(raw_ref.sz) == 2
                corrph = reshape(corrph, [raw_ref.sz(1) raw_ref.sz(2)]);
            elseif length(raw_ref.sz) == 3
                corrph = reshape(corrph, [raw_ref.sz(1) raw_ref.sz(2) raw_ref.sz(3)]);
            end
            raw_ref.fids = raw_ref.fids .* corrph;
            raw_ref.specs=fftshift(fft(raw_ref.fids,[],raw_ref.dims.t),raw_ref.dims.t);
        end
        
        % Some formats end up having subspectra in their reference scans
        % (e.g. Philips), as well as empty lines. Intercept these cases
        % here.
        if raw_ref.subspecs > 1
            raw_ref_A                       = op_takesubspec(raw_ref,1);
            [raw_ref_A]                     = op_rmempty(raw_ref_A);            % Remove empty lines
            raw_ref_B                       = op_takesubspec(raw_ref,2);
            [raw_ref_B]                     = op_rmempty(raw_ref_B);            % Remove empty lines
            raw_ref                         = op_concatAverages(raw_ref_A,raw_ref_B);
        end
        if ~raw_ref.flags.averaged
            [raw_ref,fs_A,phs_A]            = op_alignAverages(raw_ref,1,'n');  % Align averages
            raw_ref                         = op_averaging(raw_ref);            % Average
        end
        
%         % The following only for Big GABA dataset - a few Siemens datasets
%         % have been accidentally acquired with the water suppression switched
%         % on for the water reference scan. In that case, don't do ECC, but
%         % rather leave it to the phase correction in the next step.
%         if strcmp(MRSCont.vendor,'Siemens') && kk >= 37 && kk <= 42
%         else
            [raw_A,~]                   = op_eccKlose(raw_A, raw_ref);
            [raw_B,raw_ref]             = op_eccKlose(raw_B, raw_ref);        % Klose eddy current correction
%         end
        
        [raw_ref,~]                     = op_ppmref(raw_ref,4.6,4.8,4.68);  % Reference to water @ 4.68 ppm
        MRSCont.processed.ref{kk}       = raw_ref;                          % Save back to MRSCont container
    end
    
    
    %%% 3. DETERMINE POLARITY OF SPECTRUM (EG FOR MOIST WATER SUPP) %%%
    % Automate determination whether the Cr peak has positive polarity.
    % For water suppression methods like MOIST, the residual water may
    % actually have negative polarity, but end up positive in the data, so
    % that the spectrum needs to be flipped.
    raw_A_Cr    = op_freqrange(raw_A,2.8,3.2);
    % Determine the polarity of the respective peak: if the absolute of the
    % maximum minus the absolute of the minimum is positive, the polarity 
    % of the respective peak is positive; if the absolute of the maximum 
    % minus the absolute of the minimum is negative, the polarity is negative.
    polResidCr = abs(max(real(raw_A_Cr.specs))) - abs(min(real(raw_A_Cr.specs)));
    if polResidCr < 0
        raw_A = op_ampScale(raw_A,-1);
        raw_B = op_ampScale(raw_B,-1);
    end


    %%% 4. DETERMINE ON/OFF STATUS
    % Classify the two sub-spectra such that the OFF spectrum is stored to
    % field A, and the ON spectrum is stored to field B.
    [raw_A, raw_B]              = LCG_onOffClassifyMEGA(raw_A, raw_B, 'GABA');
    
    
    %%% 5. REMOVE RESIDUAL WATER %%%
    % Remove water and correct back to baseline.
    % The spectra sometimes become NaNs after filtering with too many
    % components. Loop over decreasing numbers of components here.
    [raw_A_temp,K_A,amp_A]           = op_removeWater(raw_A,[4.6 4.8],20,0.75*length(raw_A.fids),0); % Remove the residual water
    if isnan(real(raw_A_temp.fids))
        rr = 30;
        while isnan(real(raw_A_temp.fids))
            [raw_A_temp,K_A,amp_A]           = op_removeWater(raw_A,[4.6 4.8],rr,0.75*length(raw_A.fids),0); % Remove the residual water
            rr = rr-1;
        end
    end
    raw_A = raw_A_temp;
    raw_A                       = op_fddccorr(raw_A,100);                                 % Correct back to baseline

    [raw_B_temp,K_B,amp_B]           = op_removeWater(raw_B,[4.6 4.8],20,0.75*length(raw_B.fids),0); % Remove the residual water
    if isnan(real(raw_B_temp.fids))
       rr = 30;
       while isnan(real(raw_B_temp.fids))
           [raw_B_temp,K_B,amp_B]           = op_removeWater(raw_B,[4.6 4.8],rr,0.75*length(raw_B.fids),0); % Remove the residual water
            rr = rr-1;
       end
    end
    raw_B = raw_B_temp;
    raw_B                       = op_fddccorr(raw_B,100);                                 % Correct back to baseline

    
    %%% 6. BUILD SUM AND DIFF SPECTRA %%%
    % Correct the frequency axis so that Cr appears at 3.027 ppm
    [raw_A_test,f0] = op_ppmref(raw_A,2.9,3.1,3.027);
    % Fit a double-Lorentzian to the Cr-Cho area, and phase the spectrum
    % with the negative phase of that fit
    [raw_A_test,phShift] = op_phaseCrCho(raw_A_test, 1);
    % Align the ON to the OFF spectrum
    [raw_B_test,~,~] = op_alignScans(raw_B,raw_A_test,0.2,'fp');
    raw_A = raw_A_test;
    raw_B = raw_B_test;
    sum                         = op_addScans(raw_A,raw_B);
    diff1                       = op_addScans(raw_B,raw_A,1);
    [diff1_temp,K_diff1,amp_diff1]           = op_removeWater(diff1,[4.6 4.8],20,0.75*length(diff1.fids),0); % Remove the residual water
    if isnan(real(diff1_temp.fids))
        rr = 30;
        while isnan(real(diff1_temp.fids))
            [diff1_temp,K_diff1,amp_diff1]           = op_removeWater(diff1,[4.6 4.8],rr,0.75*length(diff1.fids),0); % Remove the residual water
            rr = rr-1;
        end
    end
    diff1 = diff1_temp;
    diff1                       = op_fddccorr(diff1,100);                                 % Correct back to baseline
    [sum_temp,K_sum,amp_sum]           = op_removeWater(sum,[4.6 4.8],20,0.75*length(sum.fids),0); % Remove the residual water
    if isnan(real(sum_temp.fids))
        rr = 30;
        while isnan(real(sum_temp.fids))
            [sum_temp,K_sum,amp_sum]           = op_removeWater(sum,[4.6 4.8],rr,0.75*length(sum.fids),0); % Remove the residual water
            rr = rr-1;
        end
    end
    sum = sum_temp;
    sum                         = op_fddccorr(sum,100);
    
    
    %%% 7. REFERENCE SPECTRUM CORRECTLY TO FREQUENCY AXIS 
    % Reference resulting data correctly and consistently
    [raw_A,ref_shift]           = op_ppmref(raw_A,1.8,2.2,2.008);           % Reference edit-OFF spectrum to NAA @ 2.008 ppm                                                                          
    [raw_B]                     = op_freqshift(raw_B,ref_shift);            % Apply same shift to edit-OFF
    [diff1]                     = op_freqshift(diff1,ref_shift);            % Apply same shift to diff1
    [sum]                       = op_freqshift(sum,ref_shift);              % Apply same shift to sum
    
    
    %%% 8. SAVE BACK TO MRSCONT CONTAINER
    MRSCont.processed.A{kk}     = raw_A;                                    % Save edit-OFF back to MRSCont container
    MRSCont.processed.B{kk}     = raw_B;                                    % Save edit-ON back to MRSCont container
    MRSCont.processed.diff1{kk} = diff1;                                    % Save diff1 back to MRSCont container
    MRSCont.processed.sum{kk}   = sum;                                      % Save sum back to MRSCont container


    %%% 9. GET SHORT-TE WATER DATA %%%
    if MRSCont.flags.hasWater
        % Some formats end up having subspectra in their reference scans
        % (e.g. Philips), as well as empty lines. Intercept these cases
        % here.
        raw_w                           = MRSCont.raw_w{kk};                % Get the kk-th dataset
        if raw_w.subspecs > 1
            raw_w_A                         = op_takesubspec(raw_w,1);
            [raw_w_A]                       = op_rmempty(raw_w_A);              % Remove empty lines
            raw_w_B                         = op_takesubspec(raw_w,2);
            [raw_w_A]                       = op_rmempty(raw_w_A);              % Remove empty lines
            raw_w                           = op_concatAverages(raw_w_A,raw_w_B);
        end
        if ~raw_w.flags.averaged
            [raw_w,fs_A,phs_A]              = op_alignAverages(raw_w,1,'n');    % Align averages
            raw_w                           = op_averaging(raw_w);              % Average
        end
        [raw_w,~]                       = op_eccKlose(raw_w, raw_w);        % Klose eddy current correction
        [raw_w,~]                       = op_ppmref(raw_w,4.6,4.8,4.68);    % Reference to water @ 4.68 ppm
        MRSCont.processed.w{kk}         = raw_w;                            % Save back to MRSCont container
    end
    
    
    %%% 10. QUALITY CONTROL PARAMETERS %%%
    % Calculate some spectral quality metrics here;
    MRSCont.QM.SNR.A{kk}    = op_getSNR(MRSCont.processed.A{kk}); % NAA amplitude over noise floor
    FWHM_Hz                 = op_getLW(MRSCont.processed.A{kk},1.8,2.2); % in Hz
    MRSCont.QM.FWHM.A{kk}   = FWHM_Hz./MRSCont.processed.A{kk}.txfrq*1e6; % convert to ppm
    
    if MRSCont.flags.hasRef
        MRSCont.QM.SNR.ref{kk}  = op_getSNR(MRSCont.processed.ref{kk},4.2,5.2); % water amplitude over noise floor
        FWHM_Hz                 = op_getLW(MRSCont.processed.ref{kk},4.2,5.2); % in Hz
        MRSCont.QM.FWHM.ref{kk} = FWHM_Hz./MRSCont.processed.ref{kk}.txfrq*1e6; % convert to ppm
    end
    if MRSCont.flags.hasWater
        MRSCont.QM.SNR.w{kk}    = op_getSNR(MRSCont.processed.w{kk},4.2,5.2); % water amplitude over noise floor
        FWHM_Hz                 = op_getLW(MRSCont.processed.w{kk},4.2,5.2); % in Hz
        MRSCont.QM.FWHM.w{kk}   = FWHM_Hz./MRSCont.processed.w{kk}.txfrq*1e6; % convert to ppm
    end
            
end
fprintf('... done.\n');
toc(refProcessTime);


%%% 11. SET FLAGS %%%
MRSCont.flags.avgsAligned       = 1;
MRSCont.flags.averaged          = 1;
MRSCont.flags.ECCed             = 1;
MRSCont.flags.waterRemoved      = 1;

% Close any remaining open figures
close all;

end