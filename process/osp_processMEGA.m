function [MRSCont] = osp_processMEGA(MRSCont, target)
%% [MRSCont] = osp_processMEGA(MRSCont, target)
%   This function performs the following steps to process MEGA-edited
%   (2-step) MRS data (e.g. MEGA-PRESS, MEGA-sLASER):
%       - Alignment of individual averages using robust spectral registration
%       - Averaging
%       - Removal of residual water using HSVD filtering
%       - Klose Eddy current correction (if a reference scan is provided)
%       - Automated zero-order phase correction
%       - Correct referencing of the ppm frequency axis
%
%   USAGE:
%       [MRSCont] = osp_processMEGA(MRSCont, 'target');
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       target      = String. Can be 'GABA' or 'GSH'. Default: 'GABA'
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
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
%       2019-08-20: First public version of the code.

% Close any remaining open figures
close all;
warning('off','all');

% Parse input arguments
if nargin < 2
    target = 'GABA'; % GABA editing as default
end

%% Loop over all datasets
refProcessTime = tic;
reverseStr = '';
for kk = 1:MRSCont.nDatasets
    msg = sprintf('Processing data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    %%% 1. GET RAW DATA %%%
    raw         = MRSCont.raw{kk};                                          % Get the kk-th dataset
    
    % Get sub-spectra, depending on whether they are stored as such
    if raw.subspecs == 2
        raw_A   = op_takesubspec(raw,1);                    % Get first subspectrum
        raw_B   = op_takesubspec(raw,2);                    % Get second subspectrum
    else
        raw_A   = op_takeaverages(raw,1:2:raw.averages);    % Get first subspectrum
        raw_B   = op_takeaverages(raw,2:2:raw.averages);    % Get second subspectrum
    end

    % Perform robust spectral correction with weighted averaging.
    % This can obviously only be done, if the spectra have not been 
    % pre-averaged, i.e. in some older RDA and DICOM files (which should, 
    % generally, not be used).
    if ~raw.flags.averaged
        raw_A   = op_robustSpecReg(raw_A, 'MEGA', 0);
        raw_B   = op_robustSpecReg(raw_B, 'MEGA', 0);                  
    end
        
    
    %%% 2. GET REFERENCE DATA / EDDY CURRENT CORRECTION %%%
    % If there are reference scans, perform the same operations
    if MRSCont.flags.hasRef
        raw_ref                     = MRSCont.raw_ref{kk};              % Get the kk-th dataset
        
        % Some formats end up having subspectra in their reference scans
        % (e.g. Philips), as well as empty lines. Intercept these cases
        % here.
        if raw_ref.subspecs > 1
            raw_ref_A               = op_takesubspec(raw_ref,1);
            [raw_ref_A]             = op_rmempty(raw_ref_A);            % Remove empty lines
            raw_ref_B               = op_takesubspec(raw_ref,2);
            [raw_ref_B]             = op_rmempty(raw_ref_B);            % Remove empty lines
            raw_ref                 = op_concatAverages(raw_ref_A,raw_ref_B);
        end
        if ~raw_ref.flags.averaged
            [raw_ref,~,~]           = op_alignAverages(raw_ref,1,'n');  % Align averages
            raw_ref                 = op_averaging(raw_ref);            % Average
        end
        
%         % The following IF only for Big GABA dataset - a few Siemens datasets
%         % have been accidentally acquired with the water suppression switched
%         % on for the water reference scan. In that case, don't do ECC, but
%         % rather leave it to the phase correction in step 5.
%         if strcmp(MRSCont.vendor,'Siemens') && kk >= 37 && kk <= 42
%         else
            [raw_A,~]               = op_eccKlose(raw_A, raw_ref);
            [raw_B,raw_ref]         = op_eccKlose(raw_B, raw_ref);        % Klose eddy current correction
%         end
        
        [raw_ref,~]                 = op_ppmref(raw_ref,4.6,4.8,4.68);  % Reference to water @ 4.68 ppm
        MRSCont.processed.ref{kk}   = raw_ref;                          % Save back to MRSCont container
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
    polResidCr  = abs(max(real(raw_A_Cr.specs))) - abs(min(real(raw_A_Cr.specs)));
    if polResidCr < 0
        raw_A = op_ampScale(raw_A,-1);
        raw_B = op_ampScale(raw_B,-1);
    end


    %%% 4. DETERMINE ON/OFF STATUS
    % Classify the two sub-spectra such that the OFF spectrum is stored to
    % field A, and the ON spectrum is stored to field B.
    [raw_A, raw_B]  = osp_onOffClassifyMEGA(raw_A, raw_B, target);
    
    
    %%% 5. BUILD SUM AND DIFF SPECTRA %%%
    % Correct the frequency axis so that Cr appears at 3.027 ppm
    [raw_A,~]       = op_ppmref(raw_A,2.9,3.1,3.027);
    % Fit a double-Lorentzian to the Cr-Cho area, and phase the spectrum
    % with the negative phase of that fit
    [raw_A,~]       = op_phaseCrCho(raw_A, 1);
    % Align the sub-spectra to one another by minimizing the difference
    % between the common 'reporter' signals.
    [raw_A, raw_B]  = osp_editSubSpecAlign(raw_A, raw_B, target);
    % Create the sum spectrum
    sum             = op_addScans(raw_A,raw_B);
    % Create the GABA-edited difference spectrum
    diff1           = op_addScans(raw_B,raw_A,1);
    
    
    %%% 6. REMOVE RESIDUAL WATER %%%
    % Remove water and correct back to baseline.
    % The spectra sometimes become NaNs after filtering with too many
    % components. Loop over decreasing numbers of components here.
    [raw_A_temp,~,~]           = op_removeWater(raw_A,[4.6 4.8],20,0.75*length(raw_A.fids),0); % Remove the residual water
    if isnan(real(raw_A_temp.fids))
        rr = 30;
        while isnan(real(raw_A_temp.fids))
            [raw_A_temp,~,~]   = op_removeWater(raw_A,[4.6 4.8],rr,0.75*length(raw_A.fids),0); % Remove the residual water
            rr = rr-1;
        end
    end
    raw_A   = raw_A_temp;
    raw_A   = op_fddccorr(raw_A,100);                                 % Correct back to baseline

    [raw_B_temp,~,~]           = op_removeWater(raw_B,[4.6 4.8],20,0.75*length(raw_B.fids),0); % Remove the residual water
    if isnan(real(raw_B_temp.fids))
       rr = 30;
       while isnan(real(raw_B_temp.fids))
           [raw_B_temp,~,~]    = op_removeWater(raw_B,[4.6 4.8],rr,0.75*length(raw_B.fids),0); % Remove the residual water
            rr = rr-1;
       end
    end
    raw_B   = raw_B_temp;
    raw_B   = op_fddccorr(raw_B,100);                                 % Correct back to baseline

    [diff1_temp,~,~]           = op_removeWater(diff1,[4.6 4.8],20,0.75*length(diff1.fids),0); % Remove the residual water
    if isnan(real(diff1_temp.fids))
        rr = 30;
        while isnan(real(diff1_temp.fids))
            [diff1_temp,~,~]   = op_removeWater(diff1,[4.6 4.8],rr,0.75*length(diff1.fids),0); % Remove the residual water
            rr = rr-1;
        end
    end
    diff1   = diff1_temp;
    diff1   = op_fddccorr(diff1,100);                                 % Correct back to baseline
    
    [sum_temp,~,~]           = op_removeWater(sum,[4.6 4.8],20,0.75*length(sum.fids),0); % Remove the residual water
    if isnan(real(sum_temp.fids))
        rr = 30;
        while isnan(real(sum_temp.fids))
            [sum_temp,~,~]   = op_removeWater(sum,[4.6 4.8],rr,0.75*length(sum.fids),0); % Remove the residual water
            rr = rr-1;
        end
    end
    sum     = sum_temp;
    sum     = op_fddccorr(sum,100);
    
    
    %%% 7. REFERENCE SPECTRUM CORRECTLY TO FREQUENCY AXIS 
    % Reference resulting data correctly and consistently
    [raw_A,ref_shift]   = op_ppmref(raw_A,1.8,2.2,2.008);           % Reference edit-OFF spectrum to NAA @ 2.008 ppm                                                                          
    [raw_B]             = op_freqshift(raw_B,ref_shift);            % Apply same shift to edit-OFF
    [diff1]             = op_freqshift(diff1,ref_shift);            % Apply same shift to diff1
    [sum]               = op_freqshift(sum,ref_shift);              % Apply same shift to sum
    
    
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
        raw_w                       = MRSCont.raw_w{kk};                % Get the kk-th dataset
        if raw_w.subspecs > 1
            raw_w_A                 = op_takesubspec(raw_w,1);
            [raw_w_A]               = op_rmempty(raw_w_A);              % Remove empty lines
            raw_w_B                 = op_takesubspec(raw_w,2);
            [raw_w_A]               = op_rmempty(raw_w_A);              % Remove empty lines
            raw_w                   = op_concatAverages(raw_w_A,raw_w_B);
        end
        if ~raw_w.flags.averaged
            [raw_w,~,~]             = op_alignAverages(raw_w,1,'n');    % Align averages
            raw_w                   = op_averaging(raw_w);              % Average
        end
        [raw_w,~]                   = op_eccKlose(raw_w, raw_w);        % Klose eddy current correction
        [raw_w,~]                   = op_ppmref(raw_w,4.6,4.8,4.68);    % Reference to water @ 4.68 ppm
        MRSCont.processed.w{kk}     = raw_w;                            % Save back to MRSCont container
    end
    
    
    %%% 10. QUALITY CONTROL PARAMETERS %%%
    % Calculate some spectral quality metrics here;
    MRSCont.QM.SNR.A(kk)    = op_getSNR(MRSCont.processed.A{kk}); % NAA amplitude over noise floor
    FWHM_Hz                 = op_getLW(MRSCont.processed.A{kk},1.8,2.2); % in Hz
    MRSCont.QM.FWHM.A(kk)   = FWHM_Hz./MRSCont.processed.A{kk}.txfrq*1e6; % convert to ppm
    
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


%%% 11. SET FLAGS %%%
MRSCont.flags.avgsAligned   = 1;
MRSCont.flags.averaged      = 1;
MRSCont.flags.ECCed         = 1;
MRSCont.flags.waterRemoved  = 1;

% Close any remaining open figures
close all;

end