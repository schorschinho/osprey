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

warning('off','all');

%% Loop over all datasets
refProcessTime = tic;
reverseStr = '';
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
    progressbar = waitbar(0,'Start','Name','Osprey Process');
    waitbar(0,progressbar,sprintf('Processed data from dataset %d out of %d total datasets...\n', 0, MRSCont.nDatasets))
end
fileID = fopen(fullfile(MRSCont.outputFolder, 'LogFile.txt'),'a+');
for kk = 1:MRSCont.nDatasets
    msg = sprintf('Processing data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    fprintf(fileID,[reverseStr, msg]);
    if MRSCont.flags.isGUI        
        set(progressText,'String' ,sprintf('Processing data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets));
        drawnow
    end
    
    if ((MRSCont.flags.didProcess == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'processed') && (kk > length(MRSCont.processed.A))) || ~isfield(MRSCont.ver, 'Pro') || ~strcmp(MRSCont.ver.Pro,MRSCont.ver.CheckPro))

        %%% 1. GET RAW DATA %%%
        raw                         = MRSCont.raw{kk};                                          % Get the kk-th dataset

        %%% 1B. GET MM DATA %%% 
        if MRSCont.flags.hasMM
            raw_mm                         = MRSCont.raw_mm{kk};              % Get the kk-th dataset re_mm
            if raw_mm.averages > 1 && raw_mm.flags.averaged == 0 %re_mm
                [raw_mm,~,~]               = op_alignAverages(raw_mm, 1, 'n'); %re_mm
                raw_mm                     = op_averaging(raw_mm);            % Average re_mm
            else %re_mm
                raw_mm.flags.averaged  = 1; %re_mm
                raw_mm.dims.averages   = 0; %re_mm
            end
            [raw_mm,~]                     = op_ppmref(raw_mm,4.6,4.8,4.68);  % Reference to water @ 4.68 ppm  %re_mm            
        end  %re_mm
        
        

        %%% 2a. PHANTOM-SPECIFIC PRE-PROCESSING %%%
        % If this is phantom data (assuming room temperature), we want to
        % perform a few specific pre-processing steps.
        if MRSCont.flags.isPhantom
            % First, we undo phase cycling by dividing by the first data
            % point (this is mainly experimental at this point, but has
            % proved beneficial for phase-cycled GE data).
%             for rr = 1:raw.rawAverages
%                 phi = repelem(conj(raw.fids(1,rr))./abs(raw.fids(1,rr)),size(raw.fids,1));
%                 raw.fids(:,rr) = raw.fids(:,rr) .* phi';
%                 raw.specs = fftshift(fft(raw.fids,[],1));
%             end
            
            % Next, shift the entire metabolite spectrum by 0.15 ppm.
            % This doesn't have to be completely accurate, since additional
            % referencing steps are performed in the later stages of
            % post-processing and modelling, but we want the prominent singlets
            % to appear within 0.1 ppm of their expected in-vivo positions.
            phantomShiftPPM = 0.15 * raw.txfrq*1e-6;
            raw = op_freqshift(raw, -phantomShiftPPM);

            % Finally, apply some linebroadening. High-quality in-vitro
            % data may have linewidth lower than the simulated basis set
            % data.
            raw = op_filter(raw, 2);
            if MRSCont.flags.hasRef
                raw_ref = op_filter(raw_ref, 2);
            end
        end
        
        %%% 3. FREQUENCY/PHASE CORRECTION AND AVERAGING %%%

        if raw.averages > 1 && raw.flags.averaged == 0
            % Calculate starting values for spectral registration
             [refShift_ind_ini]=op_preref(raw,'unedited');
            if ~MRSCont.flags.isPhantom
                switch MRSCont.opts.SpecReg %Pick spectral registration method (default is Robust Spectral Registration)
                    case 'RobSpecReg'
                        [raw, fs, phs, weights, driftPre, driftPost]     = op_robustSpecReg(raw, 'unedited', 0,refShift_ind_ini); % Align and average
                    case 'RestrSpecReg'
                        [raw, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw, 'unedited', 0,refShift_ind_ini,0,MRSCont.opts.fit.range(1),MRSCont.opts.fit.range(2)); % Align and average
                    case 'none'
                        [raw, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw, 'unedited', 0,refShift_ind_ini,1); % Align and average   
                end                        
            else
                [raw, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw, 'unedited', 0,refShift_ind_ini,0,0.5,4.2); % Align and average
            end
            raw.specReg.fs              = fs; % save align parameters
            raw.specReg.phs             = phs; % save align parameters
            raw.specReg.weights         = weights{1}(1,:)'; % save align parameters);
            raw.specReg.weights         = raw.specReg.weights/max(raw.specReg.weights);
        else
            raw.flags.averaged  = 1;
            raw.dims.averages   = 0;
            raw.specReg.fs              = 0; % save align parameters
            raw.specReg.phs             = 0; % save align parameters
            raw.specReg.weights         = 1; % save align parameters
            driftPre = op_measureDrift(raw);
            driftPost = driftPre;
        end

        %%% 4. GET REFERENCE DATA / EDDY CURRENT CORRECTION %%%
        % If there are reference scans, load them here to allow eddy-current
        % correction of the raw data.
        if MRSCont.flags.hasRef
            raw_ref                         = MRSCont.raw_ref{kk};              % Get the kk-th dataset
            if raw_ref.averages > 1 && raw_ref.flags.averaged == 0
                [raw_ref,~,~]               = op_alignAverages(raw_ref, 1, 'n');
                raw_ref                     = op_averaging(raw_ref);            % Average
            else
                raw_ref.flags.averaged  = 1;
                raw_ref.dims.averages   = 0;
            end
                        
            if MRSCont.flags.hasMM
                [raw_mm,~]                   = op_eccKlose(raw_mm, raw_ref);        % Klose eddy current correction
            end
            [raw,raw_ref]                   = op_eccKlose(raw, raw_ref);        % Klose eddy current correction
            [raw_ref,~]                     = op_ppmref(raw_ref,4.6,4.8,4.68);  % Reference to water @ 4.68 ppm
            MRSCont.processed.ref{kk}       = raw_ref;                          % Save back to MRSCont container
        end
        
        %%% 5. DETERMINE POLARITY OF SPECTRUM (EG FOR MOIST WATER SUPP) %%%
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
            MRSCont.raw{kk} = op_ampScale(MRSCont.raw{kk},-1);
        end


        %%% 6. REMOVE RESIDUAL WATER %%%
        % Define different water removal frequency ranges, depending on
        % whether this is phantom data
        if MRSCont.flags.isPhantom
            waterRemovalFreqRange = [4.5 5];
            fracFID = 0.2;
        else
            waterRemovalFreqRange = [4.5 4.9];
            fracFID = 0.75;
        end
        % Apply iterative water filter
        raw = op_iterativeWaterFilter(raw, waterRemovalFreqRange, 32, fracFID*length(raw.fids), 0);
        
        if MRSCont.flags.hasMM %re_mm
            raw_mm = op_iterativeWaterFilter(raw_mm, waterRemovalFreqRange, 32, fracFID*length(raw_mm.fids), 0);
        end

        %%% 7. REFERENCE SPECTRUM CORRECTLY TO FREQUENCY AXIS AND PHASE SIEMENS
        %%% DATA
        [refShift, ~] = osp_CrChoReferencing(raw);
        [raw]             = op_freqshift(raw,-refShift);            % Reference spectra by cross-correlation     
        
        if MRSCont.flags.hasMM %re_mm
            [refShift_mm, ~] = fit_OspreyReferencingMM(raw_mm);
            [raw_mm]             = op_freqshift(raw_mm,-refShift_mm);            % Reference spectra by cross-correlation
            MRSCont.processed.mm{kk}       = raw_mm;                          % Save back to MRSCont container  %re_mm
        end

        
        % Save back to MRSCont container
        if strcmp(MRSCont.vendor,'Siemens')
            % Fit a double-Lorentzian to the Cr-Cho area, and phase the spectrum
            % with the negative phase of that fit
            [raw,globalPhase]       = op_phaseCrCho(raw, 1);
            raw.specReg.phs = raw.specReg.phs - globalPhase*180/pi;
        end
        MRSCont.processed.A{kk}     = raw;


        %%% 8. GET SHORT-TE WATER DATA %%%
        if MRSCont.flags.hasWater
            raw_w                           = MRSCont.raw_w{kk};                % Get the kk-th dataset
            if raw_w.averages > 1 && raw_w.flags.averaged == 0
                [raw_w,~,~]                 = op_alignAverages(raw_w, 1, 'n');
                raw_w                       = op_averaging(raw_w);              % Average
            else
                raw_w.flags.averaged    = 1;
                raw_w.dims.averages     = 0;
            end
            [raw_w,~]                       = op_eccKlose(raw_w, raw_w);        % Klose eddy current correction
            [raw_w,~]                       = op_ppmref(raw_w,4.6,4.8,4.68);    % Reference to water @ 4.68 ppm
            
            % Apply some linebroadening, if phantom data
            if MRSCont.flags.isPhantom
                raw_w = op_filter(raw_w, 2);    
            end
            
            MRSCont.processed.w{kk}         = raw_w; % Save back to MRSCont container
        end


        %%% 9. QUALITY CONTROL PARAMETERS %%%
        % Calculate some spectral quality metrics here;
        MRSCont.QM.SNR.A(kk)    = op_getSNR(MRSCont.processed.A{kk}); % NAA amplitude over noise floor             
        MRSCont.QM.FWHM.A(kk)   = op_getLW(MRSCont.processed.A{kk},1.8,2.2); % LW in Hz
        MRSCont.QM.drift.pre.A{kk}  = driftPre;
        MRSCont.QM.drift.post.A{kk} = driftPost;
        MRSCont.QM.freqShift.A(kk)  = refShift;
        MRSCont.QM.drift.pre.AvgDeltaCr.A(kk) = mean(driftPre - 3.02);
        MRSCont.QM.drift.post.AvgDeltaCr.A(kk) = mean(driftPost - 3.02);
        MRSCont.QM.res_water_amp.A(kk) = sum(MRSCont.processed.A{kk}.watersupp.amp);
        if MRSCont.flags.hasMM
            MRSCont.QM.SNR.mm(kk)    = op_getSNR(MRSCont.processed.mm{kk},0.7,1.1); % water amplitude over noise floor
            MRSCont.QM.FWHM.mm(kk)   = op_getLW(MRSCont.processed.mm{kk},0.7,1.1); % LW in Hz
        end
        if MRSCont.flags.hasRef
            MRSCont.QM.SNR.ref(kk)  = op_getSNR(MRSCont.processed.ref{kk},4.2,5.2); % water amplitude over noise floor
            MRSCont.QM.FWHM.ref(kk) = op_getLW(MRSCont.processed.ref{kk},4.2,5.2); % LW in Hz
        end
        if MRSCont.flags.hasWater
            MRSCont.QM.SNR.w(kk)    = op_getSNR(MRSCont.processed.w{kk},4.2,5.2); % water amplitude over noise floor           
            MRSCont.QM.FWHM.w(kk)   = op_getLW(MRSCont.processed.w{kk},4.2,5.2); % LW in Hz
        end
        
        
    end

    if MRSCont.flags.isGUI        
        waitbar(kk/MRSCont.nDatasets,progressbar,sprintf('Processed data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets))
    end
end
fprintf('... done.\n');
time = toc(refProcessTime);
if MRSCont.flags.isGUI        
    set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
    pause(1);
end
fprintf(fileID,'... done.\n Elapsed time %f seconds\n',time);
fclose(fileID);

%%% 10. SET FLAGS %%%
MRSCont.flags.avgsAligned       = 1;
MRSCont.flags.averaged          = 1;
MRSCont.flags.ECCed             = 1;
MRSCont.flags.waterRemoved      = 1;
MRSCont.runtime.Proc = time;
% Close any remaining open figures
close all;

end
