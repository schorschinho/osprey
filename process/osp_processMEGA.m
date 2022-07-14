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


warning('off','all');

% Parse input arguments
if nargin < 2
    target = MRSCont.opts.editTarget{1}; % GABA editing as default
end

%% Loop over all datasets
refProcessTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end
for kk = 1:MRSCont.nDatasets
    [~] = printLog('OspreyProcess',kk,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);     
    
    if ~(MRSCont.flags.didProcess == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'processed') && (kk > length(MRSCont.processed.A))) || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp) 
        %%% 1. GET RAW DATA %%%
        raw         = MRSCont.raw{kk};                                          % Get the kk-th dataset
        
        %%% 1B. GET MM DATA %%% 
        if MRSCont.flags.hasMM
            raw_mm                         = MRSCont.raw_mm{kk};              % Get the kk-th dataset 
            if raw_mm.averages > 1 && raw_mm.flags.averaged == 0 %
                [raw_mm,~,~]               = op_alignAverages(raw_mm, 1, 'n'); %
                raw_mm                     = op_averaging(raw_mm);            % Average 
            else %re_mm
                raw_mm.flags.averaged  = 1; %re_mm
                raw_mm.dims.averages   = 0; %re_mm
            end    
        end  %re_mm
        
        if raw.averages > 1 && raw.flags.averaged == 0
        % Calculate starting values for spectral registration
        [refShift_ind_ini]=op_preref(raw,'MEGA');
        % Perform robust spectral correction with weighted averaging.
        % This can obviously only be done, if the spectra have not been 
        % pre-averaged, i.e. in some older RDA and DICOM files (which should, 
        % generally, not be used).
            if ~MRSCont.flags.isPhantom
                if MRSCont.flags.isMRSI
                    if MRSCont.opts.MoCo.lb > 0
                        raw = op_filter(raw, MRSCont.opts.MoCo.lb);
                    end
                end
                switch MRSCont.opts.SpecReg %Pick spectral registration method (default is Robust Spectral Registration)
                    case 'RobSpecReg'
                        [raw, fs, phs, weights, driftPre, driftPost]     = op_robustSpecReg(raw, 'MEGA', 0,refShift_ind_ini); % Align and average
                    case 'RestrSpecReg'
                        if isfield(MRSCont.opts,'SpecRegRange')
                            [raw, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw, 'MEGA', 0,refShift_ind_ini,0,MRSCont.opts.SpecRegRange(1),MRSCont.opts.SpecRegRange(2)); % Align and average
                        else
                            [raw, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw, 'MEGA', 0,refShift_ind_ini,0,MRSCont.opts.fit.range(1),MRSCont.opts.fit.range(2)); % Align and average
                        end
                    case 'none'
                        [raw, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw, 'MEGA', 0,refShift_ind_ini,1); % Align and average   
                end
            else
                % Next, shift the entire metabolite spectrum by 0.15 ppm.
            % This doesn't have to be completely accurate, since additional
            % referencing steps are performed in the later stages of
            % post-processing and modelling, but we want the prominent singlets
            % to appear within 0.1 ppm of their expected in-vivo positions.
            phantomShiftPPM = 0.2 * raw.txfrq*1e-6;
            raw = op_freqshift(raw, -phantomShiftPPM);
            % Finally, apply some linebroadening. High-quality in-vitro
            % data may have linewidth lower than the simulated basis set
            % data.
            raw = op_filter(raw, 6);
                switch MRSCont.opts.SpecReg %Pick spectral registration method (default is Robust Spectral Registration)
                    case 'none'
                        [raw, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw, 'MEGA', 0,refShift_ind_ini,1); % Align and average  
                    otherwise
                        [raw, fs, phs, weights, driftPre, driftPost]   = op_SpecRegFreqRestrict(raw, 'MEGA', 0,refShift_ind_ini,0,0.5,4.2);
                end
            
            end
            raw.specReg.fs              = fs; % save align parameters
            raw.specReg.phs             = phs; % save align parameters
            raw.specReg.weights         = weights; % save align parameters            
        else
            if MRSCont.flags.isMRSI
                if MRSCont.opts.MoCo.lb > 0
                    raw = op_filter(raw, MRSCont.opts.MoCo.lb);
                end
            end
                
            raw.flags.averaged  = 1;
            raw.dims.averages   = 0;
            raw.specReg.fs              = zeros(1,2); % save align parameters
            raw.specReg.phs             = zeros(1,2); % save align parameters
            raw.specReg.weights{1,1}         = ones(512,1); % save align parameters
            raw.specReg.weights{1,2}         = ones(512,1); % save align parameters
            driftPre{1} = zeros(2,1);
            driftPre{2} = zeros(2,1);
            driftPost = driftPre;
            if MRSCont.flags.isPhantom
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
            end
        end

        % Get sub-spectra, depending on whether they are stored as such
        if raw.subspecs == 2
            raw_A   = op_takesubspec(raw,1);                    % Get first subspectrum
            raw_B   = op_takesubspec(raw,2);                    % Get second subspectrum
        else
            raw_A   = op_takeaverages(raw,1:2:raw.averages);    % Get first subspectrum
            raw_B   = op_takeaverages(raw,2:2:raw.averages);    % Get second subspectrum
        end
        
        if MRSCont.flags.hasMM     
            if raw.subspecs == 2
                raw_mm_A   = op_takesubspec(raw_mm,1);                    % Get first subspectrum
                raw_mm_B   = op_takesubspec(raw_mm,2);                    % Get second subspectrum
            else
                raw_mm_A   = op_takeaverages(raw_mm,1:2:raw_mm.averages);    % Get first subspectrum
                raw_mm_B   = op_takeaverages(raw_mm,2:2:raw_mm.averages);    % Get second subspectrum
            end
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
           [raw_A,~]               = op_eccKlose(raw_A, raw_ref);
           [raw_B,raw_ref]         = op_eccKlose(raw_B, raw_ref);        % Klose eddy current correction
           

            [raw_ref,~]                 = op_ppmref(raw_ref,4.6,4.8,4.68);  % Reference to water @ 4.68 ppm
            MRSCont.processed.ref{kk}   = raw_ref;                          % Save back to MRSCont container
        end
        
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
            if ~MRSCont.flags.isMRSI
                [raw_w,~]                   = op_eccKlose(raw_w, raw_w);        % Klose eddy current correction
            else
                [raw_w,~]=op_autophase(raw_w,2,2*4.68);
            end

            
            if MRSCont.flags.hasMM   
               [raw_mm_A,~]               = op_eccKlose(raw_mm_A, raw_w);
                [raw_mm_B,~]         = op_eccKlose(raw_mm_B, raw_w);        % Klose eddy current correction
            end  %re_mm

            % Save back to MRSCont container            
            [raw_w,~]                   = op_ppmref(raw_w,4.6,4.8,4.68);    % Reference to water @ 4.68 ppm
            MRSCont.processed.w{kk}     = raw_w;      
        end
        
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
            if MRSCont.flags.hasRef
                raw_ref = op_filter(raw_ref, 2);
            end
        end


        %%% 3. DETERMINE POLARITY OF SPECTRUM (EG FOR MOIST WATER SUPP) %%%
        % Automate determination whether the Cr peak has positive polarity.
        % For water suppression methods like MOIST, the residual water may
        % actually have negative polarity, but end up positive in the data, so
        % that the spectrum needs to be flipped.
        if isfield(MRSCont,'polResidCr')
            polResidCr = MRSCont.polResidCr;
        else
            raw_A_Cr    = op_freqrange(raw_A,2.8,3.2);
            % Determine the polarity of the respective peak: if the absolute of the
            % maximum minus the absolute of the minimum is positive, the polarity 
            % of the respective peak is positive; if the absolute of the maximum 
            % minus the absolute of the minimum is negative, the polarity is negative.
            polResidCr  = abs(max(real(raw_A_Cr.specs))) - abs(min(real(raw_A_Cr.specs)));            
        end
        if polResidCr < 0
            raw_A = op_ampScale(raw_A,-1);
        end
        if isfield(MRSCont,'polResidCr')
            polResidCr = MRSCont.polResidCr;
        else
            raw_B_Cr    = op_freqrange(raw_B,2.8,3.2);
            % Determine the polarity of the respective peak: if the absolute of the
            % maximum minus the absolute of the minimum is positive, the polarity 
            % of the respective peak is positive; if the absolute of the maximum 
            % minus the absolute of the minimum is negative, the polarity is negative.
            polResidCr  = abs(max(real(raw_B_Cr.specs))) - abs(min(real(raw_B_Cr.specs)));
            if isfield(MRSCont,'polResidCr')
                polResidCr = MRSCont.polResidCr;
            end
        end
        if polResidCr < 0
            raw_B = op_ampScale(raw_B,-1);
        end
        
        if MRSCont.flags.hasMM
            raw_mm_A_Cr    = op_freqrange(raw_mm_A,3.7,4.1);
            polResidCr  = abs(max(real(raw_mm_A_Cr.specs))) - abs(min(real(raw_mm_A_Cr.specs)));
            if polResidCr < 0
                raw_mm_A = op_ampScale(raw_mm_A,-1);
            end
            raw_mm_B_Cr    = op_freqrange(raw_mm_B,3.7,4.1);
            polResidCr  = abs(max(real(raw_mm_B_Cr.specs))) - abs(min(real(raw_mm_B_Cr.specs)));
            if polResidCr < 0
                raw_mm_B = op_ampScale(raw_mm_B,-1);
            end
        end

        %%% 4. DETERMINE ON/OFF STATUS
        % Classify the two sub-spectra such that the OFF spectrum is stored to
        % field A, and the ON spectrum is stored to field B.
        if isfield(MRSCont,'switchOrder')
            switchOrder = MRSCont.switchOrder;
            temp = raw_A;
            raw_A = raw_B;
            raw_B = temp;
        else
            [raw_A, raw_B, switchOrder]  = osp_onOffClassifyMEGA(raw_A, raw_B, target);
        end
        
        % Save drift information back to container
        if ~switchOrder
            MRSCont.QM.drift.pre.A{kk} = driftPre{1};
            MRSCont.QM.drift.pre.B{kk} = driftPre{2};
            MRSCont.QM.drift.post.A{kk} = driftPost{1};
            MRSCont.QM.drift.post.B{kk} = driftPost{1};
        else
            MRSCont.QM.drift.pre.A{kk} = driftPre{2};
            MRSCont.QM.drift.pre.B{kk} = driftPre{1};
            MRSCont.QM.drift.post.A{kk} = driftPost{2};
            MRSCont.QM.drift.post.B{kk} = driftPost{1};
        end
        raw_A.specReg.fs     = raw.specReg.fs(:,1); % save align parameters
        raw_B.specReg.fs     = raw.specReg.fs(:,2);
        raw_A.specReg.phs     = raw.specReg.phs(:,1);
        raw_B.specReg.phs     = raw.specReg.phs(:,2);
        raw_A.specReg.weights 	= raw.specReg.weights{1}(1,:);
        raw_B.specReg.weights    = raw.specReg.weights{2}(1,:);
        raw_A.specReg.weights = raw_A.specReg.weights'/(max(raw_A.specReg.weights));
        raw_B.specReg.weights = raw_B.specReg.weights'/(max(raw_B.specReg.weights));
        % Generate the frequency and phase plots for the entire experiment in
        % the correct order
        fs = [raw_A.specReg.fs, raw_B.specReg.fs]';
        fs = reshape(fs, [raw.rawAverages, 1]);
        phs = [raw_A.specReg.phs, raw_B.specReg.phs]';
        phs = reshape(phs, [raw.rawAverages, 1]);
        weights = [raw_A.specReg.weights, raw_B.specReg.weights]';
        weights = reshape(weights, [raw.rawAverages, 1]);
        MRSCont.raw{kk}.specReg.fs              = fs; % save align parameters
        MRSCont.raw{kk}.specReg.phs             = phs; % save align parameters
        MRSCont.raw{kk}.specReg.weights             = weights; % save align parameters 

        %%% 5. BUILD SUM AND DIFF SPECTRA %%%
        % Correct the frequency axis so that Cr appears at 3.027 ppm
        if ~MRSCont.flags.isMRSI
            [raw_A,~]       = op_phaseCrCho(raw_A, 1);
            [raw_B,~]       = op_phaseCrCho(raw_B, 1);
            temp_spec   = op_addScans(raw_A,raw_B);  
            [refShift_SubSpecAlign, ~] = osp_CrChoReferencing(temp_spec);
            % Apply initial referencing shift
            raw_A = op_freqshift(raw_A, -refShift_SubSpecAlign);
            raw_B = op_freqshift(raw_B, -refShift_SubSpecAlign);
            % Fit a double-Lorentzian to the Cr-Cho area, and phase the spectrum
            % with the negative phase of that fit
            [raw_A,ph]       = op_phaseCrCho(raw_A, 1);
            % Align the sub-spectra to one another by minimizing the difference
            % between the common 'reporter' signals.
        else
            MRSCont.opts.UnstableWater = 1;
            refShift_SubSpecAlign =0;
        end
        
        switch MRSCont.opts.SubSpecAlignment
            case 'L1Norm'
            [raw_A, raw_B]  = osp_editSubSpecAlignLNorm(raw_A, raw_B);
            case 'L2Norm'
            [raw_A, raw_B]  = osp_editSubSpecAlign(raw_A, raw_B, target,MRSCont.opts.UnstableWater);
            otherwise
                if ~MRSCont.flags.isMRSI
                    raw_B     = op_addphase(raw_B, -ph*180/pi, 0, raw_B.centerFreq, 1); 
                end
        end
        
        if MRSCont.flags.hasMM 
            diff1_mm           = op_addScans(raw_mm_B,raw_mm_A,1);
            [diff1_mm,~]          = op_autophase(diff1_mm,0.5,1.1);
            [refShift_mm, ~] = fit_OspreyReferencingMM(diff1_mm);
            [raw_mm_A]             = op_freqshift(raw_mm_A,-refShift_mm);            % Reference spectra by cross-correlation
            [raw_mm_B]             = op_freqshift(raw_mm_B,-refShift_mm);            % Reference spectra by cross-correlation
%             switch MRSCont.opts.SubSpecAlignment
%                 case 'L1Norm'
%                 [raw_mm_A, raw_mm_B]  = osp_editSubSpecAlignLNorm(raw_mm_A, raw_mm_B);
%                  [raw_mm_A, raw_mm_B] = osp_editSubSpecAlignGauss09(raw_mm_A, raw_mm_B);
%                 case 'L2Norm'
%                 [raw_mm_A, raw_mm_B]  = osp_editSubSpecAlign(raw_mm_A, raw_mm_B, target,MRSCont.opts.UnstableWater);
%             end
         end
        

        
        % Create the sum spectrum
        Sum             = op_addScans(raw_A,raw_B);
        if switchOrder
            Sum.flags.orderswitched = 1;
        else
            Sum.flags.orderswitched = 0;
        end
        Sum.specReg.fs = fs;
        Sum.specReg.phs = phs;
        Sum.specReg.weights = weights;
        % Create the GABA-edited difference spectrum
        diff1           = op_addScans(raw_B,raw_A,1);
        if switchOrder
            diff1.flags.orderswitched = 1;
        else
            diff1.flags.orderswitched = 0;
        end
        diff1.target = target;
        diff1.specReg.fs = fs;
        diff1.specReg.phs = phs;
        diff1.specReg.weights = weights;
        
        if MRSCont.flags.hasMM 
            Sum_mm             = op_addScans(raw_mm_A,raw_mm_B);
            diff1_mm           = op_addScans(raw_mm_B,raw_mm_A,1);
            if switchOrder
            diff1_mm.flags.orderswitched = 1;
            else
                diff1_mm.flags.orderswitched = 0;
            end
        end
        
        %%% 6. REMOVE RESIDUAL WATER %%%
        % Remove water and correct back to baseline.
        % The spectra sometimes become NaNs after filtering with too many
        % components. Loop over decreasing numbers of components here.        
        % Define different water removal frequency ranges, depending on
        % whether this is phantom data
        
        
%         ZF = raw_A.sz(1); % zerofill points
%         EP = raw_A.dwelltime;  % echo position 0~1
% 
%         %SW, Datapoints, Zerofill, Lipid_Left, Lipid_Right, damp, damp_range, spe_number, echopos
%         lip1 = generator_lipid(raw.spectralwidth(1),512, ZF, 375, 460, 10, 10, ZF, EP); % generation of lipid basic matrix
%         % SW, Datapoints, Zerofill, wat_frq_range, damp, damp_range, Watnum, echoposition
%         wat1 = generator_water(raw.spectralwidth(1), 512, ZF, 40, 10, 5, 20, EP); % generation of water basic matrix
%         lip2 = generator_lipid(raw_A.spectralwidth(1),512, ZF, 240, 460, 10, 30, ZF, EP); % generation of lipid basic matrix
%         
% 
%         specs = watersup_sim(flipud(raw_A.specs(:)), real(wat1), 3);
%         raw_A.specs(:) = flipud(specs);                       
%         raw_A.fids = ifft(fftshift(raw_A.specs,1),[],1);
%         
% %         specs = watersup_sim(flipud(raw_A.specs(:)), real(lip1), 3);
% %         raw_A.specs(:) = flipud(specs);                       
% %         raw_A.fids = ifft(fftshift(raw_A.specs,1),[],1);
%         
%         specs = watersup_sim(flipud(raw_B.specs(:)), real(wat1), 3);
%         raw_B.specs(:) = flipud(specs);                       
%         raw_B.fids = ifft(fftshift(raw_B.specs,1),[],1);
%         
% %         specs = watersup_sim(flipud(raw_B.specs(:)), real(lip1), 3);
% %         raw_B.specs(:) = flipud(specs);                       
% %         raw_B.fids = ifft(fftshift(raw_B.specs,1),[],1);
%         
%         specs = watersup_sim(flipud(Sum.specs(:)), real(wat1), 3);
%         Sum.specs(:) = flipud(specs);                       
%         Sum.fids = ifft(fftshift(Sum.specs,1),[],1);
%         
% %         specs = watersup_sim(flipud(Sum.specs(:)), real(lip1), 3);
% %         Sum.specs(:) = flipud(specs);                       
% %         Sum.fids = ifft(fftshift(Sum.specs,1),[],1);
%         
%         specs = watersup_sim(flipud(diff1.specs(:)), real(wat1), 3);
%         diff1.specs(:) = flipud(specs);                       
%         diff1.fids = ifft(fftshift(diff1.specs,1),[],1);
%         
% %         specs = watersup_sim(flipud(diff1.specs(:)), real(lip2), 3);
% %         diff1.specs(:) = flipud(specs);                       
% %         diff1.fids = ifft(fftshift(diff1.specs,1),[],1);
%         
%         raw_A     = op_fddccorr(raw_A,100); % Correct back to baseline
%         raw_B     = op_fddccorr(raw_B,100); % Correct back to baseline
%         Sum     = op_fddccorr(Sum,100); % Correct back to baseline
%         diff1     = op_fddccorr(diff1,100); % Correct back to baseline
%         
%         raw_A.watersupp.amp = 0;
%         raw_B.watersupp.amp = 0;
%         Sum.watersupp.amp = 0;
%         diff1.watersupp.amp = 0;
        
        
        if MRSCont.flags.isPhantom
            waterRemovalFreqRange = [4.5 5];
            fracFID = 0.5;
        else
            waterRemovalFreqRange = [4.2 6];
            fracFID = 0.3;
        end

        % Apply iterative water filter
        raw_A = op_iterativeWaterFilter(raw_A, waterRemovalFreqRange, 16, fracFID*length(raw.fids), 0);
        raw_B = op_iterativeWaterFilter(raw_B, waterRemovalFreqRange, 16, fracFID*length(raw.fids), 0);
        diff1 = op_iterativeWaterFilter(diff1, waterRemovalFreqRange, 16, fracFID*length(raw.fids), 0);
        Sum = op_iterativeWaterFilter(Sum, waterRemovalFreqRange, 16, fracFID*length(raw.fids), 0);
        if MRSCont.flags.hasMM
            raw_mm_A = op_iterativeWaterFilter(raw_mm_A, waterRemovalFreqRange, 32, fracFID*length(raw.fids), 0);
            raw_mm_B = op_iterativeWaterFilter(raw_mm_B, waterRemovalFreqRange, 32, fracFID*length(raw.fids), 0);
            diff1_mm = op_iterativeWaterFilter(diff1_mm, waterRemovalFreqRange, 32, fracFID*length(raw.fids), 0);
            Sum_mm = op_iterativeWaterFilter(Sum_mm, waterRemovalFreqRange, 32, fracFID*length(raw.fids), 0);
        end
       

        %%% 7. REFERENCE SPECTRUM CORRECTLY TO FREQUENCY AXIS 
        % Reference resulting data correctly and consistently
        if ~MRSCont.flags.isMRSI
            [refShift_final, ~] = osp_CrChoReferencing(Sum);
            [raw_A]             = op_freqshift(raw_A,-refShift_final);            % Apply same shift to edit-OFF
            [raw_B]             = op_freqshift(raw_B,-refShift_final);            % Apply same shift to edit-OFF
            [diff1]             = op_freqshift(diff1,-refShift_final);            % Apply same shift to diff1
            [Sum]               = op_freqshift(Sum,-refShift_final);              % Apply same shift to sum
        else
            refShift_final = 0;
        end
        if MRSCont.flags.hasMM
            [diff1_mm,~]          = op_autophase(diff1_mm,0.5,1.1);
            [refShift_mm, ~] = fit_OspreyReferencingMM(diff1_mm);
            [raw_mm_A]             = op_freqshift(raw_mm_A,-refShift_mm);            % Reference spectra by cross-correlation
            [raw_mm_B]             = op_freqshift(raw_mm_B,-refShift_mm);            % Reference spectra by cross-correlation
            [diff1_mm]             = op_freqshift(diff1_mm,-refShift_mm);            % Reference spectra by cross-correlation
            [Sum_mm]             = op_freqshift(Sum_mm,-refShift_mm);            % Reference spectra by cross-correlation
            MRSCont.processed.A_mm{kk}       = raw_mm_A;                          % Save back to MRSCont container  %re_mm
            MRSCont.processed.B_mm{kk}       = raw_mm_B;                          % Save back to MRSCont container  %re_mm           
            MRSCont.processed.diff1_mm{kk}       = diff1_mm;                          % Save back to MRSCont container  %re_mm
            MRSCont.processed.sum_mm{kk}       = Sum_mm;                          % Save back to MRSCont container  %re_mm
        end


        %%% 8. SAVE BACK TO MRSCONT CONTAINER
        MRSCont.processed.A{kk}     = raw_A;                                    % Save edit-OFF back to MRSCont container
        MRSCont.processed.B{kk}     = raw_B;                                    % Save edit-ON back to MRSCont container
        MRSCont.processed.diff1{kk} = diff1;                                    % Save diff1 back to MRSCont container
        MRSCont.processed.sum{kk}   = Sum;                                      % Save sum back to MRSCont container
        


%         %%% 9. GET SHORT-TE WATER DATA %%%
%         if MRSCont.flags.hasWater
%             % Some formats end up having subspectra in their reference scans
%             % (e.g. Philips), as well as empty lines. Intercept these cases
%             % here.
%             raw_w                       = MRSCont.raw_w{kk};                % Get the kk-th dataset
%             if raw_w.subspecs > 1
%                 raw_w_A                 = op_takesubspec(raw_w,1);
%                 [raw_w_A]               = op_rmempty(raw_w_A);              % Remove empty lines
%                 raw_w_B                 = op_takesubspec(raw_w,2);
%                 [raw_w_A]               = op_rmempty(raw_w_A);              % Remove empty lines
%                 raw_w                   = op_concatAverages(raw_w_A,raw_w_B);
%             end
%             if ~raw_w.flags.averaged
%                 [raw_w,~,~]             = op_alignAverages(raw_w,1,'n');    % Align averages
%                 raw_w                   = op_averaging(raw_w);              % Average
%             end
%             if ~MRSCont.flags.isMRSI
%                 [raw_w,~]                   = op_eccKlose(raw_w, raw_w);        % Klose eddy current correction
%             else
%                 [raw_w,~]=op_autophase(raw_w,2,2*4.68);
%             end
%             [raw_w,~]                   = op_ppmref(raw_w,4.6,4.8,4.68);    % Reference to water @ 4.68 ppm
%             MRSCont.processed.w{kk}     = raw_w;                            % Save back to MRSCont container
%         end


        %%% 10. QUALITY CONTROL PARAMETERS %%%
        SubSpec = {'A','B','diff1','sum'};       
        SNRRange = {[1.8,2.2],[2.8,3.2],[2.8,3.2],[2.8,3.2]};
        if MRSCont.flags.hasRef
            SubSpec{end+1} = 'ref';
            SNRRange{end+1} = [4.2,5.2];
        end
        if MRSCont.flags.hasWater
            SubSpec{end+1} = 'w';
            SNRRange{end+1} = [4.2,5.2];
        end
        % Calculate some spectral quality metrics here;
        for ss = 1 : length(SubSpec)          
            MRSCont.QM.SNR.(SubSpec{ss})(kk)    = op_getSNR(MRSCont.processed.(SubSpec{ss}){kk},SNRRange{ss}(1),SNRRange{ss}(2));       
            MRSCont.QM.FWHM.(SubSpec{ss})(kk)   = op_getLW(MRSCont.processed.(SubSpec{ss}){kk},SNRRange{ss}(1),SNRRange{ss}(2)); % in Hz       
            if ~(strcmp(SubSpec{ss},'ref') || strcmp(SubSpec{ss},'w'))
                MRSCont.QM.freqShift.(SubSpec{ss})(kk)  = refShift_SubSpecAlign + refShift_final;       
                MRSCont.QM.res_water_amp.(SubSpec{ss})(kk) = sum(MRSCont.processed.(SubSpec{ss}){kk}.watersupp.amp);  
                if strcmp(SubSpec{ss},'diff1') || strcmp(SubSpec{ss},'sum')
                    MRSCont.QM.drift.pre.(SubSpec{ss}){kk}  = reshape([MRSCont.QM.drift.pre.A{kk}'; MRSCont.QM.drift.pre.B{kk}'], [], 1)';
                    MRSCont.QM.drift.post.(SubSpec{ss}){kk} = reshape([MRSCont.QM.drift.post.A{kk}'; MRSCont.QM.drift.post.B{kk}'], [], 1)';
                end
                MRSCont.QM.drift.pre.AvgDeltaCr.(SubSpec{ss})(kk) = mean(MRSCont.QM.drift.pre.(SubSpec{ss}){kk} - 3.02);
                MRSCont.QM.drift.post.AvgDeltaCr.(SubSpec{ss})(kk) = mean(MRSCont.QM.drift.pre.(SubSpec{ss}){kk} - 3.02);
            end
        end
      
    end
end

time = toc(refProcessTime);
[~] = printLog('done',time,MRSCont.nDatasets,progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI); 


%%% 11. SET FLAGS %%%
MRSCont.flags.avgsAligned   = 1;
MRSCont.flags.averaged      = 1;
MRSCont.flags.ECCed         = 1;
MRSCont.flags.waterRemoved  = 1;
MRSCont.runtime.Proc = time;
% Close any remaining open figures
close all;

end