function [MRSCont] = osp_processHERCULES(MRSCont,target1,target2)
%% [MRSCont] = osp_processHERCULES(MRSCont)
%   This function performs the following steps to process HERCULES-edited
%   (4-step) MRS data:
%       - Alignment of individual averages using robust spectral registration
%       - Averaging
%       - Removal of residual water using HSVD filtering
%       - Klose Eddy current correction (if a reference scan is provided)
%       - Automated zero-order phase correction
%       - Correct referencing of the ppm frequency axis
%
%   USAGE:
%       [MRSCont] = osp_processHERMES(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-08-15)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-08-20: First version of the code.


warning('off','all');

% Parse input arguments
if nargin < 2
    if (strcmp(MRSCont.opts.editTarget{1},'HERCULES1')||strcmp(MRSCont.opts.editTarget{1},'HERCULES2'))
        MRSCont.opts.editTarget = {'GABA','GSH'};
    end
    target1 = MRSCont.opts.editTarget{1}; 
    target2 = MRSCont.opts.editTarget{2};     
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
        raw         = MRSCont.raw{kk};                          % Get the kk-th dataset

        % Get sub-spectra, depending on whether they are stored as such
        if raw.subspecs == 4

            raw_A   = op_takesubspec(raw,1);                    % Get first subspectrum
            raw_B   = op_takesubspec(raw,2);                    % Get second subspectrum
            raw_C   = op_takesubspec(raw,3);                    % Get third subspectrum
            raw_D   = op_takesubspec(raw,4);                    % Get fourth subspectrum

        else

            raw_A   = op_takeaverages(raw,1:4:raw.averages);    % Get first subspectrum
            raw_B   = op_takeaverages(raw,2:4:raw.averages);    % Get second subspectrum
            raw_C   = op_takeaverages(raw,3:4:raw.averages);    % Get third subspectrum
            raw_D   = op_takeaverages(raw,4:4:raw.averages);    % Get fourth subspectrum
        end

        if raw.averages > 1 && raw.flags.averaged == 0
        % Calculate starting values for spectral registration
        [refShift_ind_ini]=op_preref(raw,'HERCULES');
        % Perform robust spectral correction with weighted averaging.
        % This can obviously only be done, if the spectra have not been 
        % pre-averaged, i.e. in some older RDA and DICOM files (which should, 
        % generally, not be used). For phantom scans the registration
        % window is restricted to a certain frequency window.
        if ~MRSCont.flags.isPhantom
            switch MRSCont.opts.SpecReg %Pick spectral registration method (default is Robust Spectral Registration)
                case 'RobSpecReg'
                    [raw_A, fs_A, phs_A, weights_A, driftPreA, driftPostA]   = op_robustSpecReg(raw_A, 'HERCULES', 0,refShift_ind_ini);
                    [raw_B, fs_B, phs_B, weights_B, driftPreB, driftPostB]   = op_robustSpecReg(raw_B, 'HERCULES', 0, refShift_ind_ini);                  
                    [raw_C, fs_C, phs_C, weights_C, driftPreC, driftPostC]   = op_robustSpecReg(raw_C, 'HERCULES', 0, refShift_ind_ini);                    
                    [raw_D, fs_D, phs_D, weights_D, driftPreD, driftPostD]   = op_robustSpecReg(raw_D, 'HERCULES', 0, refShift_ind_ini);
                case 'RestrSpecReg'
                    [~, ~, ~, ~, commuteOrder] = osp_onOffClassifyHERMES(temp_rawA, temp_rawB, temp_rawC, temp_rawD);
                    temp.raw_A = raw_A;
                    temp.raw_B = raw_B;
                    temp.raw_C = raw_C;
                    temp.raw_D = raw_D;
                    inputVars = {'raw_A', 'raw_B', 'raw_C', 'raw_D'};
                    eval(['raw_A = temp.' inputVars{commuteOrder(1)} ';']);
                    eval(['raw_B = temp.' inputVars{commuteOrder(2)} ';']);
                    eval(['raw_C = temp.' inputVars{commuteOrder(3)} ';']);
                    eval(['raw_D = temp.' inputVars{commuteOrder(4)} ';']);
                    [raw_A, fs_A, phs_A, weights_A, driftPreA, driftPostA]   = op_SpecRegFreqRestrict(raw_A, 'HERCULES', 0,refShift_ind_ini,0,0.5,3.9);
                    [raw_B, fs_B, phs_B, weights_B, driftPreB, driftPostB]   = op_SpecRegFreqRestrict(raw_B, 'HERCULES', 0, refShift_ind_ini,0,0.5,3.8);                  
                    [raw_C, fs_C, phs_C, weights_C, driftPreC, driftPostC]   = op_SpecRegFreqRestrict(raw_C, 'HERCULES', 0, refShift_ind_ini,0,1.85,3.9);                    
                    [raw_D, fs_D, phs_D, weights_D, driftPreD, driftPostD]   = op_SpecRegFreqRestrict(raw_D, 'HERCULES', 0, refShift_ind_ini,0,1.85,3.8);
                case 'none'
                    [raw_A, fs_A, phs_A, weights_A, driftPreA, driftPostA]   = op_SpecRegFreqRestrict(raw_A, 'HERCULES', 0,refShift_ind_ini,1);
                    [raw_B, fs_B, phs_B, weights_B, driftPreB, driftPostB]   = op_SpecRegFreqRestrict(raw_B, 'HERCULES', 0, refShift_ind_ini,1);                  
                    [raw_C, fs_C, phs_C, weights_C, driftPreC, driftPostC]   = op_SpecRegFreqRestrict(raw_C, 'HERCULES', 0, refShift_ind_ini,1);                    
                    [raw_D, fs_D, phs_D, weights_D, driftPreD, driftPostD]   = op_SpecRegFreqRestrict(raw_D, 'HERCULES', 0, refShift_ind_ini,1);   
            end
        else
            [~, ~, ~, ~, commuteOrder] = osp_onOffClassifyHERMES(temp_rawA, temp_rawB, temp_rawC, temp_rawD);
            temp.raw_A = raw_A;
            temp.raw_B = raw_B;
            temp.raw_C = raw_C;
            temp.raw_D = raw_D;
            inputVars = {'raw_A', 'raw_B', 'raw_C', 'raw_D'};
            eval(['raw_A = temp.' inputVars{commuteOrder(1)} ';']);
            eval(['raw_B = temp.' inputVars{commuteOrder(2)} ';']);
            eval(['raw_C = temp.' inputVars{commuteOrder(3)} ';']);
            eval(['raw_D = temp.' inputVars{commuteOrder(4)} ';']);
            % Next, shift the entire metabolite spectrum by 0.15 ppm.
            % This doesn't have to be completely accurate, since additional
            % referencing steps are performed in the later stages of
            % post-processing and modelling, but we want the prominent singlets
            % to appear within 0.1 ppm of their expected in-vivo positions.
            phantomShiftPPM = 0.15 * raw_A.txfrq*1e-6;
            raw_A = op_freqshift(raw_A, -phantomShiftPPM);
            raw_B = op_freqshift(raw_B, -phantomShiftPPM);
            raw_C = op_freqshift(raw_C, -phantomShiftPPM);
            raw_D = op_freqshift(raw_D, -phantomShiftPPM);
            % Finally, apply some linebroadening. High-quality in-vitro
            % data may have linewidth lower than the simulated basis set
            % data.
            raw_A = op_filter(raw_A, 2);
            raw_B = op_filter(raw_B, 2);
            raw_C = op_filter(raw_C, 2);
            raw_D = op_filter(raw_D, 2);
            [raw_A, fs_A, phs_A, weights_A, driftPreA, driftPostA]   = op_SpecRegFreqRestrict(raw_A, 'HERCULES', 0,refShift_ind_ini,0,0.5,3.9);
            [raw_B, fs_B, phs_B, weights_B, driftPreB, driftPostB]   = op_SpecRegFreqRestrict(raw_B, 'HERCULES', 0, refShift_ind_ini,0,0.5,3.8);                  
            [raw_C, fs_C, phs_C, weights_C, driftPreC, driftPostC]   = op_SpecRegFreqRestrict(raw_C, 'HERCULES', 0, refShift_ind_ini,0,1.85,3.9);                    
            [raw_D, fs_D, phs_D, weights_D, driftPreD, driftPostD]   = op_SpecRegFreqRestrict(raw_D, 'HERCULES', 0, refShift_ind_ini,0,1.85,3.8);
        end
            
            
        else
            raw.flags.averaged  = 1;
            raw.dims.averages   = 0;
            raw.specReg.fs              = zeros(1,2); % save align parameters
            raw.specReg.phs             = zeros(1,2); % save align parameters
            raw.specReg.weights{1}         = ones(1,1); % save align parameters
            raw.specReg.weights{2}         = ones(1,1); % save align parameters
            driftPre{1} = 0;
            driftPre{2} = 0;
            driftPost = driftPre;
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
                raw_ref_C               = op_takesubspec(raw_ref,3);
                [raw_ref_C]             = op_rmempty(raw_ref_C);            % Remove empty lines
                raw_ref_D               = op_takesubspec(raw_ref,4);
                [raw_ref_D]             = op_rmempty(raw_ref_D);            % Remove empty lines
                raw_ref                 = op_concatAverages(raw_ref_A,raw_ref_B,raw_ref_C,raw_ref_D);            
            end
            if ~raw_ref.flags.averaged
                [raw_ref]             = op_rmempty(raw_ref); 
                [raw_ref,~,~]           = op_alignAverages(raw_ref,1,'n');  % Align averages
                raw_ref                 = op_averaging(raw_ref);            % Average
            end

            % Apply Klose eddy current correction
            [raw_A,~]                   = op_eccKlose(raw_A, raw_ref);
            [raw_B,~]                   = op_eccKlose(raw_B, raw_ref);
            [raw_C,~]                   = op_eccKlose(raw_C, raw_ref);
            [raw_D,raw_ref]             = op_eccKlose(raw_D, raw_ref);

            [raw_ref,~]                 = op_ppmref(raw_ref,4.6,4.8,4.68);  % Reference to water @ 4.68 ppm
            MRSCont.processed.ref{kk}   = raw_ref;                          % Save back to MRSCont container
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
            % Finally, apply some linebroadening. High-quality in-vitro
            % data may have linewidth lower than the simulated basis set
            % data.
            if MRSCont.flags.hasRef
                raw_ref = op_filter(raw_ref, 2);
            end
        end


        %%% 3. DETERMINE POLARITY OF SPECTRUM (EG FOR MOIST WATER SUPP) %%%
        % Automate determination whether the Cr peak has positive polarity.
        % For water suppression methods like MOIST, the residual water may
        % actually have negative polarity, but end up positive in the data, so
        % that the spectrum needs to be flipped.
        raw_A_Cr    = op_freqrange(raw_A,2.8,3.2);
        polResidCr  = abs(max(real(raw_A_Cr.specs))) - abs(min(real(raw_A_Cr.specs)));
        if polResidCr < 0
            raw_A = op_ampScale(raw_A,-1);
        end

        raw_B_Cr    = op_freqrange(raw_B,2.8,3.2);
        polResidCr  = abs(max(real(raw_B_Cr.specs))) - abs(min(real(raw_B_Cr.specs)));
        if polResidCr < 0
            raw_B = op_ampScale(raw_B,-1);
        end

        raw_C_Cr    = op_freqrange(raw_C,2.8,3.2);
        polResidCr  = abs(max(real(raw_C_Cr.specs))) - abs(min(real(raw_C_Cr.specs)));
        if polResidCr < 0
            raw_C = op_ampScale(raw_C,-1);
        end

        raw_D_Cr    = op_freqrange(raw_D,2.8,3.2);
        polResidCr  = abs(max(real(raw_D_Cr.specs))) - abs(min(real(raw_D_Cr.specs)));
        if polResidCr < 0
            raw_D = op_ampScale(raw_D,-1);
        end
        %%% 4. DETERMINE ON/OFF STATUS
        % Classify the four sub-spectra such that the OFF spectrum is stored to
        % field A, and the ON spectrum is stored to field B.
        [raw_A, raw_B, raw_C, raw_D, commuteOrder] = osp_onOffClassifyHERMES(raw_A, raw_B, raw_C, raw_D);
        subSpecNames = {'A', 'B', 'C', 'D'};
        % Save drift information back to container
        eval(['MRSCont.QM.drift.pre.A{kk} = driftPre' subSpecNames{commuteOrder(1)} ';']);
        eval(['MRSCont.QM.drift.pre.B{kk} = driftPre' subSpecNames{commuteOrder(2)} ';']);
        eval(['MRSCont.QM.drift.pre.C{kk} = driftPre' subSpecNames{commuteOrder(3)} ';']);
        eval(['MRSCont.QM.drift.pre.D{kk} = driftPre' subSpecNames{commuteOrder(4)} ';']);
        eval(['MRSCont.QM.drift.post.A{kk} = driftPost' subSpecNames{commuteOrder(1)} ';']);
        eval(['MRSCont.QM.drift.post.B{kk} = driftPost' subSpecNames{commuteOrder(2)} ';']);
        eval(['MRSCont.QM.drift.post.C{kk} = driftPost' subSpecNames{commuteOrder(3)} ';']);
        eval(['MRSCont.QM.drift.post.D{kk} = driftPost' subSpecNames{commuteOrder(4)} ';']);
        % Generate the drift plot for the entire experiment in
        % the correct order
        driftPre = [MRSCont.QM.drift.pre.A{kk}, MRSCont.QM.drift.pre.B{kk}, MRSCont.QM.drift.pre.C{kk}, MRSCont.QM.drift.pre.D{kk}]';
        try
            driftPre = reshape(driftPre, [raw.averages, 1]);
        catch
            driftPre = reshape(driftPre, [raw.rawAverages, 1]);
        end
        MRSCont.QM.drift.pre.diff1{kk}  = driftPre;
        MRSCont.QM.drift.pre.diff2{kk}  = driftPre;
        MRSCont.QM.drift.pre.sum{kk}    = driftPre;
        driftPost = [MRSCont.QM.drift.post.A{kk}, MRSCont.QM.drift.post.B{kk}, MRSCont.QM.drift.post.C{kk}, MRSCont.QM.drift.post.D{kk}]';
        try
            driftPost = reshape(driftPost, [raw.averages, 1]);
        catch
            driftPost = reshape(driftPost, [raw.rawAverages, 1]);
        end
        MRSCont.QM.drift.post.diff1{kk}  = driftPost;
        MRSCont.QM.drift.post.diff2{kk}  = driftPost;
        MRSCont.QM.drift.post.sum{kk}    = driftPost;

        eval(['raw_A.specReg.fs     = fs_' subSpecNames{commuteOrder(1)} ';']); % save align parameters
        eval(['raw_B.specReg.fs     = fs_' subSpecNames{commuteOrder(2)} ';']);
        eval(['raw_C.specReg.fs   	= fs_' subSpecNames{commuteOrder(3)} ';']);
        eval(['raw_D.specReg.fs 	= fs_' subSpecNames{commuteOrder(4)} ';']);
        eval(['raw_A.specReg.phs 	= phs_' subSpecNames{commuteOrder(1)} ';']);
        eval(['raw_B.specReg.phs    = phs_' subSpecNames{commuteOrder(2)} ';']);
        eval(['raw_C.specReg.phs    = phs_' subSpecNames{commuteOrder(3)} ';']);
        eval(['raw_D.specReg.phs    = phs_' subSpecNames{commuteOrder(4)} ';']);
        eval(['raw_A.specReg.weights 	= weights_' subSpecNames{commuteOrder(1)} '{1}(1,:);']);
        eval(['raw_B.specReg.weights    = weights_' subSpecNames{commuteOrder(2)} '{1}(1,:);']);
        eval(['raw_C.specReg.weights    = weights_' subSpecNames{commuteOrder(3)} '{1}(1,:);']);
        eval(['raw_D.specReg.weights    = weights_' subSpecNames{commuteOrder(4)} '{1}(1,:);']);
        raw_A.specReg.weights = raw_A.specReg.weights'/(max(raw_A.specReg.weights));
        raw_B.specReg.weights = raw_B.specReg.weights'/(max(raw_B.specReg.weights));
        raw_C.specReg.weights = raw_C.specReg.weights'/(max(raw_C.specReg.weights));
        raw_D.specReg.weights = raw_D.specReg.weights'/(max(raw_D.specReg.weights));
        % Generate the frequency and phase plots for the entire experiment in
        % the correct order
        fs = [raw_A.specReg.fs, raw_B.specReg.fs, raw_C.specReg.fs, raw_D.specReg.fs]';
        fs = reshape(fs, [raw.rawAverages, 1]);
        phs = [raw_A.specReg.phs, raw_B.specReg.phs, raw_C.specReg.phs, raw_D.specReg.phs]';
        phs = reshape(phs, [raw.rawAverages, 1]);
        weights = [raw_A.specReg.weights, raw_B.specReg.weights, raw_C.specReg.weights, raw_D.specReg.weights]';
        weights = reshape(weights, [raw.rawAverages, 1]);
        MRSCont.raw{kk}.specReg.fs              = fs; % save align parameters
        MRSCont.raw{kk}.specReg.phs             = phs; % save align parameters
        MRSCont.raw{kk}.specReg.weights             = weights; % save align parameters

        %%% 5. BUILD SUM AND DIFF SPECTRA %%%
        % Correct the frequency axis so that Cr appears at 3.027 ppm
        temp_spec   = op_addScans(raw_A,raw_B);
        temp_spec   = op_addScans(temp_spec,raw_C);
        temp_spec   = op_addScans(temp_spec,raw_D);    
        [refShift_SubSpecAlign, ~] = osp_CrChoReferencing(temp_spec);
        % Apply initial referencing shift
        raw_A = op_freqshift(raw_A, -refShift_SubSpecAlign);
        raw_B = op_freqshift(raw_B, -refShift_SubSpecAlign);
        raw_C = op_freqshift(raw_C, -refShift_SubSpecAlign);
        raw_D = op_freqshift(raw_D, -refShift_SubSpecAlign);
        % Fit a double-Lorentzian to the Cr-Cho area, and phase the spectrum
        % with the negative phase of that fit
        [raw_A,~] = op_phaseCrCho(raw_A, 1); 
        % Align the sub-spectra to one another by minimizing the difference
        % between the common 'reporter' signals.
        [raw_A, raw_B, raw_C, raw_D] = osp_editSubSpecAlign(raw_A, raw_B, raw_C, raw_D,MRSCont.opts.UnstableWater);
        % Create the sum spectrum
        Sum     = op_addScans(raw_A,raw_B);
        Sum     = op_addScans(Sum,raw_C);
        Sum     = op_addScans(Sum,raw_D);
        Sum.commuteOrder = commuteOrder;
        Sum.specReg.fs = fs;
        Sum.specReg.phs = phs;
        Sum.specReg.weights = weights;
        % Create the GABA-edited difference spectrum
        diff1   = op_addScans(raw_B,raw_D);
        diff1   = op_addScans(diff1,raw_A,1);
        diff1   = op_addScans(diff1,raw_C,1);
        diff1.target = target1;
        diff1.commuteOrder = commuteOrder;
        diff1.specReg.fs = fs;
        diff1.specReg.phs = phs;
        diff1.specReg.weights = weights;
        % Create the GSH-edited difference spectrum
        diff2   = op_addScans(raw_C,raw_D);
        diff2   = op_addScans(diff2,raw_A,1);
        diff2   = op_addScans(diff2,raw_B,1);
        diff2.target = target2;
        diff2.commuteOrder = commuteOrder;
        diff2.specReg.fs = fs;
        diff2.specReg.phs = phs;
        diff2.specReg.weights = weights;

        %%% 6. REMOVE RESIDUAL WATER %%%
        % Remove water and correct back to baseline.
        % The spectra sometimes become NaNs after filtering with too many
        % components. Loop over decreasing numbers of components here.        
        % Define different water removal frequency ranges, depending on
        % whether this is phantom data
        if MRSCont.flags.isPhantom
            waterRemovalFreqRange = [4.5 5];
            fracFID = 0.2;
            Kinit = 32;
        else
            if ~MRSCont.flags.isPRIAM
                waterRemovalFreqRange = [4.5 4.9];
                fracFID = 0.75;
                Kinit = 32;
            else
%                 waterRemovalFreqRange = [4.3 5.5];
%                 fracFID = 0.75;
%                 Kinit = 45;
                waterRemovalFreqRange = [4.5 4.9];
                fracFID = 0.75;
                Kinit = 32;
            end
        end
        % Apply iterative water filter
        raw_A = op_iterativeWaterFilter(raw_A, waterRemovalFreqRange, Kinit, fracFID*length(raw.fids), 0);
        raw_B = op_iterativeWaterFilter(raw_B, waterRemovalFreqRange, Kinit, fracFID*length(raw.fids), 0);
        raw_C = op_iterativeWaterFilter(raw_C, waterRemovalFreqRange, Kinit, fracFID*length(raw.fids), 0);
        raw_D = op_iterativeWaterFilter(raw_D, waterRemovalFreqRange, Kinit, fracFID*length(raw.fids), 0);
        diff1 = op_iterativeWaterFilter(diff1, waterRemovalFreqRange, Kinit, fracFID*length(raw.fids), 0);
        diff2 = op_iterativeWaterFilter(diff2, waterRemovalFreqRange, Kinit, fracFID*length(raw.fids), 0);
        Sum = op_iterativeWaterFilter(Sum, waterRemovalFreqRange, Kinit, fracFID*length(raw.fids), 0);


        %%% 7. REFERENCE SPECTRUM CORRECTLY TO FREQUENCY AXIS 
        % Reference resulting data correctly and consistently
        %[raw_A, refShift]   = op_ppmref(raw_A,1.8,2.2,2.008);          % Reference OFF-OFF spectrum to NAA @ 2.008 ppm
        temp_spec = Sum;
        [refShift_final, ~] = osp_CrChoReferencing(temp_spec);% Reference OFF-OFF spectrum to NAA @ 2.008 ppm
        [raw_A]             = op_freqshift(raw_A,-refShift_final);            % Apply same shift to ON-OFF
        [raw_B]             = op_freqshift(raw_B,-refShift_final);            % Apply same shift to ON-OFF
        [raw_C]             = op_freqshift(raw_C,-refShift_final);            % Apply same shift to OFF-ON
        [raw_D]             = op_freqshift(raw_D,-refShift_final);            % Apply same shift to ON-ON
        [diff1]             = op_freqshift(diff1,-refShift_final);            % Apply same shift to diff1
        [diff2]             = op_freqshift(diff2,-refShift_final);            % Apply same shift to diff2
        [Sum]               = op_freqshift(Sum,-refShift_final);              % Apply same shift to sum

        % Add commuteOrder
        raw_A.commuteOrder = commuteOrder(1);
        raw_B.commuteOrder = commuteOrder(2);
        raw_C.commuteOrder = commuteOrder(3);
        raw_D.commuteOrder = commuteOrder(4);

        %%% 8. SAVE BACK TO MRSCONT CONTAINER
        MRSCont.processed.A{kk}     = raw_A;                                    % Save OFF-OFF back to MRSCont container
        MRSCont.processed.B{kk}     = raw_B;                                    % Save ON-OFF back to MRSCont container
        MRSCont.processed.C{kk}     = raw_C;                                    % Save OFF-ON back to MRSCont container
        MRSCont.processed.D{kk}     = raw_D;                                    % Save ON-ON back to MRSCont container
        MRSCont.processed.diff1{kk} = diff1;                                    % Save diff1 back to MRSCont container
        MRSCont.processed.diff2{kk} = diff2;                                    % Save diff2 back to MRSCont container
        MRSCont.processed.sum{kk}   = Sum;                                      % Save sum back to MRSCont container


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
                [raw_w,fs_A,phs_A]      = op_alignAverages(raw_w,1,'n');    % Align averages
                raw_w                   = op_averaging(raw_w);              % Average
            end
            [raw_w,~]                   = op_eccKlose(raw_w, raw_w);        % Klose eddy current correction
            [raw_w,~]                   = op_ppmref(raw_w,4.6,4.8,4.68);    % Reference to water @ 4.68 ppm
            MRSCont.processed.w{kk}     = raw_w;                            % Save back to MRSCont container
        end


        %%% 10. QUALITY CONTROL PARAMETERS %%%
        % Calculate some spectral quality metrics here;
        SubSpec = {'A','B','C','D','diff1','diff2','sum'};       
        SNRRange = {[1.8,2.2],[2.8,3.2],[2.8,3.2],[1.8,2.2],[2.8,3.2],[2.8,3.2],[2.8,3.2]};
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
            MRSCont.QM.SNR.(SubSpec{ss})(kk)    = op_getSNR(MRSCont.processed.(SubSpec{ss}){kk});       
            MRSCont.QM.FWHM.(SubSpec{ss})(kk)   = op_getLW(MRSCont.processed.(SubSpec{ss}){kk},SNRRange{ss}(1),SNRRange{ss}(2)); % in Hz       
            if ~(strcmp(SubSpec{ss},'ref') || strcmp(SubSpec{ss},'w'))
                MRSCont.QM.freqShift.(SubSpec{ss})(kk)  = refShift_SubSpecAlign + refShift_final;       
                MRSCont.QM.res_water_amp.(SubSpec{ss})(kk) = sum(MRSCont.processed.(SubSpec{ss}){kk}.watersupp.amp);  
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