function [MRSCont] = osp_processDW(MRSCont)
%% [MRSCont] = osp_processDW(MRSCont)
%   This function processes diffusion-weighted MRS data.
%
%   USAGE:
%       [MRSCont] = osp_processDW(MRSCont);
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

warning('off','all');

% Loop over all datasets
refProcessTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end

for kk = 1:MRSCont.nDatasets
    [~] = printLog('OspreyProcess', kk, MRSCont.nDatasets, progressText, MRSCont.flags.isGUI, MRSCont.flags.isMRSI); 
    
    if ~(MRSCont.flags.didProcess == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'processed') && (kk > length(MRSCont.processed.A)))  || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)

        %%% 1. GET RAW DATA %%%
        raw                         = MRSCont.raw{kk};                     % Get the kk-th dataset
        % Determine the dimension that stores the different
        % diffusion-weighted spectra, and determine the FID-A function that
        % we need to extract the respective dimension.
        if raw.dims.subSpecs == 0
            if raw.dims.extras == 0
                error('This dataset does not appear to have sub-spectra. ABORT!');
            else
                dimDW = raw.dims.extras;
                takeFunction = @op_takeextras;
            end
        else
            dimDW = raw.dims.subSpecs;
            takeFunction = @op_takesubspec;
        end
        
        %%% 2. GET REFERENCE DATA / EDDY CURRENT CORRECTION %%%
        % If there are reference scans, load them here to allow eddy-current
        % correction of the raw data.
        if MRSCont.flags.hasRef
            raw_ref                 = MRSCont.raw_ref{kk};                 % Get the kk-th dataset
        end
        
        %%% 3. LOOP OVER THE DIFFUSION-WEIGHTING DIMENSION
        % Determine number of diffusion-weighted spectra
        nDW = raw.sz(dimDW);
        
        for dd = 1:nDW
            % Raw data
            raw_temp        = takeFunction(raw, dd);
            % Reference data
            if MRSCont.flags.hasRef
                raw_ref_temp    = takeFunction(raw_ref, dd);
            end
            
            %%% 2a. PHANTOM-SPECIFIC PRE-PROCESSING %%%
            % If this is phantom data (assuming room temperature), we want to
            % perform a few specific pre-processing steps.
            if MRSCont.flags.isPhantom
                % Next, shift the entire metabolite spectrum by 0.15 ppm.
                % This doesn't have to be completely accurate, since additional
                % referencing steps are performed in the later stages of
                % post-processing and modelling, but we want the prominent singlets
                % to appear within 0.1 ppm of their expected in-vivo positions.
                phantomShiftPPM = 0.15 * raw_temp.txfrq*1e-6;
                raw_temp = op_freqshift(raw_temp, -phantomShiftPPM);
                
                % Finally, apply some linebroadening. High-quality in-vitro
                % data may have linewidth lower than the simulated basis set
                % data.
                raw_temp = op_filter(raw_temp, 2);
                if MRSCont.flags.hasRef
                    raw_ref_temp = op_filter(raw_ref_temp, 2);
                end
            end
        
            
            %%% 3. FREQUENCY/PHASE CORRECTION AND AVERAGING %%%
            if raw_temp.averages > 1 && raw_temp.flags.averaged == 0
                % Calculate starting values for spectral registration
                [refShift_ind_ini] = op_preref(raw_temp);
                if ~MRSCont.flags.isPhantom
                    switch MRSCont.opts.SpecReg %Pick spectral registration method (default is Robust Spectral Registration)
                        case 'RobSpecReg'
                            [raw_temp, fs, phs, weights, driftPre, driftPost]     = op_robustSpecReg(raw_temp, 'unedited', 0, refShift_ind_ini); % Align and average
                        case 'RestrSpecReg'
                            [raw_temp, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw_temp, 'unedited', 0, refShift_ind_ini,0, MRSCont.opts.fit.range(1), MRSCont.opts.fit.range(2)); % Align and average
                        case 'none'
                            [raw_temp, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw_temp, 'unedited', 0, refShift_ind_ini,1); % Align and average
                    end
                else
                    [raw_temp, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw_temp, 'unedited', 0, refShift_ind_ini, 0 , 0.5, 4.2); % Align and average
                end
                raw_temp.specReg.fs              = fs; % save align parameters
                raw_temp.specReg.phs             = phs; % save align parameters
                raw_temp.specReg.weights         = weights{1}(1,:)'; % save align parameters);
                raw_temp.specReg.weights         = raw_temp.specReg.weights/max(raw_temp.specReg.weights);
            else
                raw_temp.flags.averaged  = 1;
                raw_temp.dims.averages   = 0;
                raw_temp.specReg.fs              = 0; % save align parameters
                raw_temp.specReg.phs             = 0; % save align parameters
                raw_temp.specReg.weights         = 1; % save align parameters
                driftPre = op_measureDrift(raw_temp);
                driftPost = driftPre;
            end
            
            %%% 4. EDDY CURRENT CORRECTION %%%
            % If there are reference scans, do eddy-current correction.
            if MRSCont.flags.hasRef
                
                % Align and average (if more than one transient)
                if raw_ref_temp.averages > 1 && raw_ref_temp.flags.averaged == 0
                    [raw_ref_temp,~,~]  = op_alignAverages(raw_ref_temp, 1, 'n');
                    raw_ref_temp        = op_averaging(raw_ref_temp);            % Average
                else
                    raw_ref_temp.flags.averaged  = 1;
                    raw_ref_temp.dims.averages   = 0;
                end
                
                [raw_temp,raw_ref_temp] = op_eccKlose(raw_temp, raw_ref_temp);   % Klose eddy current correction
                [raw_ref_temp,~]        = op_ppmref(raw_ref_temp,4.6,4.8,4.68);  % Reference to water @ 4.68 ppm
                
                % Save back to MRSCont
                MRSCont.processed.ref{kk}{dd}     = raw_temp;
            
            end
        
            %%% 5. REMOVE RESIDUAL WATER %%%
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
            raw_temp = op_iterativeWaterFilter(raw_temp, waterRemovalFreqRange, 32, fracFID*length(raw_temp.fids), 0);
        
            % Save back to MRSCont
            MRSCont.processed.A{kk}{dd}     = raw_temp;
        end
        

        %%% 9. QUALITY CONTROL PARAMETERS %%%
        SubSpec = {'A'};       
        SNRRange = {[1.8,2.2]};
        if MRSCont.flags.hasMM
            SubSpec{end+1} = 'mm';
            SNRRange{end+1} = [0.7,1.1];
        end
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
            for dd = 1:nDW
                MRSCont.QM.SNR.(SubSpec{ss})(kk,dd)    = op_getSNR(MRSCont.processed.(SubSpec{ss}){kk}{dd});
                MRSCont.QM.FWHM.(SubSpec{ss})(kk,dd)   = op_getLW(MRSCont.processed.(SubSpec{ss}){kk}{dd},SNRRange{ss}(1),SNRRange{ss}(2)); % in Hz
                if ~(strcmp(SubSpec{ss},'ref') || strcmp(SubSpec{ss},'w') || strcmp(SubSpec{ss},'mm'))
                    MRSCont.QM.drift.pre.(SubSpec{ss}){kk}{dd}  = driftPre;
                    MRSCont.QM.drift.post.(SubSpec{ss}){kk}{dd} = driftPost;
                    MRSCont.QM.res_water_amp.(SubSpec{ss})(kk,dd) = sum(MRSCont.processed.(SubSpec{ss}){kk}{dd}.watersupp.amp);
                    if strcmp(SubSpec{ss},'diff1') || strcmp(SubSpec{ss},'sum')
                        MRSCont.QM.drift.pre.(SubSpec{ss}){kk}{dd}  = reshape([MRSCont.QM.drift.pre.A'; MRSCont.QM.drift.pre.B'], [], 1)';
                        MRSCont.QM.drift.post.(SubSpec{ss}){kk}{dd} = reshape([MRSCont.QM.drift.post.A'; MRSCont.QM.drift.post.B'], [], 1)';
                    end
                    MRSCont.QM.drift.pre.AvgDeltaCr.(SubSpec{ss})(kk,dd) = mean(MRSCont.QM.drift.pre.(SubSpec{ss}){kk}{dd} - 3.02);
                    MRSCont.QM.drift.post.AvgDeltaCr.(SubSpec{ss})(kk,dd) = mean(MRSCont.QM.drift.pre.(SubSpec{ss}){kk}{dd} - 3.02);
                end
            end
        end              
    end
end

time = toc(refProcessTime);
[~] = printLog('done', time, MRSCont.nDatasets, progressText, MRSCont.flags.isGUI, MRSCont.flags.isMRSI); 

%%% 10. SET FLAGS %%%
MRSCont.flags.avgsAligned       = 1;
MRSCont.flags.averaged          = 1;
MRSCont.flags.ECCed             = 1;
MRSCont.flags.waterRemoved      = 1;
MRSCont.runtime.Proc = time;
% Close any remaining open figures
close all;

end
