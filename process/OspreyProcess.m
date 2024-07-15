function [MRSCont] = OspreyProcess(MRSCont)
%% [MRSCont] = OspreyProcess(MRSCont)
%   This function pre-processes MRS data from all major vendors.
%   Data is read from the provided input filenames. It is shaped,
%   preprocessed, aligned, etc. according to the type of sequence
%   (un-edited data, MEGA-edited (ON/OFF), HERMES/HERCULES (A/B/C/D),
%   etc.).
%
%   USAGE:
%       MRSCont = OspreyProcess(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-19)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-02-19: First version of the code.

outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));

% Checking for version, toolbox, and previously run modules
[~,MRSCont.ver.CheckOsp ] = osp_CheckRunPreviousModule(MRSCont, 'OspreyProcess');


%% Data post-processing starts here
warning('off','all');

% Loop over all datasets
refProcessTime = tic;
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
else
    progressText = '';
end

for kk = 1:MRSCont.nDatasets(1) %Subject loop
    for ll = 1: 1:MRSCont.nDatasets(2) %Experiment loop
        [~] = printLog('OspreyProcess',kk,ll,MRSCont.nDatasets(1),progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);


        if ~(MRSCont.flags.didProcess == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'processed') && (kk > size(MRSCont.processed.metab,2)))  || ~strcmp(MRSCont.ver.Osp,MRSCont.ver.CheckOsp)
            metab_ll = MRSCont.opts.MultipleSpectra.metab(ll);
%%          %%% 1. GET RAW DATA %%%
            raw                         = MRSCont.raw{metab_ll,kk};                                          % Get the kk-th dataset

            % Get sequence type set up
            if raw.flags.isUnEdited
                seq = 'unedited';
            elseif raw.flags.isMEGA
                seq = 'MEGA';
            elseif raw.flags.isHERMES
                seq = 'HERMES';
            elseif raw.flags.isHERCULES
                seq = 'HERCULES';
            end

%%          %%% 1B. GET MM DATA %%%
            if MRSCont.flags.hasMM
                mm_ll = MRSCont.opts.MultipleSpectra.mm(ll);
                raw_mm                         = MRSCont.raw_mm{mm_ll,kk};              % Get the kk-th dataset re_mm
                if raw_mm.averages > 1 && raw_mm.flags.averaged == 0 %re_mm
                    [raw_mm,fs,phs]               = op_alignAverages(raw_mm, 1, 'n'); %re_mm
                    raw_mm                     = op_averaging(raw_mm);            % Average re_mm
                else %re_mm
                    raw_mm.flags.averaged  = 1; %re_mm
                    raw_mm.dims.averages   = 0; %re_mm
                    fs = 0;
                    phs = 0;
                end
                if raw_mm.flags.isUnEdited
                    [raw_mm,~]                     = op_ppmref(raw_mm,4.6,4.8,4.68);  % Reference to water @ 4.68 ppm  %re_mm
                end
                if raw_mm.flags.isUnEdited
                    raw_mm.names = {'A'};
                    raw_mm.target = {''};
                    raw_mm.specReg{1}.fs              = fs; % save align parameters
                    raw_mm.specReg{1}.phs             = phs; % save align parameters
                    raw_mm.specReg{1}.weights{1}         = ones(size(phs)); % save align parameters
                    driftPre = 0;
                end
                if raw_mm.flags.isMEGA
                    raw_mm.names = {'A','B'};
                    raw_mm.target = MRSCont.opts.editTarget';
                    raw_mm.specReg{1}.fs              = fs; % save align parameters
                    raw_mm.specReg{1}.phs             = phs; % save align parameters
                    raw_mm.specReg{1}.weights{1}         = ones(size(phs)); % save align parameters
                    raw_mm.specReg{1}.weights{2}         = ones(size(phs)); % save align parameters
                    driftPre{1} = 0;
                    driftPre{2} = 0;
                end
                if raw_mm.flags.isHERMES || raw_mm.flags.isHERCULES
                    raw_mm.names = {'A', 'B', 'C', 'D'};
                    raw_mm.target = MRSCont.opts.editTarget';
                    raw_mm.specReg{1}.fs              = fs; % save align parameters
                    raw_mm.specReg{1}.phs             = phs; % save align parameters
                    raw_mm.specReg{1}.weights{1}         = size(phs); % save align parameters
                    raw_mm.specReg{1}.weights{2}         = size(phs); % save align parameters
                    raw_mm.specReg{1}.weights{3}         = size(phs); % save align parameters
                    raw_mm.specReg{1}.weights{4}         = size(phs); % save align parameters
                    driftPre{1} = 0;
                    driftPre{2} = 0;
                    driftPre{3} = 0;
                    driftPre{4} = 0;
                end
            end  %re_mm

%%          %%% 1C. GET MM REFERENCE DATA %%%
            if MRSCont.flags.hasMMRef
                mm_ref_ll = MRSCont.opts.MultipleSpectra.mm_ref(ll);
                raw_mm_ref = MRSCont.raw_mm_ref{mm_ref_ll,kk};              % Get the kk-th dataset
                raw_mm_ref = op_combine_water_subspecs(raw_mm_ref);
            end
            
%%          %%% 1D. GET REFERENCE DATA %%%
            if MRSCont.flags.hasRef
                ref_ll = MRSCont.opts.MultipleSpectra.ref(ll);
                raw_ref = MRSCont.raw_ref{ref_ll,kk};              % Get the kk-th dataset
                
                % Combine SPECIAL sub-spectra
                if MRSCont.flags.isSPECIAL
                    raw_ref = combine_special_subspecs(raw_ref);
                end

                if raw_ref.averages > 1 && raw_ref.flags.averaged == 0
                    [raw_ref,fs,phs]                 = op_alignAverages(raw_ref, 1, 'n');
                    badAverages=find(abs(fs)>6);
                    goodAverages=find(abs(fs)<6); 
                    if ~isempty(badAverages)
                        fs(badAverages) =[];
                        phs(badAverages) =[];
                        raw_ref.AveragesToKeep = goodAverages;
                        raw_ref.fids=raw_ref.fids(:,goodAverages,:,:);
                        raw_ref.specs=fftshift(ifft(raw_ref.fids,[],raw_ref.dims.t),raw_ref.dims.t);
                        raw_ref.sz=size(raw_ref.fids);
                        raw_ref.averages=length(goodAverages);
                    end
                    raw_ref                       = op_averaging(raw_ref);              % Average
                    raw_ref.specReg{1}.fs              = fs; % save align parameters
                    raw_ref.specReg{1}.phs             = phs; % save align parameters
                else
                    raw_ref.flags.averaged    = 1;
                    raw_ref.dims.averages     = 0;
                end
            
                raw_ref = op_combine_water_subspecs(raw_ref);                
            end

%%           %%% 1E. GET SHORT-TE WATER DATA %%%
             % And do the processing
             if MRSCont.flags.hasWater
                w_ll = MRSCont.opts.MultipleSpectra.w(ll);
                raw_w                           = MRSCont.raw_w{w_ll,kk};                % Get the kk-th dataset
                
                % Combine SPECIAL sub-spectra
                if MRSCont.flags.isSPECIAL
                    raw_w = combine_special_subspecs(raw_w);
                end
                
                if raw_w.averages > 1 && raw_w.flags.averaged == 0
                    [raw_w,fs,phs]                 = op_alignAverages(raw_w, 1, 'n');
                    badAverages=find(abs(fs)>6);
                    goodAverages=find(abs(fs)<6); 
                    if ~isempty(badAverages)
                        fs(badAverages) =[];
                        phs(badAverages) =[];
                        raw_w.AveragesToKeep = goodAverages;
                        raw_w.fids=raw_w.fids(:,goodAverages,:,:);
                        raw_w.specs=fftshift(ifft(raw_w.fids,[],raw_w.dims.t),raw_w.dims.t);
                        raw_w.sz=size(raw_w.fids);
                        raw_w.averages=length(goodAverages);
                    end
                    raw_w                       = op_averaging(raw_w);              % Average
                    raw_w.specReg{1}.fs              = fs; % save align parameters
                    raw_w.specReg{1}.phs             = phs; % save align parameters
                else
                    raw_w.flags.averaged    = 1;
                    raw_w.dims.averages     = 0;
                end

                if raw_w.subspecs > 1
                    raw_w = op_combine_water_subspecs(raw_w);
                end
                if ~MRSCont.flags.isMRSI
                    [raw_w,~]                       = op_eccKlose(raw_w, raw_w);        % Klose eddy current correction
                else
                    [raw_w,~]=op_autophase(raw_w,2,2*4.68);
                end
                [raw_w,~]                       = op_ppmref(raw_w,4.6,4.8,4.68);    % Reference to water @ 4.68 ppm
                raw_w.names = {'A'};
                MRSCont.processed.w{w_ll,kk}         = raw_w; % Save back to MRSCont container
            end

%%          %%% 2a. PHANTOM-SPECIFIC PRE-PROCESSING %%%
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
                if MRSCont.flags.hasWater
                    raw_w = op_filter(raw_w, 2);
                    MRSCont.processed.w{ll,kk}         = raw_w; % Save back to MRSCont container
                end
            end

%%          %%% 2b. COMBINING SPECIAL SUB-SPECTRA %%%
            if MRSCont.flags.isSPECIAL
                raw = combine_special_subspecs(raw);
            end
        
%%          %%% 3. FREQUENCY/PHASE CORRECTION AND AVERAGING %%%
            if raw.averages > 1 && raw.flags.averaged == 0
                % Calculate starting values for spectral registration
                 [refShift_ind_ini]=op_preref(raw,seq);
                % Perform robust spectral correction with weighted averaging.
                % This can obviously only be done, if the spectra have not been
                % pre-averaged, i.e. in some older RDA and DICOM files (which should,
                % generally, not be used).
                if ~MRSCont.flags.isPhantom
                    switch MRSCont.opts.SpecReg %Pick spectral registration method (default is Robust Spectral Registration)
                        case 'ProbSpecReg'
                            [raw, fs, phs, weights, driftPre, driftPost]   = op_probabSpecReg(raw, seq, 0,refShift_ind_ini);
                        case 'RobSpecReg'
                            [raw, fs, phs, weights, driftPre, driftPost]     = op_robustSpecReg(raw, seq, 0,refShift_ind_ini); % Align and average
                        case 'RestrSpecReg'
                            if isfield(MRSCont.opts,'SpecRegRange')
                                [raw, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw, seq, 0,refShift_ind_ini,0,MRSCont.opts.SpecRegRange(1),MRSCont.opts.SpecRegRange(2)); % Align and average
                            else
                                [raw, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw, seq, 0,refShift_ind_ini,0,MRSCont.opts.fit.range(1),MRSCont.opts.fit.range(2)); % Align and average
                            end
                        case 'none'
                            [raw, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw, seq, 0,refShift_ind_ini,1); % Align and average
                    end
                else
                    switch MRSCont.opts.SpecReg %Pick spectral registration method (default is Robust Spectral Registration)
                        case 'none'
                            [raw, fs, phs, weights, driftPre, driftPost]     = op_SpecRegFreqRestrict(raw, seq, 0,refShift_ind_ini,1); % Align and average
                        otherwise
                            [raw, fs, phs, weights, driftPre, driftPost]   = op_SpecRegFreqRestrict(raw, seq, 0,refShift_ind_ini,0,0.5,4.2);
                    end
                end

                raw.specReg{1}.fs              = fs; % save align parameters
                raw.specReg{1}.phs             = phs; % save align parameters
                if raw.flags.isUnEdited
                    raw.specReg{1}.weights         = weights; % save align parameters);
                    raw.specReg{1}.weights{1}         = raw.specReg{1}.weights{1}'/max(raw.specReg{1}.weights{1}(1,:));
                end
                if raw.flags.isMEGA
                    raw.specReg{1}.weights         = weights; % save align parameters);
                    raw.specReg{1}.weights{1} = raw.specReg{1}.weights{1}'/(max(raw.specReg{1}.weights{1}(1,:)));
                    raw.specReg{1}.weights{2} = raw.specReg{1}.weights{2}'/(max(raw.specReg{1}.weights{2}(1,:)));
                end
                if raw.flags.isHERMES || raw.flags.isHERCULES
                    raw.specReg{1}.weights         = weights; % save align parameters);
                    raw.specReg{1}.weights{1} = raw.specReg{1}.weights{1}'/(max(raw.specReg{1}.weights{1}(1,:)));
                    raw.specReg{1}.weights{2} = raw.specReg{1}.weights{2}'/(max(raw.specReg{1}.weights{2}(1,:)));
                    raw.specReg{1}.weights{3} = raw.specReg{1}.weights{3}'/(max(raw.specReg{1}.weights{3}(1,:)));
                    raw.specReg{1}.weights{4} = raw.specReg{1}.weights{4}'/(max(raw.specReg{1}.weights{4}(1,:)));
                end

            else
                % If there is still an average dimension, squeeze
                if raw.dims.averages > 0
                    raw.fids = squeeze(sum(raw.fids, raw.dims.averages));
                    raw.specs=fftshift(fft(raw.fids,[],raw.dims.t),raw.dims.t);
                    if raw.dims.t>raw.dims.averages
                        raw.dims.t=raw.dims.t-1;
                    else
                        raw.dims.t=raw.dims.t;
                    end
                    if raw.dims.coils>raw.dims.averages
                        raw.dims.coils=raw.dims.coils-1;
                    else
                        raw.dims.coils=raw.dims.coils;
                    end
                    raw.dims.averages=0;
                    if raw.dims.subSpecs>raw.dims.averages
                        raw.dims.subSpecs=raw.dims.subSpecs-1;
                    else
                        raw.dims.subSpecs=raw.dims.subSpecs;
                    end
                    if raw.dims.extras>raw.dims.averages
                        raw.dims.extras=raw.dims.extras-1;
                    else
                        raw.dims.extras=raw.dims.extras;
                    end
                end
                raw.sz=size(raw.fids);

                raw.flags.averaged  = 1;
                if raw.flags.isUnEdited
                    raw.specReg{1}.fs              = 0; % save align parameters
                    raw.specReg{1}.phs             = 0; % save align parameters
                    raw.specReg{1}.weights{1}         = 1; % save align parameters
                    driftPre = op_measureDrift(raw);
                end
                if raw.flags.isMEGA
                    raw.specReg{1}.fs              = zeros(1,2); % save align parameters
                    raw.specReg{1}.phs             = zeros(1,2); % save align parameters
                    raw.specReg{1}.weights{1}         = ones(1,1); % save align parameters
                    raw.specReg{1}.weights{2}         = ones(1,1); % save align parameters
                    driftPre{1} = 0;
                    driftPre{2} = 0;
                end
                if raw.flags.isHERMES || raw.flags.isHERCULES
                    raw.specReg{1}.fs              = zeros(1,4); % save align parameters
                    raw.specReg{1}.phs             = zeros(1,4); % save align parameters
                    raw.specReg{1}.weights{1}         = ones(1,1); % save align parameters
                    raw.specReg{1}.weights{2}         = ones(1,1); % save align parameters
                    raw.specReg{1}.weights{3}         = ones(1,1); % save align parameters
                    raw.specReg{1}.weights{4}         = ones(1,1); % save align parameters
                    driftPre{1} = 0;
                    driftPre{2} = 0;
                    driftPre{3} = 0;
                    driftPre{4} = 0;
                end
                driftPost = driftPre;
            end

%%          %%% 4. GET REFERENCE DATA / EDDY CURRENT CORRECTION %%%
            % If there are reference scans, load them here to allow eddy-current
            % correction of the raw data.
            if MRSCont.flags.hasRef
                if MRSCont.flags.hasMM
                    if MRSCont.flags.hasMMRef
                        if MRSCont.opts.ECC.mm(mm_ref_ll,kk)
                            [raw_mm,raw_mm_ref]                   = op_eccKlose(raw_mm, raw_mm_ref);        % Klose eddy current correction
                        else
                            [~,raw_mm_ref]                   = op_eccKlose(raw_mm, raw_mm_ref);        % Klose eddy current correction
                        end
                        MRSCont.processed.mm_ref{mm_ref_ll,kk}       = raw_mm_ref;                          % Save back to MRSCont container
                    else
                        if MRSCont.opts.ECC.mm(mm_ref_ll,kk)
                            [raw_mm,~]                   = op_eccKlose(raw_mm, raw_ref);        % Klose eddy current correction
                        else
                            [~,raw_mm_ref]                   = op_eccKlose(raw_mm, raw_ref);        % Klose eddy current correction
                        end
                    end
                end
                if MRSCont.opts.ECC.raw(metab_ll,kk)
                    [raw,raw_ref]                   = op_eccKlose(raw, raw_ref);        % Klose eddy current correction
                else
                    try
                        [~,raw_ref]                   = op_eccKlose(raw, raw_ref);        % Klose eddy current correction
                    catch
                        [~,raw_ref]                   = op_eccKlose(raw_ref, raw_ref);        % Klose eddy current correction
                    end
                end
                [raw_ref,~]                     = op_ppmref(raw_ref,4.6,4.8,4.68);  % Reference to water @ 4.68 ppm
                MRSCont.processed.ref{metab_ll,kk}       = raw_ref;                          % Save back to MRSCont container
            end
            
%%          %%% 5. DETERMINE POLARITY OF SPECTRUM (EG FOR MOIST WATER SUPP) %%%
            % Automate determination whether the Cr peak has positive polarity.
            % For water suppression methods like MOIST, the residual water may
            % actually have negative polarity, but end up positive in the data, so
            % that the spectrum needs to be flipped.
            if ~isfield(MRSCont.opts.SubSpecAlignment, 'polResidCr')
                if ~isfield(MRSCont.opts.SubSpecAlignment, 'PreservePolarity') || ...
                    (isfield(MRSCont.opts.SubSpecAlignment, 'PreservePolarity') && ~MRSCont.opts.SubSpecAlignment.PreservePolarity)
                    raw_Cr     = op_freqrange(raw,2.8,3.2);
                    % Determine the polarity of the respective peak: if the absolute of the
                    % maximum minus the absolute of the minimum is positive, the polarity
                    % of the respective peak is positive; if the absolute of the maximum
                    % minus the absolute of the minimum is negative, the polarity is negative.
                    polResidCr = abs(max(real(raw_Cr.specs))) - abs(min(real(raw_Cr.specs)));
                    polResidCr = squeeze(polResidCr);
                    polResidCr(polResidCr<0) = -1;
                    polResidCr(polResidCr>0) = 1;
                else
                    polResidCr = 1;
                end
            else
                polResidCr = MRSCont.opts.SubSpecAlignment.polResidCr;
            end
            raw = op_ampScale(raw,polResidCr);
            MRSCont.raw{metab_ll,kk} = op_ampScale(MRSCont.raw{metab_ll,kk},polResidCr);

            % Do the same for the metabolite-nulled data
            if MRSCont.flags.hasMM
                raw_Cr39    = op_freqrange(raw_mm,3.7,4.1);
                polResidCr39  = abs(max(real(raw_Cr39.specs))) - abs(min(real(raw_Cr39.specs)));
                polResidCr39(polResidCr39<0) = -1;
                polResidCr39(polResidCr39>0) = 1;
                if polResidCr39(1) ~= polResidCr(1)
                    raw_mm = op_ampScale(raw_mm,polResidCr);
                else
                    raw_mm = op_ampScale(raw_mm,polResidCr39);
                end

            end

%%          %%% 6. DETERMINE ON/OFF STATUS FOR EDITED DATA & NAMES
            % Classify the MEGA sub-spectra such that the OFF spectrum is stored in
            % the first entry , and the ON spectrum is the second entry along the
            % subspectra dimension. For HERMES and HERCULES A,B,C,D are stored
            % as first, second, third, and fourth entry, respectively.
            % Parse input arguments
            if raw.flags.isUnEdited
                raw.names = {'A'};
                MRSCont.QM.drift.pre.A{kk} = driftPre;
                MRSCont.QM.drift.post.A{kk} = driftPost;
                raw.target = {''};
            end
            if raw.flags.isMEGA
                target = MRSCont.opts.editTarget{1}; % GABA editing as default
                if isfield(MRSCont.opts, 'Order')
                    Order = MRSCont.opts.Order;
                else
                    Order = [];
                end
                [raw, switchOrder]  = osp_onOffClassifyMEGA(raw, target,Order);
                raw.names = {'A','B'};
                raw.target = MRSCont.opts.editTarget';

                if switchOrder
                    raw.flags.orderswitched = 1;
                else
                    raw.flags.orderswitched = 0;
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
                    raw.specReg{1}.fs              = raw.specReg{1}.fs(:,[2 1]); % save align parameters
                    raw.specReg{1}.phs             = raw.specReg{1}.phs(:,[2 1]); % save align parameters
                    raw.specReg{1}.weights         = raw.specReg{1}.weights([2 1]); % save align parameters);
                end
            end
            if raw.flags.isHERMES || raw.flags.isHERCULES
                if (strcmp(MRSCont.opts.editTarget{1},'HERCULES1')||strcmp(MRSCont.opts.editTarget{1},'HERCULES2'))
                    MRSCont.opts.editTarget = {'GABA','GSH'};
                end
                target1 = MRSCont.opts.editTarget{1};
                target2 = MRSCont.opts.editTarget{2};
                if length(MRSCont.opts.editTarget) > 2
                   target3 = MRSCont.opts.editTarget{3};
                end
                if isfield(MRSCont.opts, 'Order')
                    Order = MRSCont.opts.Order;
                else
                    Order = [];
                end
                [raw, commuteOrder] = osp_onOffClassifyHERMES(raw,[target1 target2],Order);
                raw.commuteOrder = commuteOrder;
                raw.names = {'A', 'B', 'C', 'D'};
                raw.target = MRSCont.opts.editTarget';
                driftPre = driftPre(:,commuteOrder);
                driftPost = driftPost(:,commuteOrder);
                raw.specReg{1}.fs = raw.specReg{1}.fs(:,commuteOrder);
                raw.specReg{1}.phs = raw.specReg{1}.phs(:,commuteOrder);
                raw.specReg{1}.weights = raw.specReg{1}.weights(:,commuteOrder);
                % Save drift information back to container
                MRSCont.QM.drift.pre.A{kk} = driftPre{1};
                MRSCont.QM.drift.pre.B{kk} = driftPre{2};
                MRSCont.QM.drift.pre.C{kk} = driftPre{3};
                MRSCont.QM.drift.pre.D{kk} = driftPre{4};
                MRSCont.QM.drift.post.A{kk} = driftPost{1};
                MRSCont.QM.drift.post.B{kk} = driftPost{2};
                MRSCont.QM.drift.post.C{kk} = driftPost{3};
                MRSCont.QM.drift.post.D{kk} = driftPost{4};
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
                if exist('target3', 'var') % For HERMES 4 acqusitions
                    MRSCont.QM.drift.pre.diff3{kk}  = driftPre;
                end
                MRSCont.QM.drift.pre.sum{kk}    = driftPre;
                driftPost = [MRSCont.QM.drift.post.A{kk}, MRSCont.QM.drift.post.B{kk}, MRSCont.QM.drift.post.C{kk}, MRSCont.QM.drift.post.D{kk}]';
                try
                    driftPost = reshape(driftPost, [raw.averages, 1]);
                catch
                    driftPost = reshape(driftPost, [raw.rawAverages, 1]);
                end
                MRSCont.QM.drift.post.diff1{kk}  = driftPost;
                MRSCont.QM.drift.post.diff2{kk}  = driftPost;
                if exist('target3', 'var') % For HERMES 4 acqusitions
                    MRSCont.QM.drift.post.diff3{kk}  = driftPre;
                end
                MRSCont.QM.drift.post.sum{kk}    = driftPost;
            end

%%          %%% 6b. BUILD SUM AND DIFF SPECTRA %%%
            % Correct the frequency axis so that Cr appears at 3.027 ppm
            if raw.flags.isMEGA
                % Create all Hadamard combinations to generate a well defined
                % raw struct
                raw=op_HadamardScans(raw,[-1 1],'diff1');
                raw=op_HadamardScans(raw,[1 1],'sum');
                raw_no_subspec_aling = raw;
                if ~raw.flags.isMRSI
%                     [raw,~]       = op_phaseCrCho(raw, 1);
                    [refShift_SubSpecAlign, ~] = osp_XReferencing(raw,[3.03 3.22],[1 1],[1.85 4.2]);% determine frequency shift
                    if abs(refShift_SubSpecAlign) > 10 % This a huge shift. Most likley wrong and we will try it again with tNAA only
                        [refShift_SubSpecAlign, ~] = osp_XReferencing(raw,2.01,1,[1.85 4.2]);% determine frequency shift
                    end
%                     % Apply initial referencing shift
                    raw = op_freqshift(raw, -refShift_SubSpecAlign);
%                     % Fit a double-Lorentzian to the Cr-Cho area, and phase the spectrum
%                     % with the negative phase of that fit
                    [raw,ph]       = op_phaseCrCho(raw, 1);

                else
                    refShift_SubSpecAlign =0;
                end

                switch MRSCont.opts.SubSpecAlignment.mets
                    case 'L1Norm'
                        [raw]  = osp_editSubSpecAlignLNorm(raw,seq);
                    case 'L2Norm'
                        [raw]  = osp_editSubSpecAlign(raw, seq, target,MRSCont.opts.UnstableWater);
                    otherwise
                end

                if MRSCont.flags.hasMM
                    raw_mm=op_HadamardScans(raw_mm,[-1 1],'diff1');
%                     diff1_mm = op_takesubspec(raw_mm,3);
%                     [diff1_mm,~]          = op_autophase(diff1_mm,0.5,1.1);
%                     raw_mm = op_mergesubspec(raw_mm,diff1_mm);
%                     [refShift_mm, ~] = fit_OspreyReferencingMM(raw_mm);
%                     [raw_mm]             = op_freqshift(raw_mm,-refShift_mm);
                    switch MRSCont.opts.SubSpecAlignment.mm
                        case 'L1Norm'
                        [raw_mm]  = osp_editSubSpecAlignLNorm(raw_mm,seq);
                        case 'L2Norm'
                        [raw_mm]  = osp_editSubSpecAlign(raw_mm, seq, target,MRSCont.opts.UnstableWater);
                        otherwise
                    end
                    if switchOrder
                        raw_mm.flags.orderswitched = 1;
                    else
                        raw_mm.flags.orderswitched = 0;
                    end
                    raw_mm=op_HadamardScans(raw_mm,[-1 1],'diff1');
                    raw_mm=op_HadamardScans(raw_mm,[1 1],'sum');
                end
                % Create the edited difference spectrum
                raw=op_HadamardScans(raw,[-1 1],'diff1');
                % Create the sum spectrum
                raw=op_HadamardScans(raw,[1 1],'sum');
                raw.target = target;
                raw = op_freqshift(raw, refShift_SubSpecAlign);
            end

            if raw.flags.isHERMES || raw.flags.isHERCULES
                % Create all Hadamard combinations to generate a well defined
                % raw struct
                if ~exist('target3', 'var')
                    raw=op_HadamardScans(raw,[-1 1 -1 1],'diff1');
                    raw=op_HadamardScans(raw,[-1 -1 1 1],'diff2');
                else  % For HERMES 4 acqusitions
                    raw=op_HadamardScans(raw,[-1 1 -1 1],'diff1');
                    raw=op_HadamardScans(raw,[-1 -1 1 1],'diff2');
                    raw=op_HadamardScans(raw,[-1 1 1 -1],'diff3');
                end
                raw=op_HadamardScans(raw,[1 1 1 1],'sum');
                raw_no_subspec_aling = raw;

                % Correct the frequency axis so that Cr appears at 3.027 ppm
                [refShift_SubSpecAlign, ~] = osp_XReferencing(raw,[3.03 3.22],[1 1],[1.85 3.9]);% determine frequency shift
                if abs(refShift_SubSpecAlign) > 10 % This a huge shift. Most likley wrong and we will try it again with tNAA only
                    [refShift_SubSpecAlign, ~] = osp_XReferencing(raw,2.01,1,[1.85 3.9]);% determine frequency shift
                end

                % Apply initial referencing shift
                raw = op_freqshift(raw, -refShift_SubSpecAlign);
                % Fit a double-Lorentzian to the Cr-Cho area, and phase the spectrum
                % with the negative phase of that fit
                [raw,~] = op_phaseCrCho(raw, 1);
                % Align the sub-spectra to one another by minimizing the difference
                % between the common 'reporter' signals.
                switch MRSCont.opts.SubSpecAlignment.mets
                    case 'L1Norm'
                        [raw]  = osp_editSubSpecAlignLNorm(raw,seq);
                    case 'L2Norm'
                        if ~exist('target3', 'var')
                            [raw]  = osp_editSubSpecAlign(raw, seq, target1,target2,MRSCont.opts.UnstableWater);
                        else
                            if strcmp(target3,'EtOH')
                                [raw]  = osp_editSubSpecAlign(raw, seq, target1,target3,MRSCont.opts.UnstableWater);
                            else
                               [raw]  = osp_editSubSpecAlign(raw, seq, target2,target3,MRSCont.opts.UnstableWater); 
                            end
                        end
                    otherwise
                end
                % Update all Hadamard combinations after subspectra alignment
                if ~exist('target3', 'var')
                    raw=op_HadamardScans(raw,[-1 1 -1 1],'diff1');
                    raw=op_HadamardScans(raw,[-1 -1 1 1],'diff2');
                else  % For HERMES 4 acqusitions
                    raw=op_HadamardScans(raw,[-1 1 -1 1],'diff1');
                    raw=op_HadamardScans(raw,[-1 -1 1 1],'diff2');
                    raw=op_HadamardScans(raw,[-1 1 1 -1],'diff3');
                end
                raw=op_HadamardScans(raw,[1 1 1 1],'sum');
                raw = op_freqshift(raw, refShift_SubSpecAlign);
                fields.Method   = 'Addition of sub-spectra';
                if raw.flags.isHERCULES
                    fields.Details  = ['HERCULES Hadmard combination (Oeltzschner et al. 2019)'];
                else
                    fields.Details  = ['HERMES Hadmard combination (Saleh et al. 2018)'];
                end
                raw = op_add_analysis_provenance(raw,fields);
                
            end
            
%%          %%% 7. REMOVE RESIDUAL WATER %%%
            % Define different water removal frequency ranges, depending on
            % whether this is phantom data
            if MRSCont.flags.isPhantom
                waterRemovalFreqRange = [4.5 5];
                fracFID = 0.2;
            else
                waterRemovalFreqRange = [4.2 4.9];
                fracFID = 0.75;
            end
            % Apply iterative water filter
            raw = op_iterativeWaterFilter(raw, waterRemovalFreqRange, 32, fracFID*length(raw.fids), 0);
            if exist('raw_no_subspec_aling','var')
                raw_no_subspec_aling = op_iterativeWaterFilter(raw_no_subspec_aling, waterRemovalFreqRange, 32, fracFID*length(raw.fids), 0);
            end

            if MRSCont.flags.hasMM %re_mm
                raw_mm = op_iterativeWaterFilter(raw_mm, waterRemovalFreqRange, 32, fracFID*length(raw_mm.fids), 0);
            end

%%          %%% 8. REFERENCE SPECTRUM CORRECTLY TO FREQUENCY AXIS AND PHASE SIEMENS
            %%% DATA
            if MRSCont.flags.isPhantom
                [refShift, ~, ~] = osp_PhantomReferencing(raw);
            else
                [refShift, ~] = osp_XReferencing(raw,[3.03 3.22],[1 1],[1.85 4.2]);% determine frequency shift
                if abs(refShift) > 10 % This a huge shift. Most likley wrong and we will try it again with tNAA only
                    [refShift, ~] = osp_XReferencing(raw,2.01,1,[1.85 4.2]);% determine frequency shift
                end
            end
            [raw]             = op_freqshift(raw,-refShift);            % Reference spectra by cross-correlation
            if exist('raw_no_subspec_aling','var')
                [raw_no_subspec_aling]             = op_freqshift(raw_no_subspec_aling,-refShift);            % Reference spectra by cross-correlation
            end

            if MRSCont.flags.hasMM %re_mm
                if raw_mm.flags.isMEGA
                    temp = op_takesubspec(raw_mm,3);
                    [temp,~]          = op_autophase(temp,0.5,1.1);
                    raw_mm = op_mergesubspec(raw_mm,temp);
                end
                [refShift_mm, ~] = fit_OspreyReferencingMM(raw_mm);
                [raw_mm]             = op_freqshift(raw_mm,-refShift_mm);            % Reference spectra by cross-correlation
            end

            % Save back to MRSCont container
            if (strcmp(MRSCont.vendor,'Siemens') && raw.flags.isUnEdited && ~MRSCont.flags.isPhantom) || MRSCont.flags.isMRSI
                % Fit a double-Lorentzian to the Cr-Cho area, and phase the spectrum
                % with the negative phase of that fit
                [raw,globalPhase]       = op_phaseCrCho(raw, 1);
                raw.specReg{1}.phs = raw.specReg{1}.phs - globalPhase*180/pi;
            end

%%          %%% 9. QUALITY CONTROL PARAMETERS %%%
            SubSpec = raw.names;
            % Calculate some spectral quality metrics here;
            [raw,SNR] = op_get_Multispectra_SNR(raw);
            FWHM = op_get_Multispectra_LW(raw);
            MRSCont.processed.metab{metab_ll,kk}     = raw;
            if exist('raw_no_subspec_aling','var')
                MRSCont.processed_no_align.metab{metab_ll,kk}     = raw_no_subspec_aling;
            end
            for ss = 1 : length(SubSpec)
                MRSCont.QM.SNR.metab(metab_ll,kk,ss)    = SNR{ss};
                MRSCont.QM.FWHM.metab(metab_ll,kk,ss)   = FWHM(ss); % in Hz
                MRSCont.QM.freqShift.metab(metab_ll,kk,ss)  = refShift;
                MRSCont.QM.res_water_amp.metab(metab_ll,kk,ss) = sum(MRSCont.processed.metab{kk}.watersupp{ss}.amp);
                if strcmp(SubSpec{ss},'diff1') ||strcmp(SubSpec{ss},'diff2') || strcmp(SubSpec{ss},'diff3') ||strcmp(SubSpec{ss},'sum')
                    if raw.flags.isMEGA
                        MRSCont.QM.drift.pre.(SubSpec{ss}){metab_ll,kk}  = reshape([MRSCont.QM.drift.pre.A{kk}'; MRSCont.QM.drift.pre.B{kk}'], 1, [])';
                        MRSCont.QM.drift.post.(SubSpec{ss}){metab_ll,kk} = reshape([MRSCont.QM.drift.post.A{kk}'; MRSCont.QM.drift.post.B{kk}'], 1, [])';
                    else
                        MRSCont.QM.drift.pre.(SubSpec{ss}){metab_ll,kk}  = reshape([MRSCont.QM.drift.pre.A{kk}'; MRSCont.QM.drift.pre.B{kk}'; MRSCont.QM.drift.pre.C{kk}'; MRSCont.QM.drift.pre.D{kk}'], 1, [])';
                        MRSCont.QM.drift.post.(SubSpec{ss}){metab_ll,kk} = reshape([MRSCont.QM.drift.post.A{kk}'; MRSCont.QM.drift.post.B{kk}'; MRSCont.QM.drift.post.C{kk}'; MRSCont.QM.drift.post.D{kk}'], 1, [])';
                    end
                end
                MRSCont.QM.drift.pre.AvgDeltaCr.(SubSpec{ss})(metab_ll,kk) = mean(MRSCont.QM.drift.pre.(SubSpec{ss}){kk} - 3.02);
                MRSCont.QM.drift.post.AvgDeltaCr.(SubSpec{ss})(metab_ll,kk) = mean(MRSCont.QM.drift.post.(SubSpec{ss}){kk} - 3.02);
            end
            if MRSCont.flags.hasMM
                SubSpec = raw_mm.names;
                [raw_mm,SNR] = op_get_Multispectra_SNR(raw_mm,1);
                FWHM = op_get_Multispectra_LW(raw_mm);
                MRSCont.processed.mm{ll,kk}     = raw_mm;
                for ss = 1 : length(SubSpec)
                    MRSCont.QM.SNR.mm(mm_ll,kk,ss)    = SNR{ss};
                    MRSCont.QM.FWHM.mm(mm_ll,kk,ss)   = FWHM(ss);
                end
            end
            if MRSCont.flags.hasRef
                MRSCont.QM.SNR.ref(ref_ll,kk)    = op_getSNR(MRSCont.processed.ref{kk});
                MRSCont.QM.FWHM.ref(ref_ll,kk)   = op_getLW(MRSCont.processed.ref{kk},4.2,5.2);
                MRSCont.processed.ref{kk}.QC_names = {'water'};
            end
            if MRSCont.flags.hasWater
                MRSCont.QM.SNR.w(w_ll,kk)    = op_getSNR(MRSCont.processed.w{kk});
                MRSCont.QM.FWHM.w(w_ll,kk)   = op_getLW(MRSCont.processed.w{kk},4.2,5.2);
                MRSCont.processed.w{kk}.QC_names = {'water'};
            end
            if MRSCont.flags.hasMMRef
                MRSCont.QM.SNR.mm_ref(mm_ref_ll,kk)    = op_getSNR(MRSCont.processed.mm_ref{kk});
                MRSCont.QM.FWHM.mm_ref(mm_ref_ll,kk)   = op_getLW(MRSCont.processed.mm_ref{kk},4.2,5.2);
                MRSCont.processed.mm_ref{kk}.QC_names = {'water'};
            end
        end
    end
end


time = toc(refProcessTime);
[~] = printLog('done',time,MRSCont.nDatasets(1),MRSCont.nDatasets(2),progressText,MRSCont.flags.isGUI ,MRSCont.flags.isMRSI);

%% %%% 10. SET FLAGS %%%
MRSCont.flags.avgsAligned       = 1;
MRSCont.flags.averaged          = 1;
MRSCont.flags.ECCed             = 1;
MRSCont.flags.waterRemoved      = 1;
MRSCont.runtime.Proc = time;
% Close any remaining open figures
close all;

% Gather some more information from the processed data;
SubSpecNames = fieldnames(MRSCont.processed);
NoSubSpec = length(fieldnames(MRSCont.processed));
for ss = 1 : NoSubSpec
    for kk = 1 : MRSCont.nDatasets
            temp_sz(1,kk)= MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1);
            temp_sz_sw{1,kk} = ['np_sw_' num2str(round(MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1))) '_' num2str(round(MRSCont.processed.(SubSpecNames{ss}){1,kk}.spectralwidth))];
    end
    [MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint_spectralwidth,MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint_spectralwidth_ind,~]  = unique(temp_sz_sw,'Stable');
    [MRSCont.info.(SubSpecNames{ss}).max_ndatapoint,MRSCont.info.(SubSpecNames{ss}).max_ndatapoint_ind] = max(temp_sz);
end

%% %%% 11. COMBINE MULTI SPECTRA %%%
if MRSCont.nDatasets(2) > 1
    for ss = 1 : NoSubSpec
        for kk = 1:MRSCont.nDatasets(1)
            if size(MRSCont.processed.(SubSpecNames{ss}),1) > 1
                if ~isfield(MRSCont.processed.(SubSpecNames{ss}){1,kk}, 'extra_names')
                    MRSCont.processed.(SubSpecNames{ss}){1,kk}.extra_names{1}= MRSCont.opts.extras.names{1};
                end
                for ll = 2:MRSCont.nDatasets(2)
                    if isfield(MRSCont.processed.(SubSpecNames{ss}){ll,kk}, 'extra_names')
                        MRSCont.processed.(SubSpecNames{ss}){1,kk}= op_mergeextra(MRSCont.processed.(SubSpecNames{ss}){1,kk},MRSCont.processed.(SubSpecNames{ss}){ll,kk},MRSCont.processed.(SubSpecNames{ss}){ll,kk}.extra_names{1});
                    else
                        MRSCont.processed.(SubSpecNames{ss}){1,kk}= op_mergeextra(MRSCont.processed.(SubSpecNames{ss}){1,kk},MRSCont.processed.(SubSpecNames{ss}){ll,kk},MRSCont.opts.extras.names{ll});
                    end
                end
            end
        end
        MRSCont.processed.(SubSpecNames{ss})(2:end,:) = [];
    end
end

%% If DualVoxel or MRSI we want to extract y-axis scaling
% Creates y-axis range to align the process plots between datasets

if MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI
    MRSCont.plot.processed.match = 1; % Scaling between datasets is turned off by default
else
    MRSCont.plot.processed.match = 0; % Scaling between datasets is turned off by default
end
MRSCont = osp_scale_yaxis(MRSCont,'OspreyLoad');
MRSCont = osp_scale_yaxis(MRSCont,'OspreyProcess');
%% Clean up and save
% Set exit flags and reorder fields
MRSCont.flags.didProcess    = 1;
diary off

% Calculate QM 
MRSCont = calculate_QM(MRSCont);

% Write .tsv file and .json sidecar
osp_WriteBIDsTable(MRSCont.QM.tables, [outputFolder filesep 'QM_processed_spectra'])

% Optional:  Create all pdf figures
if MRSCont.opts.savePDF
    osp_plotAllPDF(MRSCont, 'OspreyProcess')
end

% Optional: write edited files to LCModel .RAW files
if MRSCont.opts.saveLCM && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    [MRSCont] = osp_saveLCM(MRSCont);
end

% Optional: write edited files to jMRUI .txt files
if MRSCont.opts.savejMRUI && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    [MRSCont] = osp_saveJMRUI(MRSCont);
end

% Optional: write edited files to vendor specific format files readable to
% LCModel and jMRUI
% SPAR/SDAT if Philips
% RDA if Siemens
if MRSCont.opts.saveVendor && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    [MRSCont] = osp_saveVendor(MRSCont);
end

% Optional: write edited files to NIfTI-MRS format
if MRSCont.opts.saveNII && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    [MRSCont] = osp_saveNII(MRSCont);
end

% Save the output structure to the output folder
% Determine output folder
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end



%% Functions for processing
function [raw] = combine_special_subspecs(raw)
% For SPECIAL-localized data, we adopt the pipeline from 
% https://github.com/CIC-methods/FID-A/blob/master/exampleRunScripts/run_specialproc_auto.m
[raw_ai, fs_temp, phs_temp] = op_alignAverages(raw, 0.4, 'y');
% fs_ai  = fs_temp;
% phs_ai = phs_temp;
[raw_ai, fs_temp, phs_temp] = op_alignISIS(raw_ai, 0.4);
% fs_ai(:,2)  = fs_ai(:,2) + fs_temp;
% phs_ai(:,2) = phs_ai(:,2) + phs_temp;
[raw_ai, fs_temp, phs_temp] = op_alignAverages(raw_ai, 0.4, 'y');
% fs_ai  = fs_ai + fs_temp;
% phs_ai = phs_ai + phs_temp;
[raw_ai, fs_temp, phs_temp] = op_alignISIS(raw_ai, 0.4);
% fs_ai(:,2)  = fs_ai(:,2) + fs_temp;
% phs_ai(:,2) = phs_ai(:,2) + phs_temp;

%Now combine the SPECIAL subspecs
raw = op_combinesubspecs(raw_ai, 'diff');

end

function MRSCont = calculate_QM(MRSCont)
% Sort processed field
MRSCont = osp_orderProcessFields(MRSCont);
% Store data quality measures in csv file
if MRSCont.flags.isUnEdited
    subspec = 1;
    name = 'A';
elseif MRSCont.flags.isMEGA
    subspec = 1;
    name = 'A';
elseif MRSCont.flags.isHERMES
    subspec = 7;
    name = 'sum';
elseif MRSCont.flags.isHERCULES
    subspec = 7;
    name = 'sum';
else
    msg = 'No flag set for sequence type!';
    fprintf(fileID,msg);
    error(msg);
end
names = {'Cr_SNR','Cr_FWHM','residual_water_ampl','freqShift'};
if MRSCont.flags.hasRef
    names = {'Cr_SNR','Cr_FWHM','water_FWHM','residual_water_ampl','freqShift'};
end

if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    if ~MRSCont.flags.hasRef
        QM = horzcat(MRSCont.QM.SNR.metab(1,:,subspec)',MRSCont.QM.FWHM.metab(1,:,subspec)',MRSCont.QM.res_water_amp.metab(1,:,subspec)',MRSCont.QM.freqShift.metab(1,:,subspec)');
    else
        QM = horzcat(MRSCont.QM.SNR.metab(1,:,subspec)',MRSCont.QM.FWHM.metab(1,:,subspec)',MRSCont.QM.FWHM.ref(1,:)',MRSCont.QM.res_water_amp.metab(1,:,subspec)',MRSCont.QM.freqShift.metab(1,:,subspec)');
    end
    MRSCont.QM.tables = array2table(QM,'VariableNames',names);

    MRSCont.QM.tables = addprop(MRSCont.QM.tables, {'VariableLongNames'}, {'variable'}); % add long name to table properties
    % Loop over field names to populate descriptive fields of table for JSON export
    for JJ = 1:length(names)
        switch names{JJ}
            case 'Cr_SNR'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'Cr_SNR'} = 'Signal to noise ratio of creatine';
                MRSCont.QM.tables.Properties.VariableDescriptions{'Cr_SNR'} = ['The maximum amplitude of the creatine peak divided by twice the standard deviation of the noise calculated from subspectrum ' name];
                MRSCont.QM.tables.Properties.VariableUnits{'Cr_SNR'} = 'arbitrary';
            case 'Cr_FWHM'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'Cr_FWHM'} = 'Full width at half maximum of creatine';
                MRSCont.QM.tables.Properties.VariableDescriptions{'Cr_FWHM'} = ['The width of the creatine peak at half the maximum amplitude calculated as the average of the FWHM of the data and the FWHM of a lorentzian fit calculated from subspectrum ' name];
                MRSCont.QM.tables.Properties.VariableUnits{'Cr_FWHM'} = 'Hz';
            case 'water_FWHM'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'water_FWHM'} = 'Full width at half maximum of reference water peak';
                MRSCont.QM.tables.Properties.VariableDescriptions{'water_FWHM'} = 'The width of the water peak at half the maximum amplitude calculated as the average of the FWHM of the data and the FWHM of a lorentzian fit';
                MRSCont.QM.tables.Properties.VariableUnits{'water_FWHM'} = 'Hz';
            case 'residual_water_ampl'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'residual_water_ampl'} = 'Residual water amplitude';
                MRSCont.QM.tables.Properties.VariableDescriptions{'residual_water_ampl'} = ['The amplitude of the water signal removed by the HSVD filter calculated from subspectrum ' name];
                MRSCont.QM.tables.Properties.VariableUnits{'residual_water_ampl'} = 'arbitrary';
            case 'freqShift'
                MRSCont.QM.tables.Properties.CustomProperties.VariableLongNames{'freqShift'} = 'Frequency shift';
                MRSCont.QM.tables.Properties.VariableDescriptions{'freqShift'} = ['Frequency shift calculated from the cross-correlation between spectrum ' name 'and the reference peaks (creatine and choline)'];
                MRSCont.QM.tables.Properties.VariableUnits{'freqShift'} = 'Hz';
        end
    end

end
end
