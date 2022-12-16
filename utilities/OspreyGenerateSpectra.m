function [MRSCont] = OspreyGenerateSpectra(nDatasets,outputFolder,saveSpec,alter,changedComb,zeroed)

if nargin < 6
    zeroed.ph0 = 0;
    zeroed.ph1 = 0;
    zeroed.gaussLB = 0;
    zeroed.lorentzLB = zeros(27,1);
    zeroed.freqShift = zeros(27,1);
    zeroed.ampl = zeros(35,1);
    zeroed.beta_j = 0;
    zeroed.lineShape = 0;
    if nargin < 5
        changedComb = 1;
        if nargin <4
            alter.Group1.ph0 = 0;
            alter.Group1.ph1 = 0;
            alter.Group1.gaussLB = 0;
            alter.Group1.lorentzLB = zeros(27,1);
            alter.Group1.freqShift = zeros(27,1);
            alter.Group1.ampl = zeros(35,1);
            alter.Group1.beta_j = zeros(14,1);
            alter.Group1.lineShape = zeros(1,29);
            alter.Group1.SNR = 0;
    
            alter.Group1.ph0_SD = 0;
            alter.Group1.ph1_SD = 0;
            alter.Group1.gaussLB_SD = 0;
            alter.Group1.lorentzLB_SD = zeros(27,1);
            alter.Group1.freqShift_SD = zeros(27,1);
            alter.Group1.ampl_SD = zeros(35,1);
            alter.Group1.beta_j_SD = zeros(14,1);
            alter.Group1.lineShape_SD = zeros(1,29);
            alter.Group1.SNR_SD = 0;
            if nargin <3
                saveSpec = 0;
                if nargin<2
                     [settingsFolder,~,~] = fileparts(which('OspreySettings.m'));
                    allFolders      = strsplit(settingsFolder, filesep);
                    outputFolder       = [strjoin(allFolders(1:end-1), filesep) 'simulated' filesep]; % parent folder (= Osprey folder)
                    if nargin<1
                        nDatasets=1;
                    end
                end
            end
        end
    end
end
refSimTime = tic;
%% Load known values and dummies
prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips.mat')));
load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Lookup_with_updated_lineShape.mat'));
load(fullfile(prior_knowledge_folder,'splineArray.mat'));
load(fullfile(prior_knowledge_folder,'ResampledBASIS.mat'));
load(fullfile(prior_knowledge_folder,'DataDummy.mat'));
load(fullfile(prior_knowledge_folder,'DummyContainer.mat'));
load(fullfile(prior_knowledge_folder,'DummyNIIMRSHeader.mat'))
MRSCont.flags.hasMMRef = 0;

%% Update nii header
nii_mrs.hdr.dim = [6, 1, 1, 1, 2048,1,1,1];
nii_mrs.hdr.pixdim(5) = dataToFit.dwelltime;
nii_mrs.hdr_ext.SpectrometerFrequency = dataToFit.Bo*42.577; %Set Bo
nii_mrs.hdr_ext.dim_5  ='DIM_DYN'; %Set DIM_DYN
nii_mrs.hdr_ext.EchoTime = 0.035; 
nii_mrs.hdr_ext.RepetitionTime = 2;
nii_mrs.hdr_ext = rmfield(nii_mrs.hdr_ext, 'dim_6');
nii_mrs.hdr_ext.Manufacturer = 'Philips';
nii_mrs.hdr_ext.ManufacturersModelName = 'Simulated';
nii_mrs.hdr_ext.DeviceSerialNumber = 'XXXXX';
nii_mrs.hdr_ext.SoftwareVersions = 'XXXXX';
nii_mrs.hdr_ext.InstitutionName = 'XXXXX';
nii_mrs.hdr_ext.InstitutionAddress = 'XXXXX';
nii_mrs.hdr_ext.RxCoil = 'Simulated';
nii_mrs.hdr_ext.SequenceName = 'Simulated Philips PRESS SVS';
nii_mrs.hdr_ext.ProtocolName = 'Simulation for MSM';
nii_mrs.hdr_ext.PatientName = 'Simulated';
nii_mrs.hdr_ext.PatientDoB = 'XXXXX';
nii_mrs.hdr_ext.ConversionMethod = 'Osprey';
nii_mrs.hdr_ext.PulseSequenceFile.Value = 'Simulated';
nii_mrs.hdr_ext.PulseSequenceFile.Description = 'Simulated';
nii_mrs.hdr_ext.IceProgramFile.Value = 'Simulated';
nii_mrs.hdr_ext.IceProgramFile.Description = 'Simulated';
%%
% Determine the ppm ranges of both the data and the basis functions.
    dwelltime = 1/dataToFit.spectralwidth; % dwelltime [s]
    npoints = 2048;
    % Calculate t and ppm arrays using the calculated parameters:
    f = [(-dataToFit.spectralwidth/2)+(dataToFit.spectralwidth/(2*npoints)):dataToFit.spectralwidth/(npoints):(dataToFit.spectralwidth/2)-(dataToFit.spectralwidth/(2*npoints))];
    ppmRangeData = f/(dataToFit.Bo*42.577);
    % Philips data assumes the center frequency to be 4.68 ppm:
    centerFreq = 4.68;
    ppmRangeData=ppmRangeData + centerFreq;

    ppmRangeData        = ppmRangeData';
    
%%
basisSet = BASIS;
% To do: Interface with interactive user input
basisSet = fit_sortBasisSet(basisSet);
metabList = fit_createMetabList(ampl_names(1:27));

% Create the modified basis set
basisSet = fit_selectMetabs(basisSet, metabList, 1);
%% Set up MRS Container
MRSCont.nDatasets = nDatasets;
if ~isempty(alter)
    GroupNames = fieldnames(alter);
    NoGroups = length(fieldnames(alter));
else
    NoGroups = 1;
end
dataToFit = op_freqrange(dataToFit,0.5,4);
for kk = 1 : MRSCont.nDatasets *NoGroups
    MRSCont.files{kk} = sprintf([outputFolder filesep 'simulated-raw' filesep 'sim-' '%03d' filesep 'ses-001' filesep 'sim-' '%03d' '_ses-001_PRESS_' '%1d' 'T_' '%d' '_TE'],kk,kk,round(dataToFit.Bo),round(dataToFit.te) );
end
MRSCont.flags.simulated = 1;
MRSCont.flags.isSERIES = 0;
MRSCont.flags.isSPECIAL = 0;
MRSCont.opts.MultipleSpectra.metab = [1];
%% Generate Distributions for all parameters accoring to each value from Big PRESS
params = {'ph0','ph1','gaussLB','lorentzLB','freqShift','ampl','beta_j','lineShape','SNR'};

for p = 1 : length(params)
    par.(params{p}) = [];
end

for gg = 1 : NoGroups
    for p = 1 : length(params)
        eval([ 'means =' params{p} '+( alter.' (GroupNames{gg}) '.' params{p} ' ./100.*' params{p} ');' ]);
        
        eval([ 'rel =' 'alter.' (GroupNames{gg}) '.rel.' params{p} '_SD;' ]);

        if rel
            eval([ 'SD = ( abs(means) .* alter.' (GroupNames{gg}) '.' params{p} '_SD/100);' ]);
        else
            eval([ 'SD =' params{p} '_SD +( alter.' (GroupNames{gg}) '.' params{p} '_SD ./100.*' params{p} '_SD);' ]);
        end
        
        if length(means) == 1
                if p ~= 3
                    par.(params{p}) = normrnd(means,SD,MRSCont.nDatasets,1);
                else
                    par.(params{p}) = abs(normrnd(means,SD,MRSCont.nDatasets,1));
                end
        else
            for pp = 1 : length(means)
                if p ~=  4 && p ~=  6 && p ~=  7 && p ~=  8 && p ~=  9
                    par.(params{p})(1 + nDatasets*(gg-1):nDatasets*gg,pp) = normrnd(means(pp),SD(pp),MRSCont.nDatasets,1);
                else
                    par.(params{p})(1 + nDatasets*(gg-1):nDatasets*gg,pp) = abs(normrnd(means(pp),SD(pp),MRSCont.nDatasets,1));
                end
            end
        end
    end      
end


% Let's replace the anti-correlated amplitudes to create more realisitc
% results.
if changedComb
    par.ampl(:,13) = (par.ampl(:,28) + par.ampl(:,29))/2; %NAA+NAAG
    par.ampl(:,14) = (par.ampl(:,28) - par.ampl(:,29))/2; %NAA+NAAG
    par.ampl(:,3) = (par.ampl(:,30) + par.ampl(:,31))/2; %Cr+PCr
    par.ampl(:,16) = (par.ampl(:,30) - par.ampl(:,31))/2; %Cr+PCr
    par.ampl(:,15) = (par.ampl(:,32) + par.ampl(:,33))/2; %PCh+GPC
    par.ampl(:,6) = (par.ampl(:,32) - par.ampl(:,33))/2; %PCh+GPC
    par.ampl(:,9) = (par.ampl(:,34) + par.ampl(:,35))/2; %Glu+Gln
    par.ampl(:,8) = (par.ampl(:,34) - par.ampl(:,35))/2; %Glu+Gln
else
    for kk = 1 : MRSCont.nDatasets * NoGroups
        if (par.ampl(kk,28) - par.ampl(kk,14)) > 0
            par.ampl(kk,13) = par.ampl(kk,28) - par.ampl(kk,14); %NAA+NAAG
        else
            par.ampl(kk,14) = par.ampl(kk,28) - par.ampl(kk,13); %NAA+NAAG
        end
        if (par.ampl(kk,30) - par.ampl(kk,16)) > 0
            par.ampl(kk,3) = par.ampl(kk,30) -  par.ampl(kk,16); %Cr+PCr
        else
            par.ampl(kk,16) = par.ampl(kk,30) -  par.ampl(kk,3); %Cr+PCr
        end
        if (par.ampl(kk,32) - par.ampl(kk,6)) > 0
            par.ampl(kk,15) = par.ampl(kk,32) - par.ampl(kk,6); %PCh+GPC
        else
            par.ampl(kk,6) = par.ampl(kk,32) - par.ampl(kk,15); %PCh+GPC
        end
        if (par.ampl(kk,34) - par.ampl(kk,8)) > 0
            par.ampl(kk,9) = par.ampl(kk,34) - par.ampl(kk,8); %Glu+Gln
        else
            par.ampl(kk,8) = par.ampl(kk,34) - par.ampl(kk,9); %Glu+Gln
        end
    end
end


par.SNR = abs(par.SNR);
%Save the results
MRSCont.in_silico.par_full = par;
%% Zeroing out
if zeroed.ph0
    par.ph0 =  zeros(MRSCont.nDatasets,1);
end
if zeroed.ph1
    par.ph1 =  zeros(MRSCont.nDatasets,1);
end
if zeroed.gaussLB
    par.gaussLB =  zeros(MRSCont.nDatasets,1);
end
if sum(zeroed.lorentzLB) > 0
    for m = 1 : 27
        if zeroed.lorentzLB(m)
            par.lorentzLB(:,m) = zeros(MRSCont.nDatasets,1);
        end
    end
end
if sum(zeroed.freqShift) > 0
    for m = 1 : 27
        if zeroed.freqShift(m)
            par.freqShift(:,m) = zeros(MRSCont.nDatasets,1);
        end
    end
end
if sum(zeroed.ampl) > 0
    for m = 1 : 35
        if zeroed.ampl(m)
            par.ampl(:,m) = zeros(MRSCont.nDatasets,1);
        end
    end
end
if zeroed.beta_j
    par.beta_j =  zeros(MRSCont.nDatasets,size(par.beta_j(1,:),2));
end
if zeroed.lineShape
    par.lineShape =  zeros(MRSCont.nDatasets,size(par.lineShape(1,:),2));
end
par.ampl = par.ampl(:,1:27);
%% Generate in vivo like spectrum
MRSCont.nDatasets = nDatasets *NoGroups;
basisSetBckp = basisSet;
for kk = 1 : MRSCont.nDatasets
    

    basisSet =basisSetBckp;
    % ... fit parameters
    nMets       = basisSet.nMets;
    nMM         = basisSet.nMM;
    nBasisFcts  = nMets + nMM; % number of basis functions
    % lineShape   = fitParams.lineShape;
    ph0         = par.ph0(kk) * pi/180; % zero-order phase correction [convert from deg to rad]
    ph1         = par.ph1(kk) * pi/180; % first-order phase correction [convert from deg/ppm to rad/ppm]
    gaussLB     = par.gaussLB(kk); % Gaussian damping [Hz]
    lorentzLB   = par.lorentzLB(kk,:); % Lorentzian damping [Hz] for each basis function
    freqShift   = par.freqShift(kk,:); % Frequency shift [Hz] for each basis function
    ampl        = par.ampl(kk,:); % Amplitudes for metabolite/MM/lipid basis functions
    beta_j      = par.beta_j(kk,:); % Amplitudes for baseline spline basis functions
    refShift    = 0; % Reference shift applied to the data during first step of fitting

    % Normalize
    if sum(par.lineShape(kk,:)) > 0 
        lineShape = par.lineShape(kk,:)/sum(par.lineShape(kk,:));
    end


    %%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
    % Run the time-domain operations on the metabolite basis functions
    % (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
    t = basisSet.t;
    for ii=1:nBasisFcts
        basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp(-1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)';    
        basisSet.fids(:,ii) = basisSet.fids(:,ii) * exp(1i*ph0);
    end
    basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

    % Run the frequency-domain operations on the basis functions
    % (first order phase correction)
    % Cut out the frequency range of the basis set
    basisSetCut = op_freqrange(basisSet,0.5,4);
    % Create a ppm vector around a pivot point (water)
    ppm_axCut = basisSetCut.ppm;
    ppm_ax = basisSet.ppm;
    pivotPoint = 4.68;
    multiplierCut = ppm_axCut - pivotPoint;
    multiplier = ppm_ax - pivotPoint;
    multiplier(multiplier < multiplierCut(1)) = 0;  
    multiplier(multiplier > multiplierCut(end)) = 0;

    % Apply the linear phase correction
    for ii=1:nBasisFcts
        basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp(1i.*ph1.*multiplier);
    end
    basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);
    
    for ii=1:nBasisFcts
        basisSet.fids(3072:end,ii) = 0;
    end
    basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

    % Apply phasing to the spline basis functions
    Bsp = [splineArrayToAdd(:,:,1) + 1i*splineArrayToAdd(:,:,2)];
    B = zeros(4096,14);
    B(884:884+1198-1,:) = Bsp;
    B = B  * exp(1i*ph0);
    B = B .* exp(1i*ph1*multiplier);
    B = [real(B)];


    %%% 3. APPLY THE LINEAR PARAMETERS %%%
    % Convolve the lineshape with the metabolite basis functions only
    % (NOT the macromolecules or lipids or baseline splines).
    A = real(basisSet.specs);
    if sum(par.lineShape(kk,:)) > 0 
        for nb = 1:basisSet.nMets
            A(:,nb) = conv(A(:,nb), lineShape, 'same');
        end
    end

    % Calculate the final baseline ans sum up everything
    baseline    = B * beta_j';
    spectrum = A * ampl' + baseline;
    specs(:) = hilbert(spectrum);
    dataToFit.fids = ifft(fftshift(specs',dataToFit.dims.t),[],dataToFit.dims.t);
    dataToFit.specs = specs';
%     figure(1), plot(real(specs)), hold on, plot(B* beta_j);
    dataToFit.sz  = size(dataToFit.fids);

    MRSCont.processed.metab{kk} = dataToFit;
    MRSCont.processed.metab{kk}.t = basisSet.t;
    MRSCont.processed.metab{kk}.ppm = basisSet.ppm';
        
    ppmRangeBasis       = MRSCont.processed.metab{kk}.ppm;
    ppmIsInDataRange    = (ppmRangeBasis < ppmRangeData(1)) & (ppmRangeBasis > ppmRangeData(end));
    if sum(ppmIsInDataRange) == 0
        ppmIsInDataRange    = (ppmRangeBasis > ppmRangeData(1)) & (ppmRangeBasis < ppmRangeData(end));
    end
    % Now resample the basis functions to match the resolution and frequency
    % range (ppm axis) of the data.
    fids_interp     = zeros(length(ppmRangeData), 1); % allocate memory
    specs_interp    = zeros(length(ppmRangeData), 1); % allocate memory

    specs_interp(:,1)      = interp1(ppmRangeBasis(ppmIsInDataRange), MRSCont.processed.metab{kk}.specs(ppmIsInDataRange,1), ppmRangeData, 'pchip', 'extrap');
    %convert back to time domain
    %if the length of Fids is odd, then you have to do a circshift of one to
    %make sure that you don't introduce a small frequency shift into the fids
    %vector.
    if mod(length(MRSCont.processed.metab{kk}.specs),2)==0
        %disp('Length of vector is even.  Doing normal conversion');
        fids_interp(:,1)   = ifft(fftshift(specs_interp(:,1),1),[],1);
    else
        %disp('Length of vector is odd.  Doing circshift by 1');
        fids_interp(:,1)   = ifft(circshift(fftshift(specs_interp(:,1),1)),[],1);
    end

%      fids_interp(:,1) = flipud(fids_interp(:,1));
%      fids_interp(end,1) = fids_interp(end-1,1);
%      specs_interp = fftshift(fft(fids_interp,[],1),1);
    % Create the output resampled basis set container
    MRSCont.processed.metab{kk}.ppm     = ppmRangeData';
    MRSCont.processed.metab{kk}.specs   = specs_interp;
    MRSCont.processed.metab{kk}.fids    = fids_interp;
    MRSCont.processed.metab{kk}.sz      = size(MRSCont.processed.metab{kk}.fids);
    MRSCont.processed.metab{kk}.n       = length(MRSCont.processed.metab{kk}.fids);
    MRSCont.processed.metab{kk}.dwelltime = dwelltime;
    % Calculate the new spectral width and dwelltime:
    dppm                        = abs(MRSCont.processed.metab{kk}.ppm(2)-MRSCont.processed.metab{kk}.ppm(1));
    ppmrange                    = abs((MRSCont.processed.metab{kk}.ppm(end)-MRSCont.processed.metab{kk}.ppm(1)))+dppm;
    MRSCont.processed.metab{kk}.spectralwidth   = ppmrange*MRSCont.processed.metab{kk}.Bo*42.577;
    MRSCont.processed.metab{kk}.dwelltime       = 1/MRSCont.processed.metab{kk}.spectralwidth;
    
    %Update nii header
    nii_mrs.hdr_ext.ConversionTime = datestr(now,30);
    MRSCont.processed.metab{kk}.nii_mrs = nii_mrs;

    % Calculate the time scale
    MRSCont.processed.metab{kk}.t = (0:MRSCont.processed.metab{kk}.dwelltime:(MRSCont.processed.metab{kk}.sz(1)-1)*MRSCont.processed.metab{kk}.dwelltime);
    
    
  % Calculate noise to match the wanted SNR value
    Cr_spec = op_freqrange(MRSCont.processed.metab{kk},2.9,3.1);
    Cr_max = max(real(Cr_spec.specs));
    noise_ampl = Cr_max/par.SNR(kk);    
    noise_real=noise_ampl*randn(MRSCont.processed.metab{kk}.sz);
    noise_imag=noise_ampl*randn(MRSCont.processed.metab{kk}.sz);
    noise=noise_real+(1i*noise_imag);
    MRSCont.in_silico.noise{kk} = noise;
    MRSCont.processed.metab{kk}.specs = MRSCont.processed.metab{kk}.specs + noise;
    MRSCont.processed.metab{kk}.fids = ifft(fftshift(MRSCont.processed.metab{kk}.specs,MRSCont.processed.metab{kk}.dims.t),[],MRSCont.processed.metab{kk}.dims.t);

    MRSCont.processed.metab{kk}.flags.isUnEdited = 1;
    MRSCont.processed.metab{kk}.flags.isMEGA = 0;
    MRSCont.processed.metab{kk}.flags.isHERMES = 0;
    MRSCont.processed.metab{kk}.flags.isHERCULES = 0;
    MRSCont.raw{kk} = MRSCont.processed.metab{kk};
           
end

%% Calculate QM & add Info
 SubSpec = {'A'};
% Calculate some spectral quality metrics here;

metab_ll = 1;
mm_ll = 1;
ref_ll = 1;
w_ll = 1;
ref_mm_ll=1;
for kk = 1 : MRSCont.nDatasets
    raw = MRSCont.processed.metab{metab_ll,kk};
    [raw,SNR] = op_get_Multispectra_SNR(raw);
    FWHM = op_get_Multispectra_LW(raw);
    MRSCont.processed.metab{metab_ll,kk} = raw;
    if exist('raw_no_subspec_aling','var')
        MRSCont.processed_no_align.metab{metab_ll,kk}     = raw_no_subspec_aling;
    end
    for ss = 1 : length(SubSpec)
        MRSCont.QM.SNR.metab(metab_ll,kk,ss)    = SNR{ss};
        MRSCont.QM.FWHM.metab(metab_ll,kk,ss)   = FWHM(ss); % in Hz
        MRSCont.QM.freqShift.metab(metab_ll,kk,ss)  = refShift;
        MRSCont.QM.res_water_amp.metab(metab_ll,kk,ss) = 0;
        if strcmp(SubSpec{ss},'diff1') ||strcmp(SubSpec{ss},'diff2') || strcmp(SubSpec{ss},'diff3') ||strcmp(SubSpec{ss},'sum')
            if raw.flags.isMEGA
                MRSCont.QM.drift.pre.(SubSpec{ss}){metab_ll,kk}  = reshape([MRSCont.QM.drift.pre.A{kk}'; MRSCont.QM.drift.pre.B{kk}'], 1, [])';
                MRSCont.QM.drift.post.(SubSpec{ss}){metab_ll,kk} = reshape([MRSCont.QM.drift.post.A{kk}'; MRSCont.QM.drift.post.B{kk}'], 1, [])';
            else
                MRSCont.QM.drift.pre.(SubSpec{ss}){metab_ll,kk}  = reshape([MRSCont.QM.drift.pre.A{kk}'; MRSCont.QM.drift.pre.B{kk}'; MRSCont.QM.drift.pre.C{kk}'; MRSCont.QM.drift.pre.D{kk}'], 1, [])';
                MRSCont.QM.drift.post.(SubSpec{ss}){metab_ll,kk} = reshape([MRSCont.QM.drift.post.A{kk}'; MRSCont.QM.drift.post.B{kk}'; MRSCont.QM.drift.post.C{kk}'; MRSCont.QM.drift.post.D{kk}'], 1, [])';
            end
        end
        MRSCont.QM.drift.pre.AvgDeltaCr.(SubSpec{ss})(metab_ll,kk) = 0;
        MRSCont.QM.drift.post.AvgDeltaCr.(SubSpec{ss})(metab_ll,kk) =0;
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


% Gather some more information from the processed data;
SubSpecNames = fieldnames(MRSCont.processed);
NoSubSpec = length(fieldnames(MRSCont.processed));
for ss = 1 : NoSubSpec
    for kk = 1 : MRSCont.nDatasets
            temp_sz(1,kk)= MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1);
            temp_sz_sw{1,kk} = ['np_sw_' num2str(MRSCont.processed.(SubSpecNames{ss}){1,kk}.sz(1)) '_' num2str(MRSCont.processed.(SubSpecNames{ss}){1,kk}.spectralwidth)];   
    end
    [MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint_spectralwidth,MRSCont.info.(SubSpecNames{ss}).unique_ndatapoint_spectralwidth_ind,~]  = unique(temp_sz_sw,'Stable');
    [MRSCont.info.(SubSpecNames{ss}).max_ndatapoint,MRSCont.info.(SubSpecNames{ss}).max_ndatapoint_ind] = max(temp_sz);
end

MRSCont.nDatasets(2) = 1;

time = toc(refSimTime);
MRSCont.runtime.Load = 0;
MRSCont.runtime.Proc = time;
%%
% Save the output structure to the output folder
% Determine output folder

outputFile      = MRSCont.outputFile;
if~isempty(outputFolder)
    MRSCont.outputFolder = outputFolder;
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end
end

for gg = 1 : NoGroups
    groups(1 + nDatasets*(gg-1):nDatasets*gg,1) = ones(nDatasets,1)*gg;  
end
statCSV = array2table(groups,'VariableNames',{'group'});
statCSV.subject = MRSCont.files';

writetable(statCSV,[MRSCont.outputFolder  filesep  'stat.csv']);
MRSCont.file_stat = [MRSCont.outputFolder  filesep  'stat.csv'];
MRSCont.flags.hasStatfile = 1;

MRSCont = calculate_QM(MRSCont);

if saveSpec
%     [MRSCont] = osp_saveLCM(MRSCont);
%     [MRSCont] = osp_saveJMRUI(MRSCont);
%     [MRSCont] = osp_saveVendor(MRSCont);
    [MRSCont] = osp_saveNII(MRSCont);
end
MRSCont.flags.didLoad = 1;
MRSCont.flags.hasMMRef = 0;

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end



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

function [spec,noise,step]=optimize_noise(spec,target)
    bckp = spec;  %0.00011    
    [spec,~]=op_addNoise(spec,0);
    SNR    = op_getSNR(spec,2.9,3.1); 
    step = 0;
    while (SNR > target)
        step = step + 1;
        [spec,~]=op_addNoise(bckp,step * 0.00001  );
        SNR    = op_getSNR(spec);       
    end
    [spec,noise]=op_addNoise(bckp,step * 0.00001);    
end
