
function [MRSCont] = OspreyGenerateSpectra(nDatasets,outputFolder,saveSpec,parameter, alter,changedComb,zeroed,shareLorentzLB,SigRange,overwrite,add2HG)
if nargin < 11 
    add2HG=0;
    if nargin < 10
        overwrite.ph0 = [];
        overwrite.ph1 = [];
        overwrite.gaussLB = [];
        overwrite.lorentzLB = [];
        overwrite.freqShift = [];
        overwrite.metAmpl = [];
        overwrite.baseAmpl = [];
        overwrite.lineShape = [];
        overwrite.noiseAmpl = [];
        if nargin < 9
            SigRange = [1.9, 2.1];
            if nargin < 8
                    shareLorentzLB = 0;
                    if nargin < 7
                        zeroed.ph0 = 0;
                        zeroed.ph1 = 0;
                        zeroed.gaussLB = 0;
                        zeroed.lorentzLB = zeros(27,1);
                        zeroed.freqShift = zeros(27,1);
                        zeroed.metAmpl = zeros(35,1);
                        zeroed.baseAmpl = 0;
                        zeroed.lineShape = 0;
                        if nargin < 6
                            changedComb = 1;
                            if nargin <5
                                alter.Group1.ph0 = 0;
                                alter.Group1.ph1 = 0;
                                alter.Group1.gaussLB = 0;
                                alter.Group1.lorentzLB = zeros(27,1);
                                alter.Group1.freqShift = zeros(27,1);
                                alter.Group1.metAmpl = zeros(35,1);
                                alter.Group1.baseAmpl = zeros(14,1);
                                alter.Group1.lineShape = zeros(1,29);
                                alter.Group1.SNR = 0;        
                                alter.Group1.ph0_SD = 0;
                                alter.Group1.ph1_SD = 0;
                                alter.Group1.gaussLB_SD = 0;
                                alter.Group1.lorentzLB_SD = zeros(27,1);
                                alter.Group1.freqShift_SD = zeros(27,1);
                                alter.Group1.metAmpl_SD = zeros(35,1);
                                alter.Group1.baseAmpl_SD = zeros(14,1);
                                alter.Group1.lineShape_SD = zeros(1,29);
                                alter.Group1.SNR_SD = 0;
                                if nargin < 4
                                    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips.mat')));
                                    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
                                    if nargin <3
                                        saveSpec = 1;
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
                    end
            end
        end
    end
end
refSimTime = tic;
%% Load known values and dummies
prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
load(fullfile(prior_knowledge_folder,'DummyData.mat'));
load(fullfile(prior_knowledge_folder,'DummyContainer.mat'));
load(fullfile(prior_knowledge_folder,'DummyNII-MRSHeader.mat'))
MRSCont.flags.hasMMRef = 0;

%% Update dummy data
% Determine the ppm ranges of both the data and the basis functions.
dwelltime = 1/dataToFit.spectralwidth; % dwelltime [s]
npoints = 2048;
dataToFit.sz(1) = npoints;
% Calculate t and ppm arrays using the calculated parameters:
f = [(-dataToFit.spectralwidth/2)+(dataToFit.spectralwidth/(2*npoints)):dataToFit.spectralwidth/(npoints):(dataToFit.spectralwidth/2)-(dataToFit.spectralwidth/(2*npoints))];
ppmRangeData = f/(dataToFit.Bo*42.577);
% Philips data assumes the center frequency to be 4.68 ppm:
centerFreq = 4.68;
ppmRangeData=ppmRangeData + centerFreq;
dataToFit.ppm       = ppmRangeData;
ppmRangeData        = ppmRangeData';
% Calculate the time scale
dataToFit.t = (0:dataToFit.dwelltime:(dataToFit.sz(1)-1)*dataToFit.dwelltime);
dataToFit.fids = zeros(dataToFit.sz);
dataToFit.specs = zeros(dataToFit.sz);
%% Resample basis set
for bb = 1 : length(parameter.basisSet)
    load(parameter.basisSet{bb});
    BASIS = recalculateBasisSpecs(BASIS);
    BASIS = fit_sortBasisSet(BASIS);
    metabList = fit_createMetabList(parameter.metabolite_names(1:27)); 
    if add2HG
        metabList.bHG = 1;
    end
    % Create the modified basis set
    BASIS = fit_selectMetabs(BASIS, metabList, 1);
    [resBASIS{bb}] = fit_resampleBasis(dataToFit, BASIS);
end

%% Update nii header
if ~parameter.indirDim.flag
    nii_mrs.hdr.dim = [6, 1, 1, 1, 2048,1,1,1];
    nii_mrs.hdr.pixdim(5) = dataToFit.dwelltime;
    nii_mrs.hdr_ext.SpectrometerFrequency = dataToFit.Bo*42.577; %Set Bo
    nii_mrs.hdr_ext.dim_5  ='DIM_DYN'; %Set DIM_DYN
    nii_mrs.hdr_ext.EchoTime = resBASIS{1}.te/1000; 
    nii_mrs.hdr_ext.RepetitionTime = 2;
    nii_mrs.hdr_ext = rmfield(nii_mrs.hdr_ext, 'dim_6');
    nii_mrs.hdr_ext.Manufacturer = 'Philips';
    nii_mrs.hdr_ext.ManufacturersModelName = 'Simulated';
    nii_mrs.hdr_ext.DeviceSerialNumber = 'XXXXX';
    nii_mrs.hdr_ext.SoftwareVersions = 'XXXXX';
    nii_mrs.hdr_ext.InstitutionName = 'XXXXX';
    nii_mrs.hdr_ext.InstitutionAddress = 'XXXXX';
    nii_mrs.hdr_ext.RxCoil = 'Simulated';
    nii_mrs.hdr_ext.SequenceName = 'Simulated Philips SVS';
    nii_mrs.hdr_ext.ProtocolName = 'Simulation for MSM';
    nii_mrs.hdr_ext.PatientName = 'Simulated';
    nii_mrs.hdr_ext.PatientDoB = 'XXXXX';
    nii_mrs.hdr_ext.ConversionMethod = 'Osprey';
    nii_mrs.hdr_ext.PulseSequenceFile.Value = 'Simulated';
    nii_mrs.hdr_ext.PulseSequenceFile.Description = 'Simulated';
    nii_mrs.hdr_ext.IceProgramFile.Value = 'Simulated';
    nii_mrs.hdr_ext.IceProgramFile.Description = 'Simulated';
else
    switch parameter.indirDim.function
        case 'T2'
            nii_mrs.hdr.dim = [6, 1, 1, 1, 2048,parameter.indirDim.length,1,1];
            nii_mrs.hdr.pixdim(5) = dataToFit.dwelltime;
            nii_mrs.hdr_ext.SpectrometerFrequency = dataToFit.Bo*42.577; %Set Bo
            nii_mrs.hdr_ext.dim_5  ='DIM_INDIRECT_0'; %Set DIM_INDIRECT_0
            nii_mrs.hdr_ext.EchoTime = resBASIS{1}.te/1000; 
            nii_mrs.hdr_ext.RepetitionTime = 2;
            nii_mrs.hdr_ext = rmfield(nii_mrs.hdr_ext, 'dim_6');
            nii_mrs.hdr_ext.Manufacturer = 'Philips';
            nii_mrs.hdr_ext.ManufacturersModelName = 'Simulated';
            nii_mrs.hdr_ext.DeviceSerialNumber = 'XXXXX';
            nii_mrs.hdr_ext.SoftwareVersions = 'XXXXX';
            nii_mrs.hdr_ext.InstitutionName = 'XXXXX';
            nii_mrs.hdr_ext.InstitutionAddress = 'XXXXX';
            nii_mrs.hdr_ext.RxCoil = 'Simulated';
            nii_mrs.hdr_ext.SequenceName = 'Simulated Philips SVS';
            nii_mrs.hdr_ext.ProtocolName = 'Simulation for MSM';
            nii_mrs.hdr_ext.PatientName = 'Simulated';
            nii_mrs.hdr_ext.PatientDoB = 'XXXXX';
            nii_mrs.hdr_ext.ConversionMethod = 'Osprey';
            nii_mrs.hdr_ext.PulseSequenceFile.Value = 'Simulated';
            nii_mrs.hdr_ext.PulseSequenceFile.Description = 'Simulated';
            nii_mrs.hdr_ext.IceProgramFile.Value = 'Simulated';
            nii_mrs.hdr_ext.IceProgramFile.Description = 'Simulated';  
            nii_mrs.hdr_ext.dim_5_info = 'Echo time series scan';
            for bb = 1 : length(parameter.basisSet)
                nii_mrs.hdr_ext.dim_5_header.EchoTime(bb,1) = resBASIS{bb}.te/1000; 
            end
        case 'averages'
            nii_mrs.hdr.dim = [6, 1, 1, 1, 2048,parameter.indirDim.length,1,1];
            nii_mrs.hdr.pixdim(5) = dataToFit.dwelltime;
            nii_mrs.hdr_ext.SpectrometerFrequency = dataToFit.Bo*42.577; %Set Bo
            nii_mrs.hdr_ext.dim_5  ='DIM_DYN'; %Set DIM_DYN
            nii_mrs.hdr_ext.EchoTime = resBASIS{1}.te/1000; 
            nii_mrs.hdr_ext.RepetitionTime = 2;
            nii_mrs.hdr_ext = rmfield(nii_mrs.hdr_ext, 'dim_6');
            nii_mrs.hdr_ext.Manufacturer = 'Philips';
            nii_mrs.hdr_ext.ManufacturersModelName = 'Simulated';
            nii_mrs.hdr_ext.DeviceSerialNumber = 'XXXXX';
            nii_mrs.hdr_ext.SoftwareVersions = 'XXXXX';
            nii_mrs.hdr_ext.InstitutionName = 'XXXXX';
            nii_mrs.hdr_ext.InstitutionAddress = 'XXXXX';
            nii_mrs.hdr_ext.RxCoil = 'Simulated';
            nii_mrs.hdr_ext.SequenceName = 'Simulated Philips SVS';
            nii_mrs.hdr_ext.ProtocolName = 'Simulation for MSM';
            nii_mrs.hdr_ext.PatientName = 'Simulated';
            nii_mrs.hdr_ext.PatientDoB = 'XXXXX';
            nii_mrs.hdr_ext.ConversionMethod = 'Osprey';
            nii_mrs.hdr_ext.PulseSequenceFile.Value = 'Simulated';
            nii_mrs.hdr_ext.PulseSequenceFile.Description = 'Simulated';
            nii_mrs.hdr_ext.IceProgramFile.Value = 'Simulated';
            nii_mrs.hdr_ext.IceProgramFile.Description = 'Simulated';  
    end
end
    

%% Set up MRS Container
MRSCont.nDatasets = nDatasets;
if ~isempty(alter)
    GroupNames = fieldnames(alter);
    NoGroups = length(fieldnames(alter));
else
    NoGroups = 1;
end
% dataToFit = op_freqrange(dataToFit,0.5,4);
for kk = 1 : MRSCont.nDatasets *NoGroups
    MRSCont.files{kk} = sprintf([outputFolder filesep 'simulated-raw' filesep 'sim-' '%03d' filesep 'ses-001' filesep 'sim-' '%03d' '_ses-001_PRESS_' '%1d' 'T_' '%d' '_TE'],kk,kk,round(dataToFit.Bo),round(dataToFit.te) );
end
MRSCont.flags.simulated = 1;
MRSCont.flags.isSERIES = 0;
MRSCont.flags.isSPECIAL = 0;
MRSCont.opts.MultipleSpectra.metab = [1];
%% Generate Distributions for all parameters accoring to each value from Big PRESS
params = {'ph0','ph1','gaussLB','lorentzLB','freqShift','metAmpl','baseAmpl','lineShape','SNR'};

for p = 1 : length(params)
    par.(params{p}) = [];
end

for gg = 1 : NoGroups
    for p = 1 : length(params)

        if p == 9
            p;
        end
        eval([ 'means = parameter.' params{p} '.mean +( alter.' (GroupNames{gg}) '.' params{p} ' ./100.*parameter.' params{p} '.mean);' ]);
        
        eval([ 'rel =' 'alter.' (GroupNames{gg}) '.rel.' params{p} '_SD;' ]);

        if rel
            eval([ 'SD = ( abs(means) .* alter.' (GroupNames{gg}) '.' params{p} '_SD/100);' ]);
        else
            eval([ 'SD =parameter.' params{p} '.SD +( alter.' (GroupNames{gg}) '.' params{p} '_SD ./100.*parameter.' params{p} '.SD);' ]);
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


par.gaussLB = sqrt(par.gaussLB);

if shareLorentzLB
    par.lorentzLB = repmat(par.lorentzLB(:,1), [1 size(par.lorentzLB,2)]);
end

% Let's replace the anti-correlated amplitudes to create more realisitc
% results.
if changedComb
    par.ampl(:,13) = (par.metAmpl(:,28) + par.metAmpl(:,29))/2; %NAA+NAAG
    par.metAmpl(:,14) = (par.metAmpl(:,28) - par.metAmpl(:,29))/2; %NAA+NAAG
    par.metAmpl(:,3) = (par.metAmpl(:,30) + par.metAmpl(:,31))/2; %Cr+PCr
    par.metAmpl(:,16) = (par.metAmpl(:,30) - par.metAmpl(:,31))/2; %Cr+PCr
    par.metAmpl(:,15) = (par.metAmpl(:,32) + par.metAmpl(:,33))/2; %PCh+GPC
    par.metAmpl(:,6) = (par.metAmpl(:,32) - par.metAmpl(:,33))/2; %PCh+GPC
    par.metAmpl(:,9) = (par.metAmpl(:,34) + par.metAmpl(:,35))/2; %Glu+Gln
    par.metAmpl(:,8) = (par.metAmpl(:,34) - par.metAmpl(:,35))/2; %Glu+Gln
else
    for kk = 1 : MRSCont.nDatasets * NoGroups
        if (par.metAmpl(kk,28) - par.metAmpl(kk,14)) > 0
            par.metAmpl(kk,13) = par.metAmpl(kk,28) - par.metAmpl(kk,14); %NAA+NAAG
        else
            par.metAmpl(kk,14) = par.metAmpl(kk,28) - par.metAmpl(kk,13); %NAA+NAAG
        end
        if (par.metAmpl(kk,30) - par.metAmpl(kk,16)) > 0
            par.metAmpl(kk,3) = par.metAmpl(kk,30) -  par.metAmpl(kk,16); %Cr+PCr
        else
            par.metAmpl(kk,16) = par.metAmpl(kk,30) -  par.metAmpl(kk,3); %Cr+PCr
        end
        if (par.metAmpl(kk,32) - par.metAmpl(kk,6)) > 0
            par.metAmpl(kk,15) = par.metAmpl(kk,32) - par.metAmpl(kk,6); %PCh+GPC
        else
            par.metAmpl(kk,6) = par.metAmpl(kk,32) - par.metAmpl(kk,15); %PCh+GPC
        end
        if (par.metAmpl(kk,34) - par.metAmpl(kk,8)) > 0
            par.metAmpl(kk,9) = par.metAmpl(kk,34) - par.metAmpl(kk,8); %Glu+Gln
        else
            par.metAmpl(kk,8) = par.metAmpl(kk,34) - par.metAmpl(kk,9); %Glu+Gln
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
if sum(zeroed.metAmpl) > 0
    for m = 1 : 35
        if zeroed.metAmpl(m)
            par.metAmpl(:,m) = zeros(MRSCont.nDatasets,1);
        end
    end
end
if zeroed.baseAmpl
    par.baseAmpl =  zeros(MRSCont.nDatasets,size(par.baseAmpl(1,:),2));
end
if zeroed.lineShape
    par.lineShape =  zeros(MRSCont.nDatasets,size(par.lineShape(1,:),2));
end


%% Overwrite parameters
if length(parameter.metabolite_names) > 1
    nBasis = 27;
    nCombs = 35;
else
    nBasis = 1;
    nCombs = 1;
end
par.metAmpl = par.metAmpl(:,1:nBasis);
if ~isempty(overwrite.ph0)
    par.ph0 =  overwrite.ph0;
end
if ~isempty(overwrite.ph1)
    par.ph1 =  overwrite.ph1;
end
if ~isempty(overwrite.gaussLB)
    par.gaussLB =  overwrite.gaussLB;
end
if ~isempty(overwrite.lorentzLB)
   par.lorentzLB = overwrite.lorentzLB;
end
if ~isempty(overwrite.freqShift)
   par.freqShift = overwrite.freqShift;
end
if ~isempty(overwrite.metAmpl)
   par.metAmpl = overwrite.metAmpl;
end
if ~isempty(overwrite.baseAmpl)
   par.baseAmpl = overwrite.baseAmpl;
end
if ~isempty(overwrite.lineShape)
   par.lineShape = overwrite.lineShape;
end
if ~isempty(overwrite.noiseAmpl)
   noiseAmpl = overwrite.noiseAmpl;
else
   noiseAmpl = [];
end

%% Add 2-HG
if add2HG
    par.lorentzLB(:,end+1) = par.lorentzLB(:,end);
    par.freqShift(:,end+1) = par.freqShift(:,end);
    par.metAmpl(:,end+1) = (par.metAmpl(:,3) + par.metAmpl(:,16)) * 0.5;
    bkpPar = par;
    par.lorentzLB(:,3) = par.lorentzLB(:,end);
    par.freqShift(:,3) = par.freqShift(:,end);
    par.metAmpl(:,3) = par.metAmpl(:,end);
    par.lorentzLB(:,4:end) = bkpPar.lorentzLB(:,3:end-1);
    par.freqShift(:,4:end) = bkpPar.freqShift(:,3:end-1);
    par.metAmpl(:,4:end) = bkpPar.metAmpl(:,3:end-1);
end

%% Apply indirect dimensions
% parameter.indirDim.parameter = 'metAmpl';
if parameter.indirDim.flag
    switch parameter.indirDim.function
        case 'T2'
            for bb = 1 : length(parameter.basisSet)
                TEs(bb) = resBASIS{bb}.te;
            end
            expectation = zeros(MRSCont.nDatasets,length(parameter.indirDim.expectation.mean));
            for ind = 1 : length(parameter.indirDim.expectation.mean)
                expectation(:,ind) = normrnd(parameter.indirDim.expectation.mean(ind),parameter.indirDim.expectation.SD(ind),MRSCont.nDatasets,1);
            end
            ampl_factors = exp(-permute(repmat(repmat(TEs,[length(parameter.indirDim.expectation.mean) 1]),[1 1 MRSCont.nDatasets]),[3 1 2]) ./ repmat(expectation, [1 1 length(parameter.basisSet)]));
            ampl_factors_norm = ampl_factors ./ ampl_factors(:,:,1);
            par.indirDim.expectation = expectation;
        case 'averages'
            ampl_factors_norm = ones([size(par.metAmpl) parameter.indirDim.length]);

    end
else
    ampl_factors_norm = ones(size(par.metAmpl));
end

if isfield(par,'baseAmpl')
    if parameter.indirDim.flag
        switch parameter.indirDim.function
            case 'T2'
                for bb = 1 : length(parameter.basisSet)
                    TEs(bb) = resBASIS{bb}.te;
                end
                expectation = ones(MRSCont.nDatasets,length(par.baseAmpl)) .* normrnd(22,1.5,MRSCont.nDatasets,1);
                baseAmpl_factors = exp(-permute(repmat(repmat(TEs,[length(par.baseAmpl) 1]),[1 1 MRSCont.nDatasets]),[3 1 2]) ./ repmat(expectation, [1 1 length(parameter.basisSet)]));
                baseAmpl_factors_norm = baseAmpl_factors ./ baseAmpl_factors(:,:,1);
                par.indirDim.expectationBase = expectation;
        end
    else
        baseAmpl_factors_norm = ones(size(par.baseAmpl));
    end
end


%Save the results
MRSCont.in_silico.par_full = par;
MRSCont.in_silico.par_full.metabolite_names = parameter.metabolite_names(1:nBasis);

%% Generate in vivo like spectrum
MRSCont.nDatasets = nDatasets *NoGroups;
basisSetBckp = resBASIS;
for kk = 1 : MRSCont.nDatasets
    clc
    fprintf('Generate spectrum %i of %i. \n', kk,MRSCont.nDatasets);
    MRSCont.processed.metab{kk}     =dataToFit;
    MRSCont.processed.metab{kk}.fids = [];
    MRSCont.processed.metab{kk}.specs = [];
    for bb = 1 : length(parameter.basisSet)

        basisSet =basisSetBckp{bb};
        % ... fit parameters
        nMets       = basisSet.nMets;
        nMM         = basisSet.nMM;
        nBasisFcts  = nMets + nMM; % number of basis functions
        % lineShape   = fitParams.lineShape;
        ph0         = par.ph0(kk); % zero-order phase
        ph1         = par.ph1(kk); % first-order phase 
        gaussLB     = par.gaussLB(kk); % Gaussian damping [Hz]
        lorentzLB   = par.lorentzLB(kk,:); % Lorentzian damping [Hz] for each basis function
        freqShift   = par.freqShift(kk,:); % Frequency shift [Hz] for each basis function
        ampl        = par.metAmpl(kk,:); % Amplitudes for metabolite/MM/lipid basis functions
        if isfield(par,'baseAmpl') && sum(par.baseAmpl(kk,:)) > 0
            beta_j      = par.baseAmpl(kk,:); % Amplitudes for baseline spline basis functions
            if length(parameter.basisSet) > 1
                beta_j = beta_j .* squeeze(baseAmpl_factors_norm(kk,:,bb));
            else
                beta_j = beta_j .* baseAmpl_factors_norm(kk,:,:);
            end
        end
        refShift    = 0; % Reference shift applied to the data during first step of fitting
        
    
        % Apply indirect factor
        if sum(size(ampl_factors_norm)) ~= sum(size(par.metAmpl))
            if length(parameter.basisSet) > 1
                ampl = ampl .* squeeze(ampl_factors_norm(kk,:,bb));
            else
                ampl = ampl .* ampl_factors_norm(kk,:,:);
            end
        else
            ampl = ampl .* ampl_factors_norm(kk,:);
        end
        
        % Normalize
        if isfield(par,'lineShape')
            if sum(par.lineShape(kk,:)) > 0 
                lineShape = par.lineShape(kk,:)/sum(par.lineShape(kk,:));
            end
        end
    
    
        %%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
        % Run the time-domain operations on the metabolite basis functions
        % (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
        ppm_ax = basisSet.ppm;
        t = basisSet.t;
        for ii=1:nBasisFcts
            basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp(-1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB*gaussLB.*t.*t)' .* exp(1i.*ph1.*ppm_ax);    
            basisSet.fids(:,ii) = basisSet.fids(:,ii) * exp(1i*ph0);
        end
        basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);
    
        % Run the frequency-domain operations on the basis functions
        % (first order phase correction)
        % Cut out the frequency range of the basis set
        % Create a ppm vector around a pivot point (water)
        % ppm_ax = basisSet.ppm;
        % pivotPoint = 4.68;
        % multiplier = ppm_ax - pivotPoint;
    
        % Apply the linear phase correction
        % for ii=1:nBasisFcts
        %     basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp(1i.*ph1.*multiplier);
        % end
        % basisSet.fids = ifft(ifftshift(basisSet.specs,1),[],1);
    
        % Apply phasing to the spline basis functions
        if isfield(par,'baseAmpl')  && sum(par.baseAmpl(kk,:)) > 0
            [splineArrayToAdd] = fit_makeSplineBasis(dataToFit, [0.5 4], 0.4, 0);
            Bsp = [splineArrayToAdd(:,:,1) + 1i*splineArrayToAdd(:,:,2)];
            B = zeros(2048,12);
            B(ppm_ax>0.5 & ppm_ax < 4,:) = Bsp;
            B = B  * exp(1i*ph0);
            B = B .* exp(1i*ph1*ppm_ax);
        end
    
    
        %%% 3. APPLY THE LINEAR PARAMETERS %%%
        % Convolve the lineshape with the metabolite basis functions only
        % (NOT the macromolecules or lipids or baseline splines).
        A = basisSet.specs;
        if isfield(par,'lineShape')
            
            if sum(par.lineShape(kk,:)) > 0
                A = real(basisSet.specs);
                for nb = 1:basisSet.nMets
                    A(:,nb) = conv(A(:,nb), lineShape, 'same');
                end
                A = hilbert(A);
            end
            
        end

        % Calculate the final baseline ans sum up everything
        if isfield(par,'baseAmpl')  && sum(par.baseAmpl(kk,:)) > 0
            baseline    = B * beta_j(2:end-1)';
        else
            baseline = 0;
        end
        if size(A,3) == size(ampl,3)
            spectrum = A * ampl' + baseline;
        else
           spectrum = repmat(A , [1 1 parameter.indirDim.length]) .* repmat(ampl, [2048 1 1]) + baseline;
        end
        if ndims(spectrum) == 3
            spectrum = squeeze(sum(spectrum,2));
        end
       
       

        % Create the output resampled basis set container
        MRSCont.processed.metab{kk}.ppm     = ppmRangeData';
        if length(parameter.basisSet) > 1
            MRSCont.processed.metab{kk}.specs(:,bb)   = (spectrum);
            MRSCont.processed.metab{kk}.fids(:,bb)    = ifft(ifftshift((spectrum),dataToFit.dims.t),[],dataToFit.dims.t);
        else
            MRSCont.processed.metab{kk}.specs   = (spectrum);
            MRSCont.processed.metab{kk}.fids    = ifft(ifftshift((spectrum),dataToFit.dims.t),[],dataToFit.dims.t);
        end
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
          
        % figure, plot(MRSCont.processed.metab{kk}.t,real(((MRSCont.processed.metab{kk}.fids(:,1)))))
        % figure, plot(MRSCont.processed.metab{kk}.ppm,real(fftshift(fft(MRSCont.processed.metab{kk}.fids(:,1),[],1),1)))

                     
    end
    if parameter.indirDim.flag
        switch parameter.indirDim.function
            case 'T2'
                seq = MRSCont.processed.metab{kk}.seq;
                MRSCont.processed.metab{kk} = rmfield(MRSCont.processed.metab{kk},'seq');
                MRSCont.processed.metab{kk}.seq = repmat({seq},[1 length(parameter.basisSet)]);       
                MRSCont.processed.metab{kk}.spectralwidth = repmat(MRSCont.processed.metab{kk}.spectralwidth,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.dwelltime = repmat(MRSCont.processed.metab{kk}.dwelltime,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.txfrq = repmat(MRSCont.processed.metab{kk}.txfrq,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.te = TEs;
                MRSCont.processed.metab{kk}.tr = repmat(MRSCont.processed.metab{kk}.tr,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.centerFreq = repmat(MRSCont.processed.metab{kk}.centerFreq,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.pointsToLeftshift = repmat(MRSCont.processed.metab{kk}.pointsToLeftshift,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.specReg = repmat(MRSCont.processed.metab{kk}.specReg,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.exp_var = TEs;
                MRSCont.processed.metab{kk}.extra_names = repmat({'Series'},[1 length(parameter.basisSet)]); 
                MRSCont.processed.metab{kk}.names = repmat(MRSCont.processed.metab{kk}.names,[1 length(parameter.basisSet)])';
                MRSCont.processed.metab{kk}.watersupp = repmat(MRSCont.processed.metab{kk}.watersupp,[1 length(parameter.basisSet)])';
                MRSCont.processed.metab{kk}.extras = length(TEs);
                MRSCont.processed.metab{kk}.dims.extras = 2;
            case 'averages'
                seq = MRSCont.processed.metab{kk}.seq;
                MRSCont.processed.metab{kk} = rmfield(MRSCont.processed.metab{kk},'seq');
                MRSCont.processed.metab{kk}.seq = repmat({seq},[1 length(parameter.basisSet)]);       
                MRSCont.processed.metab{kk}.spectralwidth = repmat(MRSCont.processed.metab{kk}.spectralwidth,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.dwelltime = repmat(MRSCont.processed.metab{kk}.dwelltime,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.txfrq = repmat(MRSCont.processed.metab{kk}.txfrq,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.te = resBASIS{bb}.te;
                MRSCont.processed.metab{kk}.tr = repmat(MRSCont.processed.metab{kk}.tr,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.centerFreq = repmat(MRSCont.processed.metab{kk}.centerFreq,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.pointsToLeftshift = repmat(MRSCont.processed.metab{kk}.pointsToLeftshift,[1 length(parameter.basisSet)]);
                MRSCont.processed.metab{kk}.specReg = repmat(MRSCont.processed.metab{kk}.specReg,[1 parameter.indirDim.length]);
                MRSCont.processed.metab{kk}.names = repmat(MRSCont.processed.metab{kk}.names,[1 length(parameter.basisSet)])';
                MRSCont.processed.metab{kk}.watersupp = repmat(MRSCont.processed.metab{kk}.watersupp,[1 length(parameter.basisSet)])';
                MRSCont.processed.metab{kk}.averages = parameter.indirDim.length;
                MRSCont.processed.metab{kk}.rawAverages = parameter.indirDim.length;
        end
    end

    % Calculate noise to match the wanted SNR value
    if isempty(noiseAmpl)
            if strcmp(parameter.indirDim.function,'averages')
                   MRSCont.processed.metab{kk}.dims.averages=2;
                    MRSCont.processed.metab{kk}.flags.averaged=0;
            end
            peak_spec = op_freqrange(MRSCont.processed.metab{kk},SigRange(1),SigRange(2));
            peak_spec = op_averaging(peak_spec);
            peak_max = max(real(peak_spec.specs(:,1)));
            noise_ampl = peak_max/par.SNR(kk);
            noise_amplitude = noise_ampl * 1/sqrt(2048) * sqrt(MRSCont.processed.metab{kk}.sz(2)); % Remove scaling from ifft in MATLAB
            [MRSCont.processed.metab{kk},noise_vector]=op_addNoise(MRSCont.processed.metab{kk},noise_amplitude);
            

            if strcmp(parameter.indirDim.function,'averages')
                       MRSCont.processed.metab{kk}.dims.averages=0;
                        MRSCont.processed.metab{kk}.flags.averaged=1;
            end
        MRSCont.in_silico.noise_amplitude(kk) = noise_amplitude;
    else
        [MRSCont.processed.metab{kk},noise_vector]=op_addNoise(MRSCont.processed.metab{kk},noiseAmpl(kk));
        MRSCont.in_silico.noise_amplitude(kk) = noiseAmpl(kk);
        if strcmp(parameter.indirDim.function,'averages')
           MRSCont.processed.metab{kk}.dims.averages=2;
           MRSCont.processed.metab{kk}.flags.averaged=0;
        end
    end
    
    MRSCont.in_silico.noise{kk} = noise_vector;
    

    MRSCont.processed.metab{kk}.flags.isUnEdited = 1;
    MRSCont.processed.metab{kk}.flags.isMEGA = 0;
    MRSCont.processed.metab{kk}.flags.isHERMES = 0;
    MRSCont.processed.metab{kk}.flags.isHERCULES = 0;
    MRSCont.raw{kk} = MRSCont.processed.metab{kk};
           
end
%% Sort extras 
MRSCont.nDatasets = [MRSCont.nDatasets , length(parameter.basisSet)];
if parameter.indirDim.flag && ~strcmp(parameter.indirDim.function,'averages')
    
    rawBckp = MRSCont.raw;
    MRSCont.raw = [];
    for kk = 1 : MRSCont.nDatasets(1)
        for bb = 1 :  length(parameter.basisSet)
            raw = rawBckp{kk};
            raw=op_takeextra(raw,bb);
            MRSCont.raw{bb,kk} = raw;
        end
    end
    MRSCont.opts.MultipleSpectra.metab = [1:MRSCont.nDatasets(1)];
    MRSCont.opts.MultipleSpectra.mm = [];
    MRSCont.opts.MultipleSpectra.mm_ref = [];
    MRSCont.opts.MultipleSpectra.ref = [];
    MRSCont.opts.MultipleSpectra.w = [];
end

%% Calculate QM & add Info
 SubSpec = {'A'};
% Calculate some spectral quality metrics here;

metab_ll = 1;
mm_ll = 1;
ref_ll = 1;
w_ll = 1;
ref_mm_ll=1;
for kk = 1 : MRSCont.nDatasets(1)
    clc
    fprintf('QM spectrum %i of %i. \n', kk,MRSCont.nDatasets(1));
    for bb = 1 : MRSCont.nDatasets(2)
        if length(parameter.metabolite_names) > 1 && ~strcmp(parameter.indirDim.function,'averages')
            raw = MRSCont.processed.metab{metab_ll,kk};
            if raw.dims.extras > 1
                raw=op_takeextra(raw,bb);
            end
            try
                [raw,SNR] = op_get_Multispectra_SNR(raw);
                
            catch

            end
            try
                FWHM = op_get_Multispectra_LW(raw);
            catch
                FWHM{1} = nan;
            end
            MRSCont.processed.metab{metab_ll,kk}.QC_names = raw.QC_names;
            if exist('raw_no_subspec_aling','var')
                MRSCont.processed_no_align.metab{metab_ll,kk}     = raw_no_subspec_aling;
            end
            for ss = 1 : length(SubSpec)
                MRSCont.QM.SNR.metab(bb,kk,ss)    = SNR{ss};
                MRSCont.QM.FWHM.metab(bb,kk,ss)   = FWHM(ss); % in Hz
                MRSCont.QM.freqShift.metab(bb,kk,ss)  = refShift;
                MRSCont.QM.res_water_amp.metab(bb,kk,ss) = 0;
                if strcmp(SubSpec{ss},'diff1') ||strcmp(SubSpec{ss},'diff2') || strcmp(SubSpec{ss},'diff3') ||strcmp(SubSpec{ss},'sum')
                    if raw.flags.isMEGA
                        MRSCont.QM.drift.pre.(SubSpec{ss}){bb,kk}  = reshape([MRSCont.QM.drift.pre.A{kk}'; MRSCont.QM.drift.pre.B{kk}'], 1, [])';
                        MRSCont.QM.drift.post.(SubSpec{ss}){bb,kk} = reshape([MRSCont.QM.drift.post.A{kk}'; MRSCont.QM.drift.post.B{kk}'], 1, [])';
                    else
                        MRSCont.QM.drift.pre.(SubSpec{ss}){bb,kk}  = reshape([MRSCont.QM.drift.pre.A{kk}'; MRSCont.QM.drift.pre.B{kk}'; MRSCont.QM.drift.pre.C{kk}'; MRSCont.QM.drift.pre.D{kk}'], 1, [])';
                        MRSCont.QM.drift.post.(SubSpec{ss}){bb,kk} = reshape([MRSCont.QM.drift.post.A{kk}'; MRSCont.QM.drift.post.B{kk}'; MRSCont.QM.drift.post.C{kk}'; MRSCont.QM.drift.post.D{kk}'], 1, [])';
                    end
                end
                MRSCont.QM.drift.pre.AvgDeltaCr.(SubSpec{ss})(bb,kk) = 0;
                MRSCont.QM.drift.post.AvgDeltaCr.(SubSpec{ss})(bb,kk) =0;
            end
        else
            for ss = 1 : length(SubSpec)
                raw = MRSCont.processed.metab{metab_ll,kk};
                
                if strcmp(parameter.indirDim.function,'averages')
                   raw.dims.extras=0;
                   raw.dims.averages=2;
                   raw.flags.averaged=0;
                   raw = op_averaging(raw); 
                                      
                end
                if raw.dims.extras > 1
                    raw=op_takeextra(raw,bb);
                end
                MRSCont.QM.SNR.metab(bb,kk,ss)    = op_getSNR(raw,SigRange(1),SigRange(2));
                MRSCont.QM.FWHM.metab(bb,kk,ss)   = op_getLW(raw,SigRange(1),SigRange(2));
                MRSCont.QM.freqShift.metab(bb,kk,ss)  = refShift;
                MRSCont.QM.res_water_amp.metab(bb,kk,ss) = 0;
                MRSCont.QM.drift.pre.AvgDeltaCr.(SubSpec{ss})(bb,kk) = 0;
                MRSCont.QM.drift.post.AvgDeltaCr.(SubSpec{ss})(bb,kk) =0;
            end
            MRSCont.processed.metab{metab_ll,kk}.QC_names = {'water'};
            if exist('raw_no_subspec_aling','var')
                MRSCont.processed_no_align.metab{metab_ll,kk}     = raw_no_subspec_aling;
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


if nDatasets > 250
   MRSCont = rmfield(MRSCont, 'raw'); 
   MRSCont = rmfield(MRSCont, 'processed');
end

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