function [MRSCont] = OspreyGenerateSpectra(nDatasets,outputFolder,saveSpec,alter,changedComb)
%% [MRSCont] = OspreyGenerateSpectra(MRSCont)
%   This function checks whether the structural image that the voxels were
%   coregistered to in OspreyCoreg has already been segmented by SPM12.
%
%   If it has not been, OspreySeg will call the SPM12 "New Segment"
%   function to perform segmentation into gray matter, white matter, and
%   CSF, and return fractional tissue volumes.
%
%   USAGE:
%       MRSCont = OspreySeg(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-08-21)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-08-21: First version of the code.
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
refSimTime = tic;
%% Load known values
load('/Users/helge/Documents/GitHub/osprey/ISMRM/MRS_BigPRESS_Philips_Lookup_with_updated_lineShape.mat');
load('/Users/helge/Documents/GitHub/osprey/ISMRM/splineArray.mat');
load('/Users/helge/Documents/GitHub/osprey/ISMRM/ResampledBASIS.mat');
load('/Users/helge/Documents/GitHub/osprey/ISMRM/DataDummy.mat');
load('/Users/helge/Documents/GitHub/osprey/ISMRM/DummyContainer.mat');
load('/Users/helge/Documents/GitHub/osprey/ISMRM/SNR_lookup.mat');

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
    MRSCont.files{kk} = [MRSCont.ospFolder filesep 'Simulated_Short_TE_brain_3T_' num2str(kk)];
end
MRSCont.flags.simulated = 1;
%% Generate Distributions for all parameters accoring to each value from Big PRESS
params = {'ph0','ph1','gaussLB','lorentzLB','freqShift','ampl','beta_j','lineShape','SNR'};

for p = 1 : length(params)
    par.(params{p}) = [];
end

for gg = 1 : NoGroups
    for p = 1 : length(params)
        eval([ 'means =' params{p} '+( alter.' (GroupNames{gg}) '.' params{p} ' ./100.*' params{p} ');' ]);
        eval([ 'SD =' params{p} '_SD +( alter.' (GroupNames{gg}) '.' params{p} '_SD ./100.*' params{p} '_SD);' ]);
        
        if length(means) == 1
            if p ==1 || p == 2
                par.(params{p}) = [par.(params{p}); means + randn(MRSCont.nDatasets,1)*SD/10];
            else
                par.(params{p}) = [par.(params{p}); means + randn(MRSCont.nDatasets,1)*SD];
            end
        else
            for pp = 1 : length(means)
                if p==7 ||p ==8
                    par.(params{p})(1 + nDatasets*(gg-1):nDatasets*gg,pp) = means(pp) + randn(MRSCont.nDatasets,1)*SD(pp)/2;
                else if p== 6
                        pd = makedist('Normal','mu',means(pp),'sigma',SD(pp)*0.9);
                        t = truncate(pd,0,inf);
                        par.(params{p})(1 + nDatasets*(gg-1):nDatasets*gg,pp) = random(t,MRSCont.nDatasets,1);
                    else
                        par.(params{p})(1 + nDatasets*(gg-1):nDatasets*gg,pp) = means(pp) + randn(MRSCont.nDatasets,1)*SD(pp);
                    end
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



%Save the results
MRSCont.in_silico.par_full = par;
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
    lineShape = par.lineShape(kk,:)/sum(par.lineShape(kk,:));


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
    Bsp = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
    B = zeros(4096,14);
    B(884:884+1198-1,:) = Bsp;
    B = B  * exp(1i*ph0);
    B = B .* exp(1i*ph1*multiplier);
    B = [real(B)];


    %%% 3. APPLY THE LINEAR PARAMETERS %%%
    % Convolve the lineshape with the metabolite basis functions only
    % (NOT the macromolecules or lipids or baseline splines).
    A = real(basisSet.specs);
    for nb = 1:basisSet.nMets
        A(:,nb) = conv(A(:,nb), lineShape, 'same');
    end

    % Calculate the final baseline ans sum up everything
    baseline    = B * beta_j';
    spectrum = A * ampl' + baseline;
    specs(:) = hilbert(spectrum);
    dataToFit.fids = ifft(fftshift(specs',dataToFit.dims.t),[],dataToFit.dims.t);
    dataToFit.specs = specs';
%     figure(1), plot(real(specs)), hold on, plot(B* beta_j);
    dataToFit.sz  = size(dataToFit.fids);

    MRSCont.processed.A{kk} = dataToFit;
    MRSCont.processed.A{kk}.t = basisSet.t;
    MRSCont.processed.A{kk}.ppm = basisSet.ppm';
        
    ppmRangeBasis       = MRSCont.processed.A{kk}.ppm;
    ppmIsInDataRange    = (ppmRangeBasis < ppmRangeData(1)) & (ppmRangeBasis > ppmRangeData(end));
    if sum(ppmIsInDataRange) == 0
        ppmIsInDataRange    = (ppmRangeBasis > ppmRangeData(1)) & (ppmRangeBasis < ppmRangeData(end));
    end
    % Now resample the basis functions to match the resolution and frequency
    % range (ppm axis) of the data.
    fids_interp     = zeros(length(ppmRangeData), 1); % allocate memory
    specs_interp    = zeros(length(ppmRangeData), 1); % allocate memory

    specs_interp(:,1)      = interp1(ppmRangeBasis(ppmIsInDataRange), MRSCont.processed.A{kk}.specs(ppmIsInDataRange,1), ppmRangeData, 'pchip', 'extrap');
    %convert back to time domain
    %if the length of Fids is odd, then you have to do a circshift of one to
    %make sure that you don't introduce a small frequency shift into the fids
    %vector.
    if mod(length(MRSCont.processed.A{kk}.specs),2)==0
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
    MRSCont.processed.A{kk}.ppm     = ppmRangeData';
    MRSCont.processed.A{kk}.specs   = specs_interp;
    MRSCont.processed.A{kk}.fids    = fids_interp;
    MRSCont.processed.A{kk}.sz      = size(MRSCont.processed.A{kk}.fids);
    MRSCont.processed.A{kk}.n       = length(MRSCont.processed.A{kk}.fids);
    MRSCont.processed.A{kk}.dwelltime = dwelltime;
    % Calculate the new spectral width and dwelltime:
    dppm                        = abs(MRSCont.processed.A{kk}.ppm(2)-MRSCont.processed.A{kk}.ppm(1));
    ppmrange                    = abs((MRSCont.processed.A{kk}.ppm(end)-MRSCont.processed.A{kk}.ppm(1)))+dppm;
    MRSCont.processed.A{kk}.spectralwidth   = ppmrange*MRSCont.processed.A{kk}.Bo*42.577;
    MRSCont.processed.A{kk}.dwelltime       = 1/MRSCont.processed.A{kk}.spectralwidth;
    % Calculate the time scale
    MRSCont.processed.A{kk}.t = (0:MRSCont.processed.A{kk}.dwelltime:(MRSCont.processed.A{kk}.sz(1)-1)*MRSCont.processed.A{kk}.dwelltime);
    
    % Calculate noise to match the wanted SNR value
    %  [MRSCont.processed.A{kk},MRSCont.in_silico.noise{kk},~]=optimize_noise(MRSCont.processed.A{kk},par.SNR(kk));
      [~,idx] = min(abs(par.SNR(kk)-SNR_lookup));
      [MRSCont.processed.A{kk},MRSCont.in_silico.noise{kk}]=op_addNoise(MRSCont.processed.A{kk},noise_sd(idx(1)));   
    
    MRSCont.raw{kk} = MRSCont.processed.A{kk};
           
end

%% Calculate QM & add Info

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
for kk = 1 : MRSCont.nDatasets
    for ss = 1 : length(SubSpec)          
        MRSCont.QM.SNR.(SubSpec{ss})(kk)    = op_getSNR(MRSCont.processed.(SubSpec{ss}){kk});       
        MRSCont.QM.FWHM.(SubSpec{ss})(kk)   = op_getLW(MRSCont.processed.(SubSpec{ss}){kk},SNRRange{ss}(1),SNRRange{ss}(2)); % in Hz       
        if ~(strcmp(SubSpec{ss},'ref') || strcmp(SubSpec{ss},'w') || strcmp(SubSpec{ss},'mm'))
             MRSCont.QM.drift.pre.(SubSpec{ss}){kk}  = 0;
            MRSCont.QM.drift.post.(SubSpec{ss}){kk} = 0;
            MRSCont.QM.freqShift.(SubSpec{ss})(kk)  = 0;       
            MRSCont.QM.res_water_amp.(SubSpec{ss})(kk) = sum(MRSCont.processed.(SubSpec{ss}){kk}.watersupp.amp);  
            if strcmp(SubSpec{ss},'diff1') || strcmp(SubSpec{ss},'sum')
                MRSCont.QM.drift.pre.(SubSpec{ss}){kk}  = reshape([MRSCont.QM.drift.pre.A'; MRSCont.QM.drift.pre.B'], [], 1)';
                MRSCont.QM.drift.post.(SubSpec{ss}){kk} = reshape([MRSCont.QM.drift.post.A'; MRSCont.QM.drift.post.B'], [], 1)';
            end
            MRSCont.QM.drift.pre.AvgDeltaCr.(SubSpec{ss})(kk) = mean(MRSCont.QM.drift.pre.(SubSpec{ss}){kk} - 3.02);
            MRSCont.QM.drift.post.AvgDeltaCr.(SubSpec{ss})(kk) = mean(MRSCont.QM.drift.pre.(SubSpec{ss}){kk} - 3.02);
        end
    end  
end

names = {'NAA_SNR','NAA_FWHM'};

QM = horzcat(MRSCont.QM.SNR.(SubSpec{1})',MRSCont.QM.FWHM.(SubSpec{1})');

MRSCont.QM.tables = array2table(QM,'VariableNames',names);
writetable(MRSCont.QM.tables,[MRSCont.outputFolder '/QM_processed_spectra.csv']);

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
if saveSpec
    [MRSCont] = osp_saveLCM(MRSCont);
    [MRSCont] = osp_saveJMRUI(MRSCont);
    [MRSCont] = osp_saveVendor(MRSCont);
end

if MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end


% function [spec,noise,step]=optimize_noise(spec,target)
%     bckp = spec;  %0.00011    
%     [spec,~]=op_addNoise(spec,0);
%     SNR    = op_getSNR(spec); 
%     step = 0;
%     while SNR > target
%         step = step + 1;
%         [spec,~]=op_addNoise(bckp,step * 0.000001  );
%         SNR    = op_getSNR(spec);       
%     end
%     [spec,noise]=op_addNoise(bckp,step * 0.000001);    
% end
