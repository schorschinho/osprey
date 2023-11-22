function [basisSim] = makeMMLipBasis(basisSet, MMLipConfig, DataToModel)
% Georg Oeltzschner, Johns Hopkins University 2023
% Creates MM and lipid basis functions according to the 'MMLipConfig'
% struct (provided by MM/lipid definition JSON), matching the metabolite
% basis functions provided in 'basisSet'

% Duplicate basis struct
basisSim        = basisSet;
basisSim.fids   = [];
basisSim.specs  = [];
basisSim.nMets  = 0;
basisSim.nMM    = 0;
basisSim.sz     = [];
basisSim.name   = {};

% Obtain area of water to calculate the one-proton scaling factor
idx = find(contains(basisSet.name,'H2O'));
integralH2O = sum(real(basisSet.specs(:,idx,1)));
oneProtonArea = integralH2O/2;

% Parse MM/lipid definition
basisSimFieldNames = fieldnames(MMLipConfig);
for ff = 1:length(basisSimFieldNames)
    basisSim.name{ff} = basisSimFieldNames{ff};
    simPPM = MMLipConfig.(basisSimFieldNames{ff}).ppm;
    simFW  = MMLipConfig.(basisSimFieldNames{ff}).fwmin;
    simAmp = MMLipConfig.(basisSimFieldNames{ff}).amp;
    temp = [];
    for ss = 1:length(simPPM)
        hzPerPPM = lookUpGyromagRatioForNucleus(DataToModel.nucleus) * basisSim.Bo;
        % convert FWHM to Gaussian sigma
        targetSigma = 2*pi*(simFW(ss) * hzPerPPM) / (2*sqrt(2*log(2)));
        freqShift = (-basisSim.centerFreq + simPPM(ss)) * hzPerPPM;
        amp = simAmp(ss);
        temp(:,ss) = amp .* exp(-1i*2*pi*freqShift.*basisSim.t)';
        temp(:,ss) = temp(:,ss) .* exp(-basisSim.t' .* basisSim.t' .* (targetSigma.^2)/2);
        temp(1,ss) = 0.5*temp(1,ss); %first-point correction (measured for only half the dwell-time)

        % Calculate actual area
        integralTemp = sum(real(fftshift(fft(temp(:,ss)))));
        scaleFac = (amp.*oneProtonArea)./integralTemp;
        temp(:,ss) = scaleFac.*temp(:,ss);
    end
    basisSim.fids(:,ff) = sum(temp,2);
end
basisSim.specs = fftshift(fft(basisSim.fids,[],1),1);
basisSim.sz = size(basisSim.fids);
basisSim.nMM = length(basisSimFieldNames);
end