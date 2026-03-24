function [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams)
%   This function applies a set of Osprey fit parameters to a water basis
%   set.
%
%   The function returns the complete fit, the residual, the water data
%   that are being fit, and the frequency axis.
%
%   USAGE:
%       [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams)
%
%   INPUTS:
%       inputData   = Struct containing all necessary data:
%                       - Original (un-zerofilled data)
%                       - Basis set (resampled to data res)
%       inputSettings = Struct containing all necessary settings:
%                       - water fit range [ppm]
%                       - scale (scaling factor applied to data prior to
%                           fitting)
%       fitParams    = Struct containing the set of Osprey model parameters
%
%   OUTPUTS:
%       ModelOutput = Struct containing:
%                       - Complete water fit
%                       - Residual
%                       - Original water data that are being fit
%                       - ppm axis of the fit range


%%% 1. UNPACK THE INPUT %%%
% ... data:
dataToFit     = inputData.dataToFit;
basisSet      = inputData.basisSet;
% ... settings:
fitRangePPMWater    = inputSettings.fitRangePPM;
scale               = inputSettings.scale;
% ... fit parameters
nBasisFcts  = basisSet.nMets;
ph0         = fitParams.ph0; % zero-order phase correction [convert from deg to rad]
ph1         = fitParams.ph1; % first-order phase correction [convert from deg/ppm to rad/ppm]
gaussLB     = fitParams.gaussLB; % Gaussian damping [Hz]
lorentzLB   = fitParams.lorentzLB; % Lorentzian damping [Hz] for water
freqShift   = fitParams.freqShift; % Frequency shift [Hz] for water
ampl        = fitParams.ampl; % Amplitude for water


%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
% Run the time-domain operations on the basis functions
% (frequency shift, Lorentzian damping, Gaussian damping)
t = basisSet.t;
for ii=1:nBasisFcts
    basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp(1i*freqShift.*t)' .* exp(-lorentzLB.*t)' .* exp(-gaussLB.*t.*t)';     
end
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% Run the frequency-domain operations on the basis functions
% (zero and first order phase correction)
f=[(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
for ii=1:nBasisFcts
    basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp(1i*ph0) .* exp(1i*ph1*2*pi.*f)';    
end
basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);

% Cut out the frequency range of the spectrum to be fit
dataToFit   = op_ampScale(dataToFit, 1/scale);
dataToFit   = op_freqrange(dataToFit, fitRangePPMWater(1), fitRangePPMWater(end));

% Cut out the frequency range of the basis set
basisSet = op_freqrange(basisSet,fitRangePPMWater(1),fitRangePPMWater(end),dataToFit.sz(1));


%%% 3. APPLY THE LINEAR PARAMETERS %%%
% Calculate the final fit
A = [real(basisSet.specs)];
completeFit = A * ampl;

% Calculate the residual
% Use only the real part to fit here
data        = [real(dataToFit.specs)]; % data
ppm         = dataToFit.ppm;
residual    = data - completeFit;


%%% 4. CREATE OUTPUT %%%
% Return a struct with all the output parameters
ModelOutput.ppm             = ppm;
ModelOutput.data            = data;
ModelOutput.completeFit     = completeFit;
ModelOutput.residual        = residual;


end
