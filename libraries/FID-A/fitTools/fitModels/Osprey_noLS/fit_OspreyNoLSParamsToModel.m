function [ModelOutput] = fit_OspreyNoLSParamsToModel(inputData, inputSettings, fitParams)
%   This function applies a set of LCModel fit parameters to a basis set of
%   metabolite basis functions, and a set of cubic baseline spline basis
%   functions.
%
%   The function returns the complete fit, the baseline, the residual,
%   individual contributions from each basis function, the data that are
%   being fit, and the frequency axis.
%
%   USAGE:
%       [ModelOutput] = fit_plotParamsToModel(inputData, inputSettings, fitParams)
%
%   INPUTS:
%       inputData   = Struct containing all necessary data:
%                       - Original (un-zerofilled data)
%                       - Basis set (resampled to zero-filled data res)
%       inputSettings = Struct containing all necessary settings:
%                       - fit range [ppm]
%                       - minimum knot spacing [ppm]
%                       - scale (scaling factor applied to data prior to
%                           fitting)
%       fitParams    = Struct containing the set of Osprey model parameters
%
%   OUTPUTS:
%       ModelOutput = Struct containing:
%                       - Complete fit
%                       - Baseline
%                       - Residual
%                       - Individual basis function contributions
%                       - Original data that are being fit
%                       - ppm axis of the fit range


%%% 1. UNPACK THE INPUT %%%
% ... data:
dataToFit     = inputData.dataToFit;
dataToFit     = op_zeropad(dataToFit, 2); % zero-fill for LCModel
basisSet      = inputData.basisSet;
if (length(fitParams.ampl) == 3)
    basisSet      = inputData.basisSet_mm;
    fitParams.freqShift = repmat(fitParams.freqShift,[basisSet.nMets+basisSet.nMM 1]);
    fitParams.lorentzLB = repmat(fitParams.lorentzLB,[basisSet.nMets+basisSet.nMM 1]);
    %dummy=fitParams.ampl;   
    %fitParams.ampl=zeros([basisSet.nMets+basisSet.nMM 1]);
    %fitParams.ampl([3 4 13])=dummy;
end
% ... settings:
fitRangePPM         = inputSettings.fitRangePPM;
minKnotSpacingPPM   = inputSettings.minKnotSpacingPPM;
scale               = inputSettings.scale;
% ... fit parameters
nMets       = basisSet.nMets;
nMM         = basisSet.nMM;
nBasisFcts  = nMets + nMM; % number of basis functions
ph0         = fitParams.ph0; % zero-order phase correction [convert from deg to rad]
gaussLB     = fitParams.gaussLB; % Gaussian damping [Hz]
asym        = fitParams.asym; % Assymetry factor
lorentzLB   = fitParams.lorentzLB; % Lorentzian damping [Hz] for each basis function
freqShift   = fitParams.freqShift; % Frequency shift [Hz] for each basis function
freqShiftmets   = fitParams.freqShiftmets; % Frequency shift [Hz] for each basis function
ampl        = fitParams.ampl; % Amplitudes for metabolite/MM/lipid basis functions
beta_j      = fitParams.beta_j; % Amplitudes for baseline spline basis functions
refShift    = fitParams.refShift; % Reference shift applied to the data during first step of fitting
ph1 = fitParams.ph1;

% Create an array of normalized cubic baseline spline basis functions.
[splineArray, ~]    = fit_makeSplineBasis(dataToFit, fitRangePPM, 1/15);

%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
f = [(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(4*basisSet.sz(1))):basisSet.spectralwidth/(2*basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(4*basisSet.sz(1)))];
asym_fwhm = 2 * gaussLB ./ (1 + exp(-asym .* f));
G = sqrt(4 .* log(2)./pi./(asym_fwhm.^2)) .*  exp(-4 * log(2) * (f./asym_fwhm).^2);
G = ifft(fftshift(G',1),[],1);
G = G(1:basisSet.sz(1));
G = G/real(G(1));

% Run the time-domain operations on the metabolite basis functions
% (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
t = basisSet.t;
for ii=1:nBasisFcts
    basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp(-1i*freqShiftmets(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* G;   
end
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% Run the frequency-domain operations on the basis functions
% (first order phase correction)
% Cut out the frequency range of the basis set
basisSet = op_freqrange(basisSet,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));


% Apply phasing to the spline basis functions
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
B = [real(B)];


%%% 3. APPLY THE LINEAR PARAMETERS %%%
% Convolve the lineshape with the metabolite basis functions only
% (NOT the macromolecules or lipids or baseline splines).
A = real(basisSet.specs);

% Calculate the final baseline
baseline    = B * beta_j;
completeFit = A * ampl + baseline;

% Calculate the residual
% Cut out the frequency range of the spectrum to be fit
% Apply initial referencing shift
dataToFit   = op_freqshift(dataToFit, -refShift);
dataToFit=op_addphase(dataToFit,ph0,ph1);
dataToFit=op_freqshift(dataToFit,freqShift);
dataToFit   = op_ampScale(dataToFit, 1/scale);
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end),length(splineArray(:,1,1)));


% Use only the real part to fit here
data        = real(dataToFit.specs); % data
ppm         = dataToFit.ppm;
residual    = data - completeFit;


%%% 4. CREATE OUTPUT %%%
% Return a struct with all the output parameters
ModelOutput.ppm             = ppm;
ModelOutput.data            = data;
ModelOutput.completeFit     = completeFit;
ModelOutput.baseline        = baseline;
ModelOutput.residual        = residual;
for kk = 1:nBasisFcts
    ModelOutput.indivMets(:,kk) = ampl(kk) * A(:,kk);
end


end
