function [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams)
%   This function applies a set of Osprey fit parameters to a basis set of
%   metabolite basis functions, and a set of cubic baseline spline basis
%   functions.
%
%   The function returns the complete fit, the baseline, the residual,
%   individual contributions from each basis function, the data that are
%   being fit, and the frequency axis.
%
%   USAGE:
%       [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams)
%
%   INPUTS:
%       inputData   = Struct containing all necessary data:
%                       - Original (un-zerofilled data)
%                       - Basis set (resampled to data res)
%       inputSettings = Struct containing all necessary settings:
%                       - fit range [ppm]
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
basisSet      = inputData.basisSet;
% ... settings:
fitRangePPM         = inputSettings.fitRangePPM;
scale               = inputSettings.scale;
% ... fit parameters
nMets       = basisSet.nMets;
nMM         = basisSet.nMM;
nBasisFcts  = nMets + nMM; % number of basis functions
ph0         = fitParams.ph0; % zero-order phase correction [convert from deg to rad]
ph1         = fitParams.ph1; % first-order phase correction [convert from deg/ppm to rad/ppm]
gaussLB     = fitParams.gaussLB; % Gaussian damping [Hz]
lorentzLB   = fitParams.lorentzLB; % Lorentzian damping [Hz] for each basis function
freqShift   = fitParams.freqShift; % Frequency shift [Hz] for each basis function
ampl        = fitParams.ampl; % Amplitudes for metabolite/MM/lipid basis functions
beta_j      = fitParams.beta_j; % Amplitudes for baseline spline basis functions
spl_pos     = fitParams.spl_pos; %spline positions
% ... is concatenated
if strcmp(inputSettings.fitStyle,'Concatenated')
    n_beta_j = length(beta_j);
    if inputSettings.flags.isMEGA
        switch inputSettings.concatenated.Subspec
            case 'diff1'
                basisSet.fids = basisSet.fids(:,:,1);
                beta_j = beta_j(1:n_beta_j/2);
            case 'sum'
                basisSet.fids = basisSet.fids(:,:,2);
                beta_j = beta_j(1+n_beta_j/2:end);
        end
    end
    if (inputSettings.flags.isHERMES || inputSettings.flags.isHERCULES)
        switch inputSettings.concatenated.Subspec
            case 'diff1'
                basisSet.fids = basisSet.fids(:,:,1);
                beta_j = beta_j(1:n_beta_j/3);
            case 'diff2'
                basisSet.fids = basisSet.fids(:,:,2);
                beta_j = beta_j(1+n_beta_j/3:2*n_beta_j/3);
            case 'sum'
                basisSet.fids = basisSet.fids(:,:,3);
                beta_j = beta_j(1+2*n_beta_j/3:end);
        end
    end
end

%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
% Run the time-domain operations on the basis functions
% (frequency shift, Lorentzian damping, Gaussian damping)
t = basisSet.t;
for ii=1:nBasisFcts
    basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp(1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)';     
end
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% Run the frequency-domain operations on the basis functions
% (zero and first order phase correction)
f=[(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
for ii=1:nBasisFcts
    basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp(1i*ph0) .* exp(1i*ph1*2*pi.*f)';    
end
basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);
% Cut out the frequency range of the basis set
appliedBasisSet = op_freqrange(basisSet, fitRangePPM(1), fitRangePPM(end));


%%% 3. APPLY THE LINEAR PARAMETERS %%%
% Set up baseline spline
base_spline = spmak(spl_pos, beta_j');
B = fnval(base_spline,1:1:appliedBasisSet.sz(1));
% Get imaginary part through Hilbert transform
B_Hilb = hilbert(B);
baseline = [real(B_Hilb)]';

% Calculate the final fit
A = [real(appliedBasisSet.specs)];
completeFit = A * ampl + baseline;

% Calculate the residual
% Cut out the frequency range of the spectrum to be fit
dataToFit   = op_ampScale(dataToFit, 1/scale);
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end));
% Use only the real part to fit here
data        = [real(dataToFit.specs)]; % data
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
