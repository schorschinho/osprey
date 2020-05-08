function F = fit_nlinwaterOsprey(x, dataToFit, basisSet, fitRangePPM)
%% F = fit_nlinwaterOsprey(x, dataToFit, basisSet, fitRangePPM)
%   This function initializes and runs fitting of an unsuppressed water
%   signal with a simulated water basis function. Parameters that are
%   applied to the basis function are zero-order phase shift, frequency
%   shift, two parameters for Gaussian/Lorentzian lineshape contribution,
%   and amplitude.
%
%   The minimization problem is solved with a Levenberg-Marquardt
%   non-linear least-squares (NLLS) optimization for the nonlinear 
%   parameters (all except amplitude), with a fast non-negative linear 
%   least-squares (FNNLS) algorithm used at every iteration of the NLLS.
%
%   USAGE:
%       [x,ampl]    = fit_waterModelNLLS(x, dataToFit, basisSet, fitRangePPM);
%
%   INPUTS:
%       data        = Vector containing the concatenated real and imaginary
%                     values of an acquired water signal.
%       basis       = Vector containing the concatenated real and imaginary
%                     values of a simulated water basis function.
%       x0          = Vector providing initial values to the non-linear
%                     least-squares solver.
%                     x0(1) = zero-order phase correction
%                     x0(2) = frequency shift
%                     x0(3) = parameter for Gaussian/Lorentzian lineshape
%                     x0(4) = parameter for Gaussian/Lorentzian lineshape
%
%   OUTPUTS:
%       x           = Optimized non-linear parameters
%       ampl        = Optimized amplitude.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-04-09)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code uses the 'fnnls' implementation of a fast non-negative
%       least-squares optimization algorithm by Rasmus Bro.
%       https://www.mathworks.com/matlabcentral/fileexchange/3388-nnls-and-constrained-regression
%       Bro et al., Journal of Chemometrics 11:393-401 (1997)
%
%   HISTORY:
%       2019-04-09: First version of the code.

% 1. Cut out the frequency range of the spectrum to be fit
dataToFit = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end));
% Turn complex-valued problem into real-valued problem
data = [real(dataToFit.specs); imag(dataToFit.specs)];

% 2. Extract parameters
nBasisFcts  = 1;
zeroPhase   = x(1); % zero-order phase correction [rad]
firstPhase  = x(2); % first-order phase correction [rad]
gaussLB     = x(3); % Gaussian damping [Hz^2]
lorentzLB   = x(4); % Lorentzian damping [Hz] for each basis function
freqShift   = x(5); % Frequency shift [Hz] for each basis function

% 3. Run the time-domain operations on the basis functions
% (frequency shift, Lorentzian damping, Gaussian damping)
t = basisSet.t;
for ii=1:nBasisFcts
    basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp(1i*freqShift.*t)' .* exp(-lorentzLB.*t)' .* exp(-gaussLB.*t.*t)';    
end
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% 4. Run the frequency-domain operations on the basis functions
% (zero and first order phase correction)
f=[(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
for ii=1:nBasisFcts
    basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp(1i*zeroPhase) .* exp(1i*firstPhase*2*pi.*f)';    
end
basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);
% Cut out the frequency range of the basis set
basisSet = op_freqrange(basisSet,fitRangePPM(1),fitRangePPM(end));


% 6. Run the first iteration of the non-negative linear solver
A = [real(basisSet.specs); imag(basisSet.specs)];
b = data;
[ampl]  = fnnls(A'*A,A'*b);

% 8. The functional to be used by the non-linear least squares solver is
% the fit residual:
F = b - (A*ampl);

end 