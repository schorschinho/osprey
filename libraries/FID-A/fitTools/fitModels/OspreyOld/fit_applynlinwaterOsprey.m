function [appliedBasisSet] = fit_applynlinwaterOsprey(basisSet, x, fitRangePPM)

% 1. Extract parameters
zeroPhase   = x(1); % zero-order phase correction [rad]
firstPhase  = x(2); % first-order phase correction [rad]
gaussLB     = x(3); % Gaussian damping [Hz^2]
lorentzLB   = x(4); % Lorentzian damping [Hz]
freqShift   = x(5); % Frequency shift [Hz]

% 2. Run the time-domain operations on the basis functions
% (frequency shift, Lorentzian damping, Gaussian damping)
t = basisSet.t;
basisSet.fids = basisSet.fids .* exp(1i*freqShift.*t)' .* exp(-lorentzLB.*t)' .* exp(-gaussLB.*t.*t)';
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% 3. Run the frequency-domain operations on the basis functions
% (zero and first order phase correction)
f=[(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
basisSet.specs = basisSet.specs .* exp(1i*zeroPhase) .* exp(1i*firstPhase*2*pi.*f)';
basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);
% Cut out the frequency range of the basis set
appliedBasisSet = op_freqrange(basisSet, fitRangePPM(1), fitRangePPM(end));

end 
