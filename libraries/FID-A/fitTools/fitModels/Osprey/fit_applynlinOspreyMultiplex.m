function [appliedBasisSet, B] = fit_applynlinOspreyMultiplex(basisSet, x, spl_pos, fitRangePPM)

% 1. Extract parameters
nMets       = basisSet.nMets;
nMM         = basisSet.nMM;
nMultiplex  = size(basisSet.fids, 3);
nBasisFcts  = nMets + nMM; % number of basis functions

% 1. Extract parameters
zeroPhase   = x(1); % zero-order phase correction [rad]
firstPhase  = x(2); % first-order phase correction [rad]
gaussLB     = x(3); % Gaussian damping [Hz^2]
lorentzLB   = x(4:nBasisFcts+3); % Lorentzian damping [Hz] for each basis function
freqShift   = x(nBasisFcts+4:2*nBasisFcts+3); % Frequency shift [Hz] for each basis function
% We want to construct a cubic spline (i.e. of order 4); therefore a
% spline with k knots requires k-4 coefficients.
beta_j      = x(2*nBasisFcts+4:end-(nMultiplex-1));
beta_j      = reshape(beta_j, [length(beta_j)/nMultiplex nMultiplex]);
% Additional frequency shift for each dataset
addFreqShift = x(end-(nMultiplex-2):end);

for rr = 1:nMultiplex
    % 2. Run the time-domain operations on the basis functions
    % (frequency shift, Lorentzian damping, Gaussian damping)
    t = basisSet.t;
    for ii=1:nBasisFcts
        basisSet.fids(:,ii,rr) = basisSet.fids(:,ii,rr) .* exp(1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)';
        if rr > 1
            % Add the same additional frequency shift (same for all
            % datasets), if multiplexed
            basisSet.fids(:,ii,rr) = basisSet.fids(:,ii,rr) .* exp(1i*addFreqShift(rr-1).*t)';
        end
    end
    basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);
    
    % 3. Run the frequency-domain operations on the basis functions
    % (zero and first order phase correction)
    f=[(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
    for ii=1:nBasisFcts
        basisSet.specs(:,ii,:) = basisSet.specs(:,ii,:) .* exp(1i*zeroPhase) .* exp(1i*firstPhase*2*pi.*f)';
    end
    basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);
    % Cut out the frequency range of the basis set
    appliedBasisSet = op_freqrange(basisSet,fitRangePPM(1),fitRangePPM(end));
    
    % 4. Set up baseline spline
    base_spline{rr} = spmak(spl_pos, beta_j(:,nMultiplex)');
    B_real(:,rr) = fnval(base_spline{rr},1:1:appliedBasisSet.sz(1));
    % Get imaginary part through Hilbert transform
    B_Hilb(:,rr) = hilbert(B_real(:,rr));
    B(:,rr) = [real(B_Hilb(:,rr)); -imag(B_Hilb(:,rr))]';
   

end

end 
