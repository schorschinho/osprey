function F = fit_nlinOspreyMultiplex(x, dataToFitInput, basisSet, spl_pos, fitRangePPM)
%   This function receives the current set of non-linear parameters during
%   every iteration of the non-linear least-squares optimization. The
%   parameters are applied to the basis set. 
%
%   Subsequently, a fast non-negative linear least-squares (FNNLS) algorithm
%   is applied to determine the optimal amplitude to the current set of
%   non-linear parameters.
%
%   This function returns the difference between the input data and
%   the current optimized model. This difference is used by the NLLS
%   solver.
%
%   USAGE:
%       F           = fit_nlinOspreyMultiplex(x, dataToFit, basisSet, spl_pos, fitRangePPM);
%
%   INPUTS:
%       data        = Vector containing the concatenated real and imaginary
%                     values of an acquired signal.
%       basis       = Vector containing the concatenated real and imaginary
%                     values of a simulated basis set.
%       x           = Vector providing initial values to the non-linear
%                     least-squares solver.
%
%   OUTPUTS:
%       F           = Difference between data and model. This is the
%                     objective to be least-squares-optimized by the 
%                     non-linear solver.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-06-04)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code uses the 'fnnls' implementation of a fast non-negative
%       least-squares optimization algorithm by Rasmus Bro.
%       https://www.mathworks.com/matlabcentral/fileexchange/3388-nnls-and-constrained-regression
%       Bro et al., Journal of Chemometrics 11:393-401 (1997)
%
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-06-04: First version of the code.


% 0. Loop over the number of data sets to be multiplex-fitted
nMultiplex = length(dataToFitInput);
nBasisFcts  = basisSet.nMets + basisSet.nMM; % number of basis functions

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
    % 2. Cut out the frequency range of the spectrum to be fit
    dataToFit = op_freqrange(dataToFitInput{rr}, fitRangePPM(1), fitRangePPM(end));
    % Turn complex-valued problem into real-valued problem
    data{rr} = [real(dataToFit.specs); imag(dataToFit.specs)];
    
    % 3. Run the time-domain operations on the basis functions
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
    
    % 4. Run the frequency-domain operations on the basis functions
    % (zero and first order phase correction)
    f=[(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
    for ii=1:nBasisFcts
        basisSet.specs(:,ii,:) = basisSet.specs(:,ii,:) .* exp(1i*zeroPhase) .* exp(1i*firstPhase*2*pi.*f)';
    end
    basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);
    % Cut out the frequency range of the basis set
    basisSet = op_freqrange(basisSet,fitRangePPM(1),fitRangePPM(end));
    
    % 5. Set up baseline spline
    base_spline{rr} = spmak(spl_pos, beta_j(:,nMultiplex)');
    B_real(:,rr) = fnval(base_spline{rr},1:1:basisSet.sz(1));
    % Get imaginary part through Hilbert transform
    B_Hilb(:,rr) = hilbert(B_real(:,rr));
    B(:,rr) = [real(B_Hilb(:,rr)); -imag(B_Hilb(:,rr))]';
   

end

% 6. Run the first iteration of the non-negative linear solver
% Concatenate data together
A = [];
b = [];
for rr = 1:nMultiplex
    A = [A; real(basisSet.specs(:,:,rr)); imag(basisSet.specs(:,:,rr))];
    b = [b; data{rr} - B(:,rr)];
end
[ampl]  = fnnls(A'*A,A'*b);

% 7. Augment the problem with soft constraints (see Wilson et al., MRM
% 2011)
% Run with the augmented data and basis set
[augmA, augmb] = fit_createSoftConstrOsprey(basisSet, A, b, ampl);
A_aug = [A; augmA];
b_aug = [b; augmb];
[ampl]  = fnnls(A_aug'*A_aug,A_aug'*b_aug);
    
% 8. The functional to be used by the non-linear least squares solver is
% the fit residual:
dataMultiplex = [];
baseMultiplex = [];
for rr = 1:nMultiplex
    dataMultiplex = [dataMultiplex; data{rr}];
    baseMultiplex = [baseMultiplex; B(:,rr)];
end
F = dataMultiplex - (A*ampl + baseMultiplex); % use if using MATLAB LevMar algorithms

end
