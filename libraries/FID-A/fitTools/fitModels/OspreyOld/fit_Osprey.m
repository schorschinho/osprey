% fit_Osprey.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [fitParams] = fit_Osprey(dataToFit, resBasisSet, fitOpts);
% 
% DESCRIPTION:
% This is the function describing the Osprey fitting model, which is
% performed with the resampled basis set resBasisSet on the data contained
% in the FID-A data structure dataToFit, using the fit options in the
% structure fitOpts.
% 
% OUTPUTS:
% fitParams   = Structure containing all necessary fit parameters
%               x = concatenated vector with zero- and first-order phase
%                   shift, Gaussian linebroadening, Lorentzian
%                   linebroadening, frequency shifts, and cubic spline 
%                   baseline coefficients
%               spl_pos = positions (in points) of the baseline spline
%                   knots
%               ampl = amplitudes of the basis functions
%
% INPUTS:
% dataToFit   = FID-A data structure
% basisSet    = FID-A basis set container
% fitOpts     = Structure containing fit options

function [fitParams, resBasisSet] = fit_Osprey(dataToFit, basisSet, fitOpts)

%%% 1. PREPARE BASIS SET %%%
% Resample basis set to match data resolution and frequency range
resBasisSet             = fit_resampleBasis(dataToFit, basisSet);

%%% 1. EXTRACT OPTIONS AND PREPARE FIT %%%
% Extract ppm fit range
fitRangePPM             = fitOpts.range;
% Initialize the baseline spline parameters
knotSpacingPPM          = fitOpts.bLineKnotSpace;
% Determine locations of baseline knots (in ppm)
fitOpts.knotLocations   = (fitRangePPM(1):knotSpacingPPM:fitRangePPM(end));
% Determine the number of baseline knots
N_b                     = length(fitOpts.knotLocations);
% Determine the knot spacing and their positions (in data pts)
fitRangePts             = dataToFit.ppm >= fitRangePPM(1) & dataToFit.ppm <= fitRangePPM(2);
nFitRangePts            = sum(fitRangePts);
knot_spacing            = (nFitRangePts-1) / (N_b-1); % knot spacing in pts
spl_pos                 = [1:knot_spacing:nFitRangePts]; % vector with knot position in pts


%%% 2. SET AND GET STARTING VALUES %%%
% Set the starting values for the non-linear parameters to be passed on
% to the Levenberg-Marquardt NLLS solving algorithm. Note that the
% amplitude parameters do not need to be initialized, as the
% non-negative linear Lawson-Hanson solver does not require starting
% values to be provided.
nBasisFcts  = resBasisSet.nMets + resBasisSet.nMM; % number of basis functions
ph0   = 0; % zero-order phase correction [rad]
ph1  = 0; % first-order phase correction [rad]
gaussLB     = op_getLW(dataToFit,1.8,2.2); % Gaussian damping [Hz^2]
lorentzLB   = 0.01 * ones(nBasisFcts,1); % Lorentzian damping [Hz] for each basis function
freqShift   = zeros(nBasisFcts,1); % Frequency shift [Hz] for each basis function
% We want to construct a cubic spline (i.e. of order 4); therefore a
% spline with k knots requires k-4 coefficients.
beta_j      = zeros(length(spl_pos) - 4, 1);
% Concatenate together into one large x0 vector.
x0          = [ph0; ph1; gaussLB; lorentzLB; freqShift; beta_j];


%%% 3. CALL NON-LINEAR SOLVER %%%
% Run the non-linear solver to optimize the non-linear parameters. At each
% iteration of the non-linear solver, the linear parameters are calculated
% with a non-negative linear solver.

% a. MATLAB-internal Levenberg-Marquardt (does not allow constraints)
% lb = [];
% ub = [];
% options = optimset('lsqnonlin');
% options = optimset(options, 'Display', 'off', 'FunValCheck', 'on', 'Algorithm', 'levenberg-marquardt');
% [x,~] = lsqnonlin(@(x) fit_metabModelNLLS(x, dataToFit, resBasisSet, spl_pos, fitRangePPM), ...
%     x0,lb,ub,options);

% b. Levenberg-Marquardt implementation allowing parameter constraints
% (c) Alexander Drentler
% https://www.mathworks.com/matlabcentral/fileexchange/53449-levenberg-marquardt-toolbox
opts.Display            = 'notify';
opts.Jacobian           = 'limit';
%   lets get rid of all stopping criteria
opts.FooTol             = NaN; 
opts.RelFooTol          = NaN;
opts.FooTol             = NaN;
opts.RelTolX            = NaN;
opts.TolX               = 1e-10;
%   note that we need to square the acceptance tolerance
abstol                  = 1e-8;
opts.AccTol             = abstol^2;
opts.MaxDamping         = 1e20;
opts.Broyden_updates    = 'off';
opts.FinDiffRelStep=eps^(1/3);
opts.TypicalX = 1;
opts.ScaleProblem = 'Jacobian';

% Set constraints
nMets   = resBasisSet.nMets;
nMM     = resBasisSet.nMM;
lb      = zeros(size(x0));
ub      = zeros(size(x0));
lb(1)   = -2*pi; 
ub(1)   = 2*pi; % Zero order phase shift [rad]
lb(2)   = -pi/4; 
ub(2)   = pi/4; % First order phase shift [rad]
lb(3)   = 0; 
ub(3)   = 5000; % Gaussian linebroadening [Hz^2]
lb(4:nMets+3) = 0; 
ub(4:nMets+3) = 10; % Lorentzian linebroadening [Hz] - Metabolites
lb(nMets+4:nMets+nMM+3) = 0; 
ub(nMets+4:nMets+nMM+3) = 50; % Lorentzian linebroadening [Hz] - MM/Lipids
lb(nMets+nMM+4:2*nMets+nMM+3) = -4; 
ub(nMets+nMM+4:2*nMets+nMM+3) = 4; % Frequency shift [Hz] - Metabolites
lb(2*nMets+nMM+4:2*nMets+2*nMM+3) = -6.5;
ub(2*nMets+nMM+4:2*nMets+2*nMM+3) = 6.5; % Frequency shift [Hz] - MM/Lipids
lb(2*nMets+2*nMM+4:end) = -1; 
ub(2*nMets+2*nMM+4:end) = 1; % Baseline coefficients

% Run non-linear solver
[x,~,~,~,~,~,~,~] = LevenbergMarquardt(@(x) fit_Osprey_Model(x, dataToFit, resBasisSet, spl_pos, fitRangePPM),x0,lb,ub,opts);


%%% 4. PERFORM FINAL COMPUTATION OF LINEAR PARAMETERS %%%
% After the non-linear optimization is finished, we need to perform the
% final evaluation of the linear parameters (i.e. the amplitudes).
% Apply the parameters from the final iteration of the non-linear solver to
% the basis set to create the phased&broadened basis functions, as well as
% the baseline.
[ampl, x] = fit_Osprey_FinalLinear(dataToFit, resBasisSet, x, spl_pos, fitRangePPM);



% Store all fit parameters in a structure
fitParams.ampl      = ampl;
fitParams.ph0       = x(1);
fitParams.ph1       = x(2);
fitParams.gaussLB   = x(3);
fitParams.lorentzLB = x(4:nBasisFcts+3);
fitParams.freqShift = x(nBasisFcts+4:2*nBasisFcts+3);
fitParams.beta_j    = x(2*nBasisFcts+4:end);
fitParams.spl_pos   = spl_pos;

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% MODEL FUNCTION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = fit_Osprey_Model(x, dataToFit, basisSet, spl_pos, fitRangePPM)
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
%       F           = fit_nlinOsprey(x, dataToFit, basisSet, spl_pos, fitRangePPM);
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

% 1. Cut out the frequency range of the spectrum to be fit
dataToFit = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end));
% Turn complex-valued problem into real-valued problem
data = [real(dataToFit.specs); imag(dataToFit.specs)];

% 2. Extract parameters
nBasisFcts  = basisSet.nMets + basisSet.nMM; % number of basis functions
ph0         = x(1); % zero-order phase correction [rad]
ph1         = x(2); % first-order phase correction [rad]
gaussLB     = x(3); % Gaussian damping [Hz^2]
lorentzLB   = x(4:nBasisFcts+3); % Lorentzian damping [Hz] for each basis function
freqShift   = x(nBasisFcts+4:2*nBasisFcts+3); % Frequency shift [Hz] for each basis function
% We want to construct a cubic spline (i.e. of order 4); therefore a
% spline with k knots requires k-4 coefficients.
beta_j      = x(2*nBasisFcts+4:end);

% 3. Run the time-domain operations on the basis functions
% (frequency shift, Lorentzian damping, Gaussian damping)
t = basisSet.t;
for ii=1:nBasisFcts
    basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp(1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)';    
end
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% 4. Run the frequency-domain operations on the basis functions
% (zero and first order phase correction)
f=[(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
for ii=1:nBasisFcts
    basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp(1i*ph0) .* exp(1i*ph1*2*pi.*f)';    
end
basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);
% Cut out the frequency range of the basis set
basisSet = op_freqrange(basisSet,fitRangePPM(1),fitRangePPM(end));

% 5. Set up baseline spline
base_spline = spmak(spl_pos, beta_j');
B = fnval(base_spline,1:1:basisSet.sz(1));
% Get imaginary part through Hilbert transform
B_Hilb = hilbert(B);
B = [real(B_Hilb) -imag(B_Hilb)]';

% 6. Run the first iteration of the non-negative linear solver
A = [real(basisSet.specs); imag(basisSet.specs)];
b = data - B;
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
F = data - (A*ampl + B); % use if using MATLAB LevMar algorithms

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FINAL LINEAR ITERATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ampl, x] = fit_Osprey_FinalLinear(dataToFit, basisSet, x, spl_pos, fitRangePPM)

% 1. Extract parameters
nMets       = basisSet.nMets;
nMM         = basisSet.nMM;
nBasisFcts  = nMets + nMM; % number of basis functions
ph0         = x(1); % zero-order phase correction [rad]
ph1         = x(2); % first-order phase correction [rad]
gaussLB     = x(3); % Gaussian damping [Hz^2]
lorentzLB   = x(4:nBasisFcts+3); % Lorentzian damping [Hz] for each basis function
freqShift   = x(nBasisFcts+4:2*nBasisFcts+3); % Frequency shift [Hz] for each basis function
% We want to construct a cubic spline (i.e. of order 4); therefore a
% spline with k knots requires k-4 coefficients.
beta_j      = x(2*nBasisFcts+4:end);

% 2. Run the time-domain operations on the basis functions
% (frequency shift, Lorentzian damping, Gaussian damping)
t = basisSet.t;
for ii=1:nBasisFcts
    basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp(1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)';     
end
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% 3. Run the frequency-domain operations on the basis functions
% (zero and first order phase correction)
f=[(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
for ii=1:nBasisFcts
    basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp(1i*ph0) .* exp(1i*ph1*2*pi.*f)';    
end
basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);
% Cut out the frequency range of the basis set
appliedBasisSet = op_freqrange(basisSet, fitRangePPM(1), fitRangePPM(end));

% 4. Set up baseline spline
base_spline = spmak(spl_pos, beta_j');
B = fnval(base_spline,1:1:appliedBasisSet.sz(1));
% Get imaginary part through Hilbert transform
B_Hilb = hilbert(B);
B = [real(B_Hilb) -imag(B_Hilb)]';

% Run the final iteration of the non-negative linear solver
% Cut out the frequency range of the spectrum to be fit
dataToFit = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end));
% Turn complex-valued problem into real-valued problem
data = [real(dataToFit.specs); imag(dataToFit.specs)];
A = [real(appliedBasisSet.specs); imag(appliedBasisSet.specs)];
b = data - B;
[ampl]  = fnnls(A'*A,A'*b);

% Augment the problem with soft constraints (see Wilson et al., MRM
% 2011)
% Run with the augmented data and basis set
[augmA, augmb] = fit_createSoftConstrOsprey(appliedBasisSet, A, b, ampl);
A_aug = [A; augmA];
b_aug = [b; augmb];
[ampl]  = fnnls(A_aug'*A_aug,A_aug'*b_aug);


end 
