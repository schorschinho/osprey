% fit_OspreyMultiplex.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [fitParams] = fit_OspreyMultiplex(dataToFit, resBasisSet, fitOpts);
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

function [fitParams] = fit_OspreyMultiplex(dataToFit, resBasisSet, fitOpts)

% If only one dataset is provided, turn into cell
nMultiplex = length(dataToFit);
if nMultiplex == 1
    dataToFit = {dataToFit};
end

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
fitRangePts             = dataToFit{1}.ppm >= fitRangePPM(1) & dataToFit{1}.ppm <= fitRangePPM(2);
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
zeroPhase   = 0; % zero-order phase correction [rad]
firstPhase  = 0; % first-order phase correction [rad]
gaussLB     = op_getLW(dataToFit{1},1.8,2.2); % Gaussian damping [Hz^2]
lorentzLB   = 0.01 * ones(nBasisFcts,1); % Lorentzian damping [Hz] for each basis function
freqShift   = zeros(nBasisFcts,1); % Frequency shift [Hz] for each basis function
% We want to construct a cubic spline (i.e. of order 4); therefore a
% spline with k knots requires k-4 coefficients.
beta_j      = zeros(length(spl_pos) - 4, 1);
% We want one cubic spline for each dataset:
beta_j      = repmat(beta_j, [nMultiplex 1]);
% Concatenate together into one large x0 vector.
x0          = [zeroPhase; firstPhase; gaussLB; lorentzLB; freqShift; beta_j];
% Add additional frequency shift for each dataset
addFreqShift = zeros(nMultiplex - 1, 1);
x0           = [x0; addFreqShift];

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

% Set constraints
nMets   = resBasisSet.nMets;
nMM     = resBasisSet.nMM;
lb      = zeros(size(x0));
ub      = zeros(size(x0));
lb(1,:)   = -2*pi; 
ub(1,:)   = 2*pi; % Zero order phase shift [rad]
lb(2,:)   = -pi/4; 
ub(2,:)   = pi/4; % First order phase shift [rad]
lb(3,:)   = 0; 
ub(3,:)   = 5000; % Gaussian linebroadening [Hz^2]
lb(4:nMets+3,:) = 0; 
ub(4:nMets+3,:) = 10; % Lorentzian linebroadening [Hz] - Metabolites
lb(nMets+4:nMets+nMM+3,:) = 0; 
ub(nMets+4:nMets+nMM+3,:) = 50; % Lorentzian linebroadening [Hz] - MM/Lipids
lb(nMets+nMM+4:2*nMets+nMM+3,:) = -4; 
ub(nMets+nMM+4:2*nMets+nMM+3,:) = 4; % Frequency shift [Hz] - Metabolites
lb(2*nMets+nMM+4:2*nMets+2*nMM+3,:) = -6.5;
ub(2*nMets+nMM+4:2*nMets+2*nMM+3,:) = 6.5; % Frequency shift [Hz] - MM/Lipids
lb(2*nMets+2*nMM+4:end-(nMultiplex-1),:) = -1; 
ub(2*nMets+2*nMM+4:end-(nMultiplex-1),:) = 1; % Baseline coefficients
lb(end-(nMultiplex-2):end,:) = -1;
ub(end-(nMultiplex-2):end,:) = 1; % additional frequency shift [Hz] - subspectra

% Run non-linear solver
[x,~,~,~,~,~,~,~] = LevenbergMarquardt(@(x) fit_nlinOspreyMultiplex(x, dataToFit, resBasisSet, spl_pos, fitRangePPM),x0,lb,ub,opts);


%%% 4. PERFORM FINAL COMPUTATION OF LINEAR PARAMETERS %%%
% After the non-linear optimization is finished, we need to perform the
% final evaluation of the linear parameters (i.e. the amplitudes).
% Apply the parameters from the final iteration of the non-linear solver to
% the basis set to create the phased&broadened basis functions, as well as
% the baseline.
[appliedBasisSet, B] = fit_applynlinOspreyMultiplex(resBasisSet, x, spl_pos, fitRangePPM);

% Run the final iteration of the non-negative linear solver
% Concatenate data together
A = [];
b = [];
for rr = 1:nMultiplex
    % Cut out the frequency range of the spectrum to be fit
    dataToFit{rr} = op_freqrange(dataToFit{rr}, fitRangePPM(1), fitRangePPM(end));
    % Turn complex-valued problem into real-valued problem
    data{rr} = [real(dataToFit{rr}.specs); imag(dataToFit{rr}.specs)];
    A = [A; real(appliedBasisSet.specs(:,:,rr)); imag(appliedBasisSet.specs(:,:,rr))];
    b = [b; data{rr} - B(:,rr)];
end
[ampl]  = fnnls(A'*A,A'*b);

% Augment the problem with soft constraints (see Wilson et al., MRM
% 2011)
% Run with the augmented data and basis set
[augmA, augmb] = fit_createSoftConstrOsprey(appliedBasisSet, A, b, ampl);
A_aug = [A; augmA];
b_aug = [b; augmb];
[ampl]  = fnnls(A_aug'*A_aug,A_aug'*b_aug);

% Store all fit parameters in a structure
fitParams.ampl      = ampl;
fitParams.ph0       = x(1);
fitParams.ph1       = x(2);
fitParams.gaussLB   = x(3);
fitParams.lorentzLB = x(4:nBasisFcts+3);
fitParams.freqShift = x(nBasisFcts+4:2*nBasisFcts+3);
fitParams.beta_j    = x(2*nBasisFcts+4:end-(nMultiplex-1));
fitParams.addFreqShift = x(end-(nMultiplex-2):end);
fitParams.spl_pos   = spl_pos;
fitParams.refShift = 0;
fitParams.refFWHM = 0;
fitParams.prelimParams = fitParams;

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
