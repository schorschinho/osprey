% fit_waterOsprey.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [fitParams] = fit_waterOsprey(dataToFit, resBasisSet, fitOpts);
% 
% DESCRIPTION:
% This is the function describing the Osprey water fitting model, which is
% performed with the resampled water basis set resBasisSet on the data contained
% in the FID-A water data structure dataToFit, using the fit options in the
% structure fitOpts.
% 
% OUTPUTS:
% fitParams   = Structure containing all necessary fit parameters
%               x = concatenated vector with zero- and first-order phase
%                   shift, Gaussian linebroadening, Lorentzian
%                   linebroadening, and frequency shift
%               ampl = amplitude of the water basis function
%
% INPUTS:
% dataToFit   = FID-A data structure
% basisSet    = FID-A basis set container
% fitOpts     = Structure containing fit options

function [fitParams] = fit_waterOsprey(dataToFit, resBasisSet, fitOpts)


%%% 1. EXTRACT OPTIONS AND PREPARE FIT %%%
% Extract ppm fit range
fitRangePPMWater = fitOpts.rangeWater;


%%% 2. SET STARTING VALUES %%%
% Set the starting values for the non-linear parameters to be passed on
% to the Levenberg-Marquardt NLLS solving algorithm. Note that the
% amplitude parameter does not need to be initialized, as the
% non-negative linear Lawson-Hanson solver does not require starting
% values to be provided.
zeroPhase   = 0; % zero-order phase correction [rad]
firstPhase  = 0; % first-order phase correction [rad]
gaussLB     = op_getLW(dataToFit,4.5,4.9); % Gaussian damping [Hz^2]
lorentzLB   = 0.01; % Lorentzian damping [Hz] for water
freqShift   = 0; % Frequency shift [Hz] for water
% Concatenate together into one large x0 vector.
x0          = [zeroPhase; firstPhase; gaussLB; lorentzLB; freqShift];


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
lb      = zeros(size(x0));
ub      = zeros(size(x0));
lb(1)   = -2*pi; 
ub(1)   = 2*pi; % Zero order phase shift [rad]
lb(2)   = -pi/4; 
ub(2)   = pi/4; % First order phase shift [rad]
lb(3)   = 0; 
ub(3)   = 5000; % Gaussian linebroadening [Hz^2]
lb(4) = 0; 
ub(4) = 50; % Lorentzian linebroadening [Hz]
lb(5) = -15; % Frequency shift [Hz]
ub(5) = 15; % Frequency shift [Hz]

% Run non-linear solver
[x,~,~,~,~,~,~,~] = LevenbergMarquardt(@(x) fit_nlinwaterOsprey(x, dataToFit, resBasisSet, fitRangePPMWater),x0,lb,ub,opts);


%%% 4. PERFORM FINAL COMPUTATION OF LINEAR PARAMETERS %%%
% After the non-linear optimization is finished, we need to perform the
% final evaluation of the linear parameter (i.e. the water amplitude).
% Apply the parameters from the final iteration of the non-linear solver to
% the basis set to create the phased&broadened water basis function.
[appliedBasisSet] = fit_applynlinwaterOsprey(resBasisSet, x, fitRangePPMWater);

% Run the final iteration of the non-negative linear solver
% Cut out the frequency range of the spectrum to be fit
dataToFit = op_freqrange(dataToFit, fitRangePPMWater(1), fitRangePPMWater(end));
% Turn complex-valued problem into real-valued problem
data = [real(dataToFit.specs); imag(dataToFit.specs)];
A = [real(appliedBasisSet.specs); imag(appliedBasisSet.specs)];
b = data;
[ampl]  = fnnls(A'*A,A'*b);

% Store all fit parameters in a structure
fitParams.ampl      = ampl;
fitParams.ph0       = x(1);
fitParams.ph1       = x(2);
fitParams.gaussLB   = x(3);
fitParams.lorentzLB = x(4);
fitParams.freqShift = x(5);

end

