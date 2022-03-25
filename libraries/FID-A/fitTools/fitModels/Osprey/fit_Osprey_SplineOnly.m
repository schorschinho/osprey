function [fitParamsStep3,out] = fit_Osprey_SplineOnly(dataToFit, minKnotSpacingPPM, fitRangePPM)
%% [fitParamsStep1] = fit_Osprey_PrelimReduced(dataToFit, minKnotSpacingPPM, fitRangePPM)
%   Performs the first step of the LCModel preliminary analysis
%   analogous to the LCModel algorithm. The algorithm is described in:
%       S.W. Provencher, "Estimation of metabolite concentrations from
%       localized in vivo NMR spectra", Magn Reson Med 30(6):672-679 (1993)
%
%   During the first step, the input spectrum is fit using a reduced basis
%   set (Cr, Glu, Ins, GPC, NAA) and simplified model, using a common
%   frequency shift and common Gaussian and Lorentzian linebroadening
%   for all basis functions of the reduced basis set.
%
%   In addition, there is an unregularized baseline contribution. It seems
%   like this initial baseline has a fixed number of splines (18?), but I 
%   need to investigate this further. For now, set the minKnotSpacingPPM
%   parameter (= DKNTMN in LCModel) to 0.16.
%
%   Input:
%       dataToFit       = FID-A data structure
%       basisSet        = FID-A basis set container
%       minKnotSpacing  = Scalar: minimum baseline knot spacing 
%                         (this is the DKNTMN parameter in LCModel)
%       fitRangePPM     = 2-element vector: fit range [ppm]
%                         (this is the range over which the difference 
%                         between spectrum and model is minimized)
%   Output:
%       fitParamsStep1  = Fit parameters:
%                         - amplitudes of basis functions
%                         - zero-order phase       [deg]
%                         - first-order phase      [deg/ppm]
%                         - Gaussian LB            [Hz]
%                         - Lorentzian LB          [Hz]
%                         - global frequency shift [Hz]
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2020-01-14)
%       goeltzs1@jhmi.edu
%
%   History:
%       2020-01-14: First version of the code.
%


%%% 0. SET EXPECTATION VALUES FOR LORENTZIAN LINEBROADENING AND FREQ SHIFTS
% See LCModel manual, Section 11.14 (p. 149)
exT2 = 2.0;     % [Hz] - expectation value for Lorentzian linebroadening
SDT2 = 0.4;     % [Hz] - standard deviation for Lorentzian linebroadening
SDSH = 0.004;   % [Hz] - standard deviation for frequency shifts
scalingT2 = sqrt(dataToFit.txfrq*1e-6 / 85.15); % scaling factor to account for T2 decrease with field strength


%%% 1. CREATE REDUCED BASIS SET %%%
% Create the spline basis functions for the given resolution, fit range,
% and knot spacing parameter.
[splineArray]       = fit_makeSplineBasis(dataToFit, fitRangePPM, minKnotSpacingPPM);
nSplines            = size(splineArray,2);


%%% 2. SET AND GET STARTING VALUES %%%
% Amplitudes for each spline basis
% function:
ampl        = [zeros(nSplines,1)];

% Concatenate together into one large starting value vector.
x0          = [ampl];


%%% 3. CALL NON-LINEAR SOLVER %%%
% Run the non-linear solver to optimize the non-linear parameters. At each
% iteration of the non-linear least squares solver, the linear parameters
% are calculated separately (similar to the VARiable PROjection algorithm).
%
% We're using a Levenberg-Marquardt implementation for the non-linear 
% problem that allows us to impose hard box constraints on the non-linear
% parameters, keeping them within reasonable limits.
% (c) Alexander Drentler
% https://www.mathworks.com/matlabcentral/fileexchange/53449-levenberg-marquardt-toolbox

% Pack everything up into structs to pass it on to the solver.
% ... data:
inputData.dataToFit     = dataToFit;
inputData.splineArray   = splineArray;
% ... settings:
inputSettings.fitRangePPM   = fitRangePPM;

% Set the hard box constraints for the parameters
lb_beta_j       = -Inf*ones(nSplines,1);
ub_beta_j       = +Inf*ones(nSplines,1); % Baseline spline amplitudes

% Concatenate together into LB/UB vectors
lb  = [lb_beta_j];
ub  = [ub_beta_j];

% Set up and run the non-linear solver.
opts.Display    = 'off';
opts.TolFun     = 1e-6;
opts.TolX       = 1e-6;
opts.MaxIter    = 400;
[x] = LevenbergMarquardt(@(x) fit_Osprey_PrelimReduced_Model(x, inputData, inputSettings),x0,lb,ub,opts);


%%% 4. PERFORM FINAL COMPUTATION OF LINEAR PARAMETERS %%%
% After the non-linear optimization is finished, we need to perform the
% final evaluation of the linear parameters (i.e. the amplitudes and
% baseline parameters).
[fitParamsFinal] = fit_Osprey_PrelimReduced_finalLinear(x, inputData, inputSettings);


%%% 5. CREATE OUTPUT %%%
out = dataToFit;
out.specs = complex(fitParamsFinal.spline(1:size(fitParamsFinal.spline,1)/2),fitParamsFinal.spline(size(fitParamsFinal.spline,1)/2+1:end));
out.fids = ifft(fftshift(out.specs,1),[],1);
% Return the fit parameters from the final linear computation to be used in
% the next LCModel analysis step (i.e. the fit with the full basis set):
fitParamsStep3  = fitParamsFinal;
%dummy=fitParamsStep1.ampl;
%fitParamsStep1.ampl=zeros([])
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% MODEL FUNCTION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = fit_Osprey_PrelimReduced_Model(x, inputData, inputSettings)
%   This function receives the current set of non-linear parameters during
%   every iteration of the non-linear least-squares optimization. The
%   parameters are applied to the basis set.
%
%   Then, a limited-memory Broyden-Fletcher-Goldfarb-Shanno solver allowing
%   for boxed constraints is applied to determine the optimal amplitudes of 
%   the linear parameters to the current set of non-linear parameters.
%
%   The linear parameters are usually the amplitudes for the basis
%   functions and cubic baseline splines.
%
%   This function returns the difference between the input data and
%   the current optimized model. This difference is used by the NLLS
%   solver.
%
%   USAGE:
%       F = fit_LCModel_PrelimReduced_Model(x, inputData, inputSettings)
%
%   INPUTS:
%       x           = Vector providing the last set of parameters coming
%                       out of the non-linear solver
%       inputData   = Struct containing all necessary data to prepare the
%                       final linear solving step (data, basis set...).
%       inputSettings = Struct containing all necessary settings.
%
%   OUTPUTS:
%       F           = Difference between data and model. This is the
%                       objective to be least-squares-optimized by the 
%                       non-linear solver.


%%% 1. UNPACK THE INPUT %%%
% ... data:
dataToFit     = inputData.dataToFit;
splineArray   = inputData.splineArray;
% ... settings:
fitRangePPM   = inputSettings.fitRangePPM;
% ... fit parameters
nSplines    = size(splineArray,2); % number of spline basis functions
ampl        = x(1:end); % Amplitudes for each basis function


%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%

% Apply phasing to the spline basis functions
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
B = [real(B); imag(B)];


%%% 3. SET UP AND CALL SOLVER FOR LINEAR PARAMETERS %%%
% To calculate the linearly occurring amplitude parameters for the
% metabolite/MM/lipid basis functions and the baseline basis functions, we
% call the linear L-BFGS-B algorithm.
% (c) Stephen Becker
% (https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper)

% Set up the equation system to be solved
% Concatenate the metabolite/MM/lipid basis functions and the baseline basis
% functions 
% Cut out the data over the fit range, and turn complex into real problem
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end),size(splineArray,1));
data        = [real(dataToFit.specs); imag(dataToFit.specs)];
b           = data;

% The function we want to minimize is the sum of squares of the residual
fcn     = @(x) norm( B*x - b)^2;
AtA     = B'*B; Ab = B'*b;
grad    = @(x) 2*( AtA*x - Ab );

% Define bounds. The lower bounds for the metabolite/MM/lipid basis
% functions are zero. All other parameters are supposed to be unbound.
l = [-inf*ones(nSplines,1)];
u = [inf*ones(nSplines,1)];

% Prepare the function wrapper
fun     = @(x)fminunc_wrapper( x, fcn, grad);
% Request very high accuracy:
opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10);
opts.printEvery     = 0;

% Run the algorithm:
% Feed initial guess from the input parameters
opts.x0 = ampl;
[ampl, ~, ~] = lbfgsb(fun, l, u, opts );


%%% 4. CREATE OBJECTIVE FUNCTION
% The objective function to be minimized by the non-linear least squares
% solver is the fit residual:
F = data - (B*ampl);


end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FINAL LINEAR ITERATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fitParamsFinal] = fit_Osprey_PrelimReduced_finalLinear(x, inputData, inputSettings)
%   This function is applied after the final iteration of the non-linear
%   solver has returned the final set of non-linear parameters.
%
%   At this point, the linear solver has to be run one last time to
%   estimate the final set of linear parameters.
%
%   The function returns all model parameters.
%
%   USAGE:
%       fitParamsFinal = fit_LCModel_PrelimReduced_finalLinear(x, inputData, inputSettings)
%
%   INPUTS:
%       x           = Vector providing the last set of parameters coming
%                     out of the non-linear solver
%       inputData   = Struct containing all necessary data to prepare the
%                     final linear solving step (data, basis set...).
%       inputSettings = Struct containing all necessary settings.
%
%   OUTPUTS:
%       fitParamsFinal = Set of final fit parameters


%%% 1. UNPACK THE INPUT %%%
% ... data:
dataToFit     = inputData.dataToFit;
splineArray   = inputData.splineArray;
% ... settings:
fitRangePPM   = inputSettings.fitRangePPM;
% ... fit parameters
nSplines    = size(splineArray,2); % number of spline basis functions
ampl        = x(1:end); % Amplitudes for each basis function


%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
% Apply phasing to the spline basis functions
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
B = [real(B); imag(B)];


%%% 3. SET UP AND CALL SOLVER FOR LINEAR PARAMETERS %%%
% To calculate the linearly occurring amplitude parameters for the
% metabolite/MM/lipid basis functions and the baseline basis functions, we
% call the linear L-BFGS-B algorithm.
% (c) Stephen Becker
% (https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper)

% Cut out the data over the fit range, and turn complex into real problem
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end),size(splineArray,1));
data        = [real(dataToFit.specs); imag(dataToFit.specs)];
b           = data;

% The function we want to minimize is the sum of squares of the residual
fcn     = @(x) norm( B*x - b)^2;
AtA     = B'*B; Ab = B'*b;
grad    = @(x) 2*( AtA*x - Ab );

% Define bounds. The lower bounds for the metabolite/MM/lipid basis
% functions are zero. All other parameters are supposed to be unbound.
l = [-inf*ones(nSplines,1)];
u = [inf*ones(nSplines,1)];

% Prepare the function wrapper
fun     = @(x)fminunc_wrapper( x, fcn, grad);
% Request very high accuracy:
opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10);
opts.printEvery     = 0;

% Run the algorithm:
% Feed initial guess from the input parameters
opts.x0 = ampl;
[ampl, ~, ~] = lbfgsb(fun, l, u, opts );


%%% 4. CREATE OUTPUT %%%
% Return the final fit parameters
fitParamsFinal.beta_j       = ampl(1:end);
fitParamsFinal.spline       = B*ampl;

% % Plot (comment out if not debugging)
% figure(99)
% plot(data); hold;
% plot(B*ampl(size(A,2)+1:end)); plot(data - (B*ampl) + 1.1*max(data));
% title('Preliminary Analysis with reduced basis set');
% hold;


end 

