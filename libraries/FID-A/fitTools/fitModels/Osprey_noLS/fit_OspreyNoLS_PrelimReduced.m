function [fitParamsStep1] = fit_OspreyNoLS_PrelimReduced(dataToFit, basisSet, fitRangePPM,refShift,refFWHM)
%% [fitParamsStep1] = fit_Osprey_PrelimReduced(dataToFit, basisSet, fitRangePPM)
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


%%% 1. CREATE BASIS SET %%%
fitMM = 0;
metabList = fit_createMetabList({'Ala','Asp','Cr','CrCH2' ...
                     ,'GABA','GPC','GSH','Glc','Gln','Glu' ...
                     ,'Ins','Lac','NAA','NAAG','PCh','PCr' ...
                     ,'Scyllo','Tau'});
reducedBasisSet     = fit_selectMetabs(basisSet, metabList, fitMM);
nMets =reducedBasisSet.nMets;
% Create the spline basis functions for the given resolution, fit range,
% and knot spacing parameter.
[splineArray, ~]    = fit_makeSplineBasis(dataToFit, fitRangePPM, 1/15);
nSplines            = size(splineArray,2);

%%% Calculate ED 
target_ed = 7;
D = diff(eye(nSplines), 2);

lambda = calc_lambda_from_ed(splineArray(:,:,1), D, (target_ed*(fitRangePPM(2)-fitRangePPM(1))));



%%% 2. SET AND GET STARTING VALUES %%%
% Set the starting values for all parameters.
ph0  = 0; % zero-order phase correction [deg]
% gaussLB     = 5; % Common Gaussian dampening [Hz]
gaussLB     = refFWHM * dataToFit.txfrq*1e-6*0.95; % Common Gaussian dampening [Hz]
freqShift   = 0; % Common frequency shift [Hz]
% Amplitudes for each metabolite basis function and each spline basis
% function:
ampl        = [zeros(nMets,1); zeros(nSplines,1)];

% Concatenate together into one large starting value vector.
x0          = [ph0; gaussLB; freqShift; ampl];


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
inputData.resBasisSet   = reducedBasisSet;
inputData.splineArray   = splineArray;
inputData.lambda        = lambda;
inputData.D             = D;
% ... settings:
inputSettings.fitRangePPM   = fitRangePPM;
% Get an estimate for the standard deviation of the noise
% This is to calculate a Q factor, but need to normalize this noise value
% here - leave alone for now
if max(dataToFit.ppm > 10)
    noiseRange = (dataToFit.ppm > 9);
else
    noiseRange = (dataToFit.ppm < -2);
end
dataNoise   = real(dataToFit.specs(noiseRange,1));
dataNoise   = detrend(dataNoise);
stdNoise    = std(dataNoise);
inputData.stdNoise      = stdNoise;


% Set the hard box constraints for the parameters
lb_ph0          = -6.28; 
ub_ph0          = +6.28; % Zero order phase shift [deg]
% lb_ph1          = -10; 
% ub_ph1          = +10;  % First order phase shift [deg/ppm]
lb_gaussLB      = 0; 
ub_gaussLB      = 150; % Gaussian dampening [Hz]
lb_freqShift    = -11.95; 
ub_freqShift    = 11.95-refShift;    % Frequency shift [Hz]
lb_amplMets     = zeros(nMets,1);
ub_amplMets     = Inf*ones(nMets,1); % Metabolite amplitudes
lb_beta_j       = -Inf*ones(nSplines,1);
ub_beta_j       = +Inf*ones(nSplines,1); % Baseline spline amplitudes

% Concatenate together into LB/UB vectors
lb  = [lb_ph0; lb_gaussLB; lb_freqShift; lb_amplMets; lb_beta_j];
ub  = [ub_ph0; ub_gaussLB; ub_freqShift; ub_amplMets; ub_beta_j];

% Set up and run the non-linear solver.
opts.Display    = 'off';
opts.TolFun     = 1e-6;
opts.TolX       = 1e-6;
opts.MaxIter    = 400;
opts.Jacobian   = 'off';
opts.ScaleProblem='none';
opts.Broyden_updates=2;  
[x] = LevenbergMarquardt(@(x) fit_Osprey_PrelimReduced_Model(x, inputData, inputSettings),x0,lb,ub,opts);


%%% 4. PERFORM FINAL COMPUTATION OF LINEAR PARAMETERS %%%
% After the non-linear optimization is finished, we need to perform the
% final evaluation of the linear parameters (i.e. the amplitudes and
% baseline parameters).
[fitParamsFinal] = fit_Osprey_PrelimReduced_finalLinear(x, inputData, inputSettings);


%%% 5. CREATE OUTPUT %%%
% Return the fit parameters from the final linear computation to be used in
% the next LCModel analysis step (i.e. the fit with the full basis set):
fitParamsStep1  = fitParamsFinal;


end

function [ED] = calc_ed_from_lambda(spline_basis, deriv_mat, lambda)

%   inv_mat = inv((spline_basis).' * spline_basis + lambda * (deriv_mat).' * deriv_mat);
%   H       = inv_mat * (spline_basis.' * spline_basis); 
  inv_mat = (spline_basis).' * spline_basis + lambda * (deriv_mat).' * deriv_mat;
  H       = inv_mat \ (spline_basis.' * spline_basis); 
  ED = trace(H);
end


function [F] = ed_obj_fn(x,spline_basis, deriv_mat, target_ed)

  ed = calc_ed_from_lambda(spline_basis, deriv_mat, x);   
  F = ((ed - target_ed) ^ 2);
end

function [lambda] = calc_lambda_from_ed(spline_basis, deriv_mat, target_ed)
    x0 = 1;
    options.MaxFunEvals = 400;
    options.MaxIter = 400;
    lambda = fminsearch(@(x)ed_obj_fn(x,spline_basis, deriv_mat,target_ed),x0,options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% MODEL FUNCTION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F] = fit_Osprey_PrelimReduced_Model(x, inputData, inputSettings)
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
resBasisSet   = inputData.resBasisSet;
splineArray   = inputData.splineArray;
stdNoise      = inputData.stdNoise;
lambda        = inputData.lambda;
D             = inputData.D;
% ... settings:
fitRangePPM   = inputSettings.fitRangePPM;
% ... fit parameters
nMets       = resBasisSet.nMets; % number of metabolite basis functions
nSplines    = size(splineArray,2); % number of spline basis functions
ph0   = x(1); % zero-order phase correction [convert from deg to rad]
gaussLB     = x(2); % Gaussian dampening [Hz]
freqShift   = x(3); % Common frequency shift [Hz]
ampl        = x(4:end); % Amplitudes for each basis function


%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
% Run the time-domain operations on the metabolite basis functions
% (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
beta = -(pi*gaussLB/2).^2/log(0.5);
t = resBasisSet.t;
for ii=1:nMets
    resBasisSet.fids(:,ii) = resBasisSet.fids(:,ii) .* exp(-beta.*t.*t)';    
end
resBasisSet.specs = fftshift(fft(resBasisSet.fids,[],1),1);
resBasisSet = op_freqrange(resBasisSet,fitRangePPM(1),fitRangePPM(end),size(splineArray,1));

% Apply phasing to the spline basis functions
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
if strcmp(dataToFit.valued,'real')
    B = [real(B); sqrt(lambda)*D];
else
    B = [real(B); imag(B); sqrt(lambda)*D];
end

%%% 3. SET UP AND CALL SOLVER FOR LINEAR PARAMETERS %%%
% To calculate the linearly occurring amplitude parameters for the
% metabolite/MM/lipid basis functions and the baseline basis functions, we
% call the linear L-BFGS-B algorithm.
% (c) Stephen Becker
% (https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper)

% Set up the equation system to be solved
% Concatenate the metabolite/MM/lipid basis functions and the baseline basis
% functions 
if strcmp(dataToFit.valued,'real')
    A   = [real(resBasisSet.specs); zeros(size(D,1),resBasisSet.sz(2))];
else
    A   = [real(resBasisSet.specs); imag(resBasisSet.specs); zeros(size(D,1),resBasisSet.sz(2))];
end

AB  = [A B];
% Cut out the data over the fit range, and turn complex into real problem
dataToFit=op_addphase(dataToFit,ph0);
dataToFit=op_freqshift(dataToFit,freqShift);
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end),size(splineArray,1));
if strcmp(dataToFit.valued,'real')
    data        = [real(dataToFit.specs); zeros(size(D,1),1)];
else
    data        = [real(dataToFit.specs); imag(dataToFit.specs); zeros(size(D,1),1)];
end

b           = data;

% The function we want to minimize is the sum of squares of the residual
fcn     = @(x) norm( AB*x - b)^2;
AtA     = AB'*AB; Ab = AB'*b;
grad    = @(x) 2*( AtA*x - Ab );

% Define bounds. The lower bounds for the metabolite/MM/lipid basis
% functions are zero. All other parameters are supposed to be unbound.
l = [zeros(nMets,1);    -inf*ones(nSplines,1)];
u = [inf*ones(nMets,1);  inf*ones(nSplines,1)];

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
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
if strcmp(dataToFit.valued,'real')
     B = [real(B)];
     A   = [real(resBasisSet.specs)];
     data        = [real(dataToFit.specs)];
     AB  = [A B];
    F = (data - (AB*ampl))./(stdNoise.^2);
else
    B = [real(B); imag(B);];
     A   = [real(resBasisSet.specs);imag(resBasisSet.specs);];
     data        = [real(dataToFit.specs);imag(dataToFit.specs);];
     AB  = [A B];
    F = (data - (AB*ampl))./(stdNoise.^2);
end



%%% 5. CALCULATE ANALYTIC JACOBIAN 
j = sqrt(-1); % i

%Computation of the Jacobian
t = resBasisSet.t;

% The size is MxN with M (number of estimated parameters) and N (number of
% points).
%Store the indices of the different partial derivatives 
AmplCol = 4 : (4 + nMets-1);
SplineAmplCol = AmplCol(end) + 1 : (AmplCol(end) + size(splineArray,2));

nparams = 3;
% nparams = 3 + length(AmplCol) + length(SplineAmplCol);

% dimensions and number of matrix entries
[npoints,~] = size(resBasisSet.fids); %number of points and number of basis functions
J = zeros(npoints,nparams);


%derivative wrt ph0
% J(:,1) = j * completeFit;
J(:,1) = j * data;

%derivative wrt gaussLB
for ii = 1:nMets
    J(:,2) = J(:,2) +  fftshift(fft(-resBasisSet.fids(:,ii).*(t.^2)',[],1),1)*ampl(ii);             
end

%derivative wrt freqShift
% for ii=1:nMets
%     J(:,3) = J(:,3) + fftshift(fft(-j*resBasisSet.fids(:,ii).*t',[],1),1)*ampl(ii);          
% end
data_fid = ifft(ifftshift(data,1),[],1);
J(:,3) = fftshift(fft(data_fid*2*j*pi.*t',[],1),1);          

% derivative  wrt basis set  amplitudes 
for ii=1:nMets
    col = AmplCol(ii);
    J(:,col) = fftshift(fft(resBasisSet.fids(:,ii),[],1),1);    
end


% derivative wrt spline amplitudes
% for ii=1:length(SplineAmplCol)
%     col = SplineAmplCol(ii);
%     J(:,col) = J(:,col)+B(:,ii);
% end

if strcmp(dataToFit.valued,'real')
    J = [real(J);];
else    
    J = [real(J); imag(J);];
end
J = J./(stdNoise.^2);
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
resBasisSet   = inputData.resBasisSet;
splineArray   = inputData.splineArray;
stdNoise      = inputData.stdNoise;
lambda        = inputData.lambda;
D             = inputData.D;
% ... settings:
fitRangePPM   = inputSettings.fitRangePPM;

% ... fit parameters
nMets       = resBasisSet.nMets; % number of metabolite basis functions
nSplines    = size(splineArray,2); % number of spline basis functions
ph0         = x(1); % zero-order phase correction [convert from deg to rad]
gaussLB     = x(2); % Gaussian dampening [Hz]
freqShift   = x(3); % Common frequency shift [Hz]
ampl        = x(4:end); % Amplitudes for each basis function


%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
% Run the time-domain operations on the metabolite basis functions
% (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
beta = -(pi*gaussLB/2).^2/log(0.5);

t = resBasisSet.t;
for ii=1:nMets
    resBasisSet.fids(:,ii) = resBasisSet.fids(:,ii) .* exp(-beta.*t.*t)';     
end
resBasisSet.specs = fftshift(fft(resBasisSet.fids,[],1),1);
resBasisSet = op_freqrange(resBasisSet,fitRangePPM(1),fitRangePPM(end),size(splineArray,1));

% Apply phasing to the spline basis functions
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
if strcmp(dataToFit.valued,'real')
    B = [real(B); sqrt(lambda)*D];
else
    B = [real(B); imag(B); sqrt(lambda)*D];
end



%%% 3. SET UP AND CALL SOLVER FOR LINEAR PARAMETERS %%%
% To calculate the linearly occurring amplitude parameters for the
% metabolite/MM/lipid basis functions and the baseline basis functions, we
% call the linear L-BFGS-B algorithm.
% (c) Stephen Becker
% (https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper)

% Set up the equation system to be solved
% Concatenate the metabolite/MM/lipid basis functions and the baseline basis
% functions
if strcmp(dataToFit.valued,'real')
    A   = [real(resBasisSet.specs); zeros(size(D,1),resBasisSet.sz(2))];
else
    A   = [real(resBasisSet.specs); imag(resBasisSet.specs); zeros(size(D,1),resBasisSet.sz(2))];
end
AB  = [A B];
% Cut out the data over the fit range, and turn complex into real problem
dataToFit=op_addphase(dataToFit,ph0);
dataToFit=op_freqshift(dataToFit,freqShift);
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end),size(splineArray,1));
if strcmp(dataToFit.valued,'real')
    data        = [real(dataToFit.specs); zeros(size(D,1),1)];
else
    data        = [real(dataToFit.specs); imag(dataToFit.specs); zeros(size(D,1),1)];
end
b           = data;

% The function we want to minimize is the sum of squares of the residual
fcn     = @(x) norm( AB*x - b)^2;
AtA     = AB'*AB; Ab = AB'*b;
grad    = @(x) 2*( AtA*x - Ab );

% Define bounds. The lower bounds for the metabolite/MM/lipid basis
% functions are zero. All other parameters are supposed to be unbound.
l = [zeros(nMets,1);    -inf*ones(nSplines,1)];
u = [inf*ones(nMets,1);  inf*ones(nSplines,1)];

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
fitParamsFinal.ampl         = ampl(1:nMets);
fitParamsFinal.ph0          = x(1);
fitParamsFinal.gaussLB      = x(2);
fitParamsFinal.freqShift    = x(3);
fitParamsFinal.beta_j       = ampl(nMets+1:end);
fitParamsFinal.ph1       = 0;

% Plot (comment out if not debugging)
% figure
% plot(data); hold;
% plot(AB*ampl);
% plot(B*ampl(size(A,2)+1:end)); plot(data - (AB*ampl) + 1.1*max(data));
% for rr = 1:nMets
%     plot(ampl(rr)*A(:,rr));
% end
% title('Preliminary Analysis with reduced basis set');
% hold;


end 

