function [fitParamsStep2] = fit_OspreyNoLS_PrelimKnotSpacing(dataToFit, basisSet, fitRangePPM, fitParamsStep1)
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
% See LCModel manual, Section 11.14 (p. 149)
bl_n = 50;

%%% 1. CREATE BASIS SET %%%

fitMM = 0;
% metabList = fit_createMetabList({'Cr','CrCH2','GPC','GSH','Gln','Glu','Ins','NAA','PCh','PCr'});
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

%%% Calculate ED start and end values for the ED optimization
max_ed = 7;
D = diff(eye(nSplines), 2);


% lambda_space = logspace(-8,10,10000);
% ED_lookup = zeros(10000,1);
% for j = 1 : length(lambda_space)
%     [ED_lookup(j)] = calc_ed_from_lambda(splineArray(:,:,1), D, lambda_space(j));
% end
% 
% [~,i] = min(abs(ED_lookup - (max_ed*(fitRangePPM(2)-fitRangePPM(1)))));
% lambda_start = lambda_space(i);
% [~,i] = min(abs(ED_lookup - 0));
% lambda_end = lambda_space(i);

lambda_start = calc_lambda_from_ed(splineArray(:,:,1), D, (max_ed*(fitRangePPM(2)-fitRangePPM(1))));
lambda_end = calc_lambda_from_ed(splineArray(:,:,1), D, 0);


lambda_vec  = logspace(log10(lambda_start), log10(lambda_end),bl_n);
optim_model_vec = zeros(bl_n,1);
ed_vec          = zeros(bl_n,1);
residual_vec = zeros(bl_n,1);

inputData.dataToFit     = dataToFit;
inputData.splineArray   = splineArray;
inputData.resBasisSet   = reducedBasisSet;
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

inputData.ph0 = fitParamsStep1.ph0; % zero-order phase correction [convert from deg to rad]
inputData.gaussLB = fitParamsStep1.gaussLB; % Gaussian dampening [Hz]
inputData.freqShift = fitParamsStep1.freqShift; % Common frequency shift [Hz]
x0 = horzcat(fitParamsStep1.ampl', zeros(1,nSplines-length(fitParamsStep1.beta_j)),fitParamsStep1.beta_j'); % Amplitudes for each basis function
x0 = x0';


for n = 1 : bl_n
    ed_vec(n) = calc_ed_from_lambda(splineArray(:,:,1), D,lambda_vec(n));
    inputData.lambda        = lambda_vec(n);
    inputData.D             = D;
    x =x0;
    % Calculate residual 
    [fit] = fit_Osprey_PrelimReduced_finalLinear(x, inputData, inputSettings);
    residual_vec(n) = sum((fit.residual(1: size(splineArray,1))).^2);   
%     residual_vec(n) = sum((fit.residual(1: 2* size(splineArray,1))).^2);   
    fits{n} = fit;    
end
% Calcualte mAIC
optim_model_vec = log(residual_vec) + 30*2 * ed_vec / size(splineArray,1); % This one is supposed to be used accroding to Martin's implementation
% optim_model_vec = log(residual_vec)+ 1.5*ed_vec/max(ed_vec);  % This one seems to work better
[~,min_lambda] = min(optim_model_vec);


%%% 5. CREATE OUTPUT %%%
% Return the fit parameters from the final linear computation to be used in
% the next LCModel analysis step (i.e. the fit with the full basis set):
fitParamsStep2.EDpPPM  = ed_vec(min_lambda)/(fitRangePPM(2)-fitRangePPM(1));
fitParamsStep2.lambda  = lambda_vec(min_lambda);
fitParamsStep2.optim_model_vec  = optim_model_vec;
fitParamsStep2.ampl  = fits{n}.ampl;
fitParamsStep2.idx  = min_lambda;
% fitParamsStep2.optim_model_vec_paper  = optim_model_vec_paper;
fitParamsStep2.fits  = fits;


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
%%%%%%%%%%%% FINAL LINEAR ITERATION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fit] = fit_Osprey_PrelimReduced_finalLinear(x, inputData, inputSettings)
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
lambda        = inputData.lambda;
D             = inputData.D;
% ... settings:
fitRangePPM   = inputSettings.fitRangePPM;

% ... fit parameters
nMets       = resBasisSet.nMets; % number of metabolite basis functions
nSplines    = size(splineArray,2); % number of spline basis functions
ph0 = inputData.ph0; % zero-order phase correction [convert from deg to rad]
gaussLB = inputData.gaussLB; % Gaussian dampening [Hz]
freqShift = inputData.freqShift; % Common frequency shift [Hz]
ampl        = x(1:end); % Amplitudes for each basis function


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
fit.residual =data - (AB*ampl);
fit.data = data;
fit.fit = AB*ampl;
fit.baseline = B*ampl(size(A,2)+1:end);
fit.ampl = ampl;


end 