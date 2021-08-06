function [fitParamsStep2] = fit_OspreyNoLSPrelimStep2(dataToFit, resBasisSet, minKnotSpacingPPM, fitRangePPM, fitParamsStep1, fitParamsStep2)
%% [fitParamsStep2] = fit_OspreyPrelimStep2(dataToFit, resBasisSet, minKnotSpacingPPM, fitRangePPM, fitParamsStep1, refFWHM)
%   Performs the second step of the LCModel preliminary analysis
%   according to the LCModel algorithm. The algorithm is described in:
%       S.W. Provencher, "Estimation of metabolite concentrations from
%       localized in vivo NMR spectra", Magn Reson Med 30(6):672-679 (1993)
%
%   During the second step, the input spectrum is fit using the full basis
%   set allowing individual frequency shifts and Lorentzian dampening
%   for each metabolite basis function. Additionally, a lineshape
%   convolution is applied to account for deviations from the ideal
%   Voigtian lineshape (= Gaussian/Lorentzian convolution).

%   In addition, pre-defined macromolecule and lipid basis functions are 
%   added, as well as an unregularized baseline. MM/lipids and the baseline
%   are phased with the same phasing parameters as the metabolites, but the
%   lineshape convolution is not applied to them.
%
%   Input:
%       dataToFit       = FID-A data structure
%       basisSet        = FID-A basis set container
%       minKnotSpacing  = Scalar: minimum baseline knot spacing 
%                           (this is the DKNTMN parameter in LCModel)
%       fitRangePPM     = 2-element vector: fit range [ppm]
%                           (this is the range over which the difference 
%                           between spectrum and model is minimized)
%       fitParamsStep1  = Fit parameters from the first preliminary
%                           analysis
%       refFWHM         = Preliminary linewidth estimate (in ppm) that is
%                           used to determine the width of the lineshape
%                           function that the basis functions are
%                           subsequently convolved with
%   Output:
%       fitParamsStep2  = Fit parameters:
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

%%% 1. SET UP EXPECTATION VALUES AND SPLINE BASIS FUNCTIONS %%%
% Create an array of normalized cubic baseline spline basis functions.
[splineArray, ~]    = fit_makeSplineBasis(dataToFit, fitRangePPM, 1/15);
nSplines            = size(splineArray,2);
D = diff(eye(nSplines), 2);

%%% 2. SET AND GET STARTING VALUES %%%
% Set the starting values for the non-linear parameters to be passed on
% to the Levenberg-Marquardt NLLS solving algorithm. Note that the
% amplitude parameters do not need to be initialized, as the
% non-negative linear Lawson-Hanson solver does not require starting
% values to be provided.
nBasisFcts  = resBasisSet.nMets + resBasisSet.nMM; % number of basis functions
ph0         = fitParamsStep1.ph0; % zero-order phase correction [deg]
gaussLB     = fitParamsStep1.gaussLB; % Gaussian dampening [Hz]
freqShift   = fitParamsStep1.freqShift; % Frequency shift [Hz] for each basis function
asym        = 0;
lorentzLB   = 0.001 * ones(nBasisFcts,1); % Lorentzian dampening [Hz] for each basis function
freqShift_mets   = zeros(nBasisFcts,1); % Frequency shift [Hz] for each basis function
ampl        = zeros(nBasisFcts+nSplines,1); % Amplitude parameters for basis functions and baseline


% Concatenate all initial guesses together into one large x0 vector.
x0          = [ph0; gaussLB; freqShift; asym; lorentzLB;  freqShift_mets;];


%%% 3. CLL NON-LINEAR SOLVER %%%
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
inputData.resBasisSet   = resBasisSet;
inputData.splineArray   = splineArray;
inputData.lambda        = fitParamsStep2.lambda;
inputData.D             = D;
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
% ... settings:
inputSettings.fitRangePPM   = fitRangePPM;


% Set the hard box constraints for the parameters
nMets   = resBasisSet.nMets;
nMM     = resBasisSet.nMM;
lb_ph0              = -45; 
ub_ph0              = +45; % Zero order phase shift [deg]
lb_gaussLB          = 1e-10; 
ub_gaussLB          = 1.5 * fitParamsStep1.gaussLB; % Gaussian dampening [Hz]
lb_freqShift   = -11.95; 
ub_freqShift   = +11.95; % Frequency shift [Hz]
lb_asym        = -0.25;
ub_asym        = 0.25;
lb_lorentzLB_mets   =  zeros(nMets, 1); 
ub_lorentzLB_mets   = (2 * ones(nMets, 1)); % Lorentzian dampening [Hz] - Metabolites
lb_lorentzLB_MM     = zeros(nMM, 1); 
ub_lorentzLB_MM     =  (2 * ones(nMM, 1)); % Lorentzian dampening [Hz] - MM/Lipids
lb_freqShift_mets   = -1 * ones(nMets,1); 
ub_freqShift_mets   = +1 * ones(nMets,1); % Frequency shift [Hz] - Metabolites
lb_freqShift_MM     = -1 * ones(nMM,1); 
ub_freqShift_MM     = +1 * ones(nMM,1); % Frequency shift [Hz] - MM/Lipids
lb_ampl             = -Inf * ones(nMets+nMM+size(splineArray,2),1); 
ub_ampl             = +Inf * ones(nMets+nMM+size(splineArray,2),1); % Amplitude for metabolite and spline basis functions


% Concatenate together into LB/UB vectors
lb = [lb_ph0; lb_gaussLB; lb_freqShift; lb_asym; lb_lorentzLB_mets; lb_lorentzLB_MM;  lb_freqShift_mets; lb_freqShift_MM;];
ub = [ub_ph0; ub_gaussLB; ub_freqShift; ub_asym; ub_lorentzLB_mets; ub_lorentzLB_MM;  ub_freqShift_mets; ub_freqShift_MM;];


% Set up and run the non-linear solver.
opts.Display    = 'off';
opts.TolFun     = 1e-20;
opts.TolX       = 1e-20;
opts.MaxIter    = 5000;
opts.Jacobian   = 'on';
opts.ScaleProblem='none';
opts.Broyden_updates=2;  
[x] = LevenbergMarquardt(@(x) fit_Osprey_PrelimStep2_Model(x, inputData, inputSettings), x0, lb, ub, opts);


%%% 4. PERFORM FINAL COMPUTATION OF LINEAR PARAMETERS %%%
% After the non-linear optimization is finished, we need to perform the
% final evaluation of the linear parameters (i.e. the amplitudes and
% baseline parameters).
[fitParamsFinal] = fit_Osprey_PrelimStep2_finalLinear(x, inputData, inputSettings);
% fitParamsFinal.Res = Res;
% fitParamsFinal.J = J;
fitParamsFinal.stdNoise = stdNoise;

%%% 5. CREATE OUTPUT %%%
% Return the fit parameters from the final linear computation to be used in
% the next LCModel analysis step (i.e. the fit with the full basis set):
fitParamsStep2  = fitParamsFinal;


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% MODEL FUNCTION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LCModel preliminary analysis step 2 model
function [F,J] = fit_Osprey_PrelimStep2_Model(x, inputData, inputSettings)
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
%       F = fit_LCModel_PrelimStep2_Model(x, inputData, inputSettings)
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
nMets       = resBasisSet.nMets;
nMM         = resBasisSet.nMM;
nBasisFcts  = nMets + nMM; % number of basis functions
nSplines    = size(splineArray,2); % number of spline basis functions
ph0         = x(1); % zero-order phase correction [convert from deg to rad]
gaussLB     = (x(2)); % Gaussian dampening [Hz^2]
freqShift   = x(3); % Frequency shift [Hz] for each basis function
asym   = x(4); % Frequency shift [Hz] for each basis function
lorentzLB   = (x(5:nBasisFcts+4)); % Lorentzian dampening [Hz] for each basis function
freqShift_mets   = x(nBasisFcts+5:2*nBasisFcts+4); % Frequency shift [Hz] for each basis function
% ampl        = x(2*nBasisFcts+5:3*nBasisFcts+4+nSplines); % Amplitudes





%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
% Run the time-domain operations on the metabolite basis functions
% (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
resBasisSet_h = resBasisSet;
t_h = resBasisSet.t;
t = resBasisSet.t;
% Calculate the Asymmetric pseudo-Voigt lineshape
% Gaussian linshape function
f = [(-resBasisSet.spectralwidth/2)+(resBasisSet.spectralwidth/(4*resBasisSet.sz(1))):resBasisSet.spectralwidth/(2*resBasisSet.sz(1)):(resBasisSet.spectralwidth/2)-(resBasisSet.spectralwidth/(4*resBasisSet.sz(1)))];
asym_fwhm = 2 * gaussLB ./ (1 + exp(-asym .* f));
G  = sqrt(4 .* log(2)./pi./(asym_fwhm.^2)) .*  exp(-4 * log(2) * (f./asym_fwhm).^2);
G = ifft(fftshift(G',1),[],1);
G = G(1:resBasisSet.sz(1));
G = G/real(G(1));

% for small fwhm values the td curve can remain constant when using the ifft
% method above - this is a fix
if (length(unique(real(G)))== 1) && (gaussLB ~= 0)       
    G = exp(-gaussLB.*t.*t);
    G = G/real(G(1));
    G = G';
end

for ii=1:nBasisFcts
    resBasisSet.fids(:,ii) = resBasisSet.fids(:,ii) .* exp(-1i*freqShift_mets(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* G;      
end
resBasisSet.specs = fftshift(fft(resBasisSet.fids,[],1),1);
resBasisSet = op_freqrange(resBasisSet,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));


% Apply phasing to the spline basis functions
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
if strcmp(dataToFit.valued,'real')
    B = [real(B); sqrt(lambda)*D];
else
    B = [real(B); imag(B); sqrt(lambda)*D];
end


%%% 3. SET UP THE LINEAR SOLVER %%%
% To calculate the linearly occurring amplitude parameters for the
% metabolite/MM/lipid basis functions and the baseline basis functions, we
% call the linear L-BFGS-B algorithm.
% (c) Stephen Becker
% (https://www.mathworks.com/matlabcentral/fileexchange/35104-lbfgsb-l-bfgs-b-mex-wrapper)

% Convolve the lineshape with the metabolite basis functions only
% (NOT the macromolecules or lipids or baseline splines).
if strcmp(dataToFit.valued,'real')
    A   = [real(resBasisSet.specs); zeros(size(D,1),resBasisSet.sz(2))];
else
    A   = [real(resBasisSet.specs); imag(resBasisSet.specs); zeros(size(D,1),resBasisSet.sz(2))];
end



% Concatenate the metabolite/MM/lipid basis functions and the baseline basis
% functions 
AB = [A B];
% Cut out the data over the fit range, and use real part only
dataToFit=op_addphase(dataToFit,ph0);
dataToFit=op_freqshift(dataToFit,freqShift);
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end),length(splineArray(:,1,1)));
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
l = [zeros(nMets+nMM,1);    -inf*ones(nSplines,1)];
u = [inf*ones(nMets+nMM,1);  inf*ones(nSplines,1)];

% Prepare the function wrapper
fun     = @(x)fminunc_wrapper( x, fcn, grad);
% Request very high accuracy:
opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10);
opts.printEvery     = 0;

% Run the algorithm:
% Feed initial guess from the input parameters
if exist('ampl.mat')
    load('ampl.mat');    
    opts.x0 = ampl;
end
[ampl, ~, ~] = lbfgsb(fun, l, u, opts );
save('ampl.mat','ampl');

%%% 4. ADD SOFT CONSTRAINTS ON AMPLITUDES %%%
% % To impose soft constraints on the amplitudes, we can augment the problem
% % with additional rows in the equation system. This is done in the function
% % fit_createSoftConstrOsprey.
% % (see Wilson et al., MRM 2011)
% [augmA, augmb] = fit_createSoftConstrOsprey(resBasisSet, AB, b, ampl);
% A_augB  = [AB; augmA];
% b_aug   = [b; augmb];
% 
% % Now, run the L-BFGS-B algorithm again with the augmented equation system
% % The function we want to minimize is the sum of squares of the residual
% fcn     = @(x) norm( A_augB*x - b_aug)^2;
% AtA     = A_augB'*A_augB; A_augb = A_augB'*b_aug;
% grad    = @(x) 2*( AtA*x - A_augb );
% % Prepare the function wrapper
% fun     = @(x)fminunc_wrapper(x, fcn, grad);
% % Run the algorithm:
% % Feed initial guess from the input parameters
% opts.x0 = ampl;
% [ampl, ~, ~] = lbfgsb( fun, l, u, opts );


%%% 4. CREATE OBJECTIVE FUNCTION
% The objective function to be minimized by the non-linear least squares
% solver is the fit residual.
%
% The LCModel algorithm wants to minimize not only the sum of squares
% between the data and the model, but imposes regularization on several
% parameters.
% Therefore, we have to reconstruct the fit spectrum from the parameters we
% derived, and then calculate the sum of squares ourselves, in addition to
% the penalty terms.
% For the last step of the preliminary analysis, this is done without any
% regularization at all.

% Return the loss function
% F = 1/sqrt(stdNoise) * SOS + regB + penaltyLBFS; 
% This would be the classic LCModel function to be minimized
% For the preliminary step, just return the functional without any regularization

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

% Model function
A = resBasisSet.specs;
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
AB = [A B];


completeFit = AB*ampl;


%Computation of the Jacobian
t = resBasisSet.t;

% The size is MxN with M (number of estimated parameters) and N (number of
% points).
%Store the indices of the different partial derivatives 

gaussLBCol = 2 * ones(1,nMets + nMM);
freqShiftCol = 3 * ones(1,nMets + nMM);
asymCol = freqShiftCol(end) + 1;
lorentzLBCol = asymCol(end) + 1 : (asymCol(1) + nMets + nMM);
freqShiftColmets = lorentzLBCol(end) + 1 : (lorentzLBCol(end) + nMets + nMM);
AmplCol = freqShiftColmets(end) + 1 : (freqShiftColmets(end) + nMets + nMM);
SplineAmplCol = AmplCol(end) + 1 : (AmplCol(end) + size(splineArray,2));

nparams = length(x);

% dimensions and number of matrix entries
[npoints,~] = size(resBasisSet.fids); %number of points and number of basis functions
J = zeros(npoints,nparams);


%derivative wrt ph0
% J(:,1) = j * completeFit;
J(:,1) = pi/180 * j * data;

%derivative wrt gaussLB
% for ii = 1:nBasisFcts
%     col = gaussLBCol(ii);
%     J(:,col) = J(:,col) +  fftshift(fft(-resBasisSet.fids(:,ii).*(t.^2)',[],1),1)*ampl(ii);             
% end
% Create a numeric solution for this one
resBasisSet_bkp = resBasisSet_h;
if gaussLB ~= 0
    gaussLB_h = gaussLB + eps;
else
    gaussLB_h = eps;
end
asym_fwhm_h = 2 * gaussLB_h ./ (1 + exp(-asym .* f));
G_h = sqrt(4 .* log(2)./pi./(asym_fwhm_h.^2)) .*  exp(-4 * log(2) * (f./asym_fwhm_h).^2);
G_h = ifft(fftshift(G_h',1),[],1);
G_h = G_h(1:resBasisSet_h.sz(1));
G_h = G_h/real(G_h(1));

% for small fwhm values the td curve can remain constant when using the ifft
% method above - this is a fix
if (length(unique(real(G)))== 1) && (gaussLB ~= 0)    
    G_h = exp(-gaussLB.*t_h.*t_h);
    G_h = G_h/real(G_h(1));
    G_h = G_h';
end

for ii=1:nBasisFcts
    resBasisSet_h.fids(:,ii) = resBasisSet_h.fids(:,ii) .* exp(-1i*freqShift_mets(ii).*t_h)' .* exp(-lorentzLB(ii).*t_h)' .* G_h;      
end
resBasisSet_h.specs = fftshift(fft(resBasisSet_h.fids,[],1),1);
resBasisSet_h = op_freqrange(resBasisSet_h,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
Ah = resBasisSet_h.specs;

ABh = [Ah B];
completeFit_h = ABh*ampl;
num_diff = (completeFit_h-completeFit)/gaussLB_h;
col = gaussLBCol(1);
J(:,col) = num_diff;  
resBasisSet_h = resBasisSet_bkp;
%derivative wrt global freqShift 
% for ii=1:nBasisFcts
%     col = freqShiftCol(ii);
%     J(:,col) = J(:,col) + fftshift(fft(-j*resBasisSet.fids(:,ii).*t',[],1),1)*ampl(ii);          
% end
data_fid = ifft(ifftshift(data,1),[],1);
J(:,3) = fftshift(fft(-data_fid*2*j*pi.*t',[],1),1);

%derivative wrt asymmetry
% Create a numeric solution for this one

if asym ~= 0
    asym_h = asym + eps;
else
    asym_h = eps;
end
asym_fwhm_h = 2 * gaussLB ./ (1 + exp(-asym_h .* f));
G_h = sqrt(4 .* log(2)./pi./(asym_fwhm_h.^2)) .*  exp(-4 * log(2) * (f./asym_fwhm_h).^2);
% G_h = sqrt(4 .* log(2)./pi) .*  exp(-4 * log(2) * (f./asym_fwhm_h).^2);
G_h = ifft(fftshift(G_h',1),[],1);
G_h = G_h(1:resBasisSet_h.sz(1));
G_h = G_h/real(G_h(1));

% for small fwhm values the td curve can remain constant when using the ifft
% method above - this is a fix
if (length(unique(real(G)))== 1) && (gaussLB ~= 0)    
    G_h = exp(-gaussLB.*t_h.*t_h);
    G_h = G_h/real(G_h(1));
    G_h = G_h';
end

for ii=1:nBasisFcts
    resBasisSet_h.fids(:,ii) = resBasisSet_h.fids(:,ii) .* exp(-1i*freqShift_mets(ii).*t_h)' .* exp(-lorentzLB(ii).*t_h)' .* G_h;      
end
resBasisSet_h.specs = fftshift(fft(resBasisSet_h.fids,[],1),1);
resBasisSet_h = op_freqrange(resBasisSet_h,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
Ah = resBasisSet_h.specs;

ABh = [Ah B];
completeFit_h = ABh*ampl;
num_diff = (completeFit_h-completeFit)/asym_h;
col = asymCol(1);

J(:,col) = num_diff;          
resBasisSet_h = resBasisSet_bkp;

%derivative wrt lorentzLB
resBasisSet_bkp = resBasisSet_h;
for ii = 1:nBasisFcts
    resBasisSet_h.fids(:,ii) = -resBasisSet_h.fids(:,ii).*t_h';         
end
resBasisSet_h.specs = fftshift(fft(resBasisSet_h.fids,[],1),1);
resBasisSet_h = op_freqrange(resBasisSet_h,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
for ii = 1:nBasisFcts
    col = lorentzLBCol(ii);
    J(:,col) = J(:,col) +  fftshift(fft(resBasisSet_h.fids(:,ii),[],1),1)*ampl(ii);         
end
resBasisSet_h = resBasisSet_bkp;

%derivative wrt freqShift mets
resBasisSet_bkp = resBasisSet_h;
for ii = 1:nBasisFcts
    resBasisSet_h.fids(:,ii) = -j*resBasisSet_h.fids(:,ii).*t_h';         
end
resBasisSet_h.specs = fftshift(fft(resBasisSet_h.fids,[],1),1);
resBasisSet_h = op_freqrange(resBasisSet_h,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
for ii = 1:nBasisFcts
    col = freqShiftColmets(ii);
    J(:,col) = J(:,col) +  fftshift(fft(resBasisSet_h.fids(:,ii),[],1),1)*ampl(ii);         
end

% % derivative  wrt basis set  amplitudes 
% for ii=1:nBasisFcts
%     col = AmplCol(ii);
%     J(:,col) = fftshift(fft(resBasisSet.fids(:,ii),[],1),1);    
% end
% 
% 
% % derivative wrt spline amplitudes
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

function [fitParamsFinal] = fit_Osprey_PrelimStep2_finalLinear(x, inputData, inputSettings)
%   This function is applied after the final iteration of the non-linear
%   solver has returned the final set of non-linear parameters.
%
%   At this point, the linear solver has to be run one last time to
%   estimate the final set of linear parameters.
%
%   The function returns all model parameters.
%
%   USAGE:
%       fitParamsFinal = fit_finalLinearSolver(x, inputData, inputSettings)
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
nMets       = resBasisSet.nMets;
nMM         = resBasisSet.nMM;
nBasisFcts  = nMets + nMM; % number of basis functions
nSplines    = size(splineArray,2); % number of spline basis functions
ph0         = x(1); % zero-order phase correction [convert from deg to rad]
gaussLB     = x(2); % Gaussian dampening [Hz^2]
freqShift   = x(3); % Frequency shift [Hz] for each basis function
asym   = x(4); % Frequency shift [Hz] for each basis function
lorentzLB   = (x(5:nBasisFcts+4)); % Lorentzian dampening [Hz] for each basis function
freqShift_mets   = x(nBasisFcts+5:2*nBasisFcts+4); % Frequency shift [Hz] for each basis function
% ampl        = x(2*nBasisFcts+5:3*nBasisFcts+4+nSplines); % Amplitudes




%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
f = [(-resBasisSet.spectralwidth/2)+(resBasisSet.spectralwidth/(4*resBasisSet.sz(1))):resBasisSet.spectralwidth/(2*resBasisSet.sz(1)):(resBasisSet.spectralwidth/2)-(resBasisSet.spectralwidth/(4*resBasisSet.sz(1)))];
asym_fwhm = 2 * gaussLB ./ (1 + exp(-asym .* f));
G = sqrt(4 .* log(2)./pi./(asym_fwhm.^2)) .*  exp(-4 * log(2) * (f./asym_fwhm).^2);
G = ifft(fftshift(G',1),[],1);
G = G(1:resBasisSet.sz(1));
G = G/real(G(1));

% Run the time-domain operations on the metabolite basis functions
% (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
t = resBasisSet.t;
for ii=1:nBasisFcts
    resBasisSet.fids(:,ii) = resBasisSet.fids(:,ii) .* exp(-1i*freqShift_mets(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* G;      
end
resBasisSet.specs = fftshift(fft(resBasisSet.fids,[],1),1);
resBasisSet = op_freqrange(resBasisSet,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));

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
% Convolve the lineshape with the metabolite basis functions only
% (NOT the macromolecules or lipids or baseline splines).
if strcmp(dataToFit.valued,'real')
    A   = [real(resBasisSet.specs); zeros(size(D,1),resBasisSet.sz(2))];
else
    A   = [real(resBasisSet.specs); imag(resBasisSet.specs); zeros(size(D,1),resBasisSet.sz(2))];
end

% Concatenate the metabolite/MM/lipid basis functions and the baseline basis
% functions 
AB = [A B];
% Cut out the data over the fit range, and use real part only
dataToFit=op_addphase(dataToFit,ph0);
dataToFit=op_freqshift(dataToFit,freqShift);
dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end),length(splineArray(:,1,1)));
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
l = [zeros(nMets+nMM,1);    -inf*ones(nSplines,1)];
u = [inf*ones(nMets+nMM,1);  inf*ones(nSplines,1)];

% Prepare the function wrapper
fun     = @(x)fminunc_wrapper( x, fcn, grad);
% Request very high accuracy:
opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10);
opts.printEvery     = 0;

% Run the algorithm:
% Feed initial guess from the input parameters
if exist('ampl.mat')
    load('ampl.mat');    
    opts.x0 = ampl;
end
[ampl, ~, ~] = lbfgsb(fun, l, u, opts );


%%% 4. ADD SOFT CONSTRAINTS ON AMPLITUDES %%%
% % To impose soft constraints on the amplitudes, we can augment the problem
% % with additional rows in the equation system. This is done in the function
% % fit_createSoftConstrOsprey.
% % (see Wilson et al., MRM 2011)
% [augmA, augmb] = fit_createSoftConstrOsprey(resBasisSet, AB, b, ampl);
% A_augB  = [AB; augmA];
% b_aug   = [b; augmb];
% 
% % Now, run the L-BFGS-B algorithm again with the augmented equation system
% % The function we want to minimize is the sum of squares of the residual
% fcn     = @(x) norm( A_augB*x - b_aug)^2;
% AtA     = A_augB'*A_augB; A_augb = A_augB'*b_aug;
% grad    = @(x) 2*( AtA*x - A_augb );
% % Prepare the function wrapper
% fun     = @(x)fminunc_wrapper(x, fcn, grad);
% % Run the algorithm:
% % Feed initial guess from the input parameters
% opts.x0 = ampl;
% [ampl, ~, ~] = lbfgsb( fun, l, u, opts );


%%% 5. CREATE OUTPUT %%%
% Return the final fit parameters
fitParamsFinal.ampl         = ampl(1:size(A,2));
fitParamsFinal.ph0          = x(1);
fitParamsFinal.gaussLB      = x(2);
fitParamsFinal.freqShift    = x(3);
fitParamsFinal.asym         = x(4);
fitParamsFinal.lorentzLB    = (x(5:nBasisFcts+4)); % Lorentzian dampening [Hz] for each basis function
fitParamsFinal.freqShiftmets   = x(nBasisFcts+5:2*nBasisFcts+4); % Frequency shift [Hz] for each basis function
fitParamsFinal.beta_j       = ampl(size(A,2)+1:end);
fitParamsFinal.ph1       = 0;

% Plot (comment out if not debugging)
% figure
% plot(data); hold;
% plot(AB*ampl);
% plot(B*ampl(size(A,2)+1:end)); plot(data - (AB*ampl) + 1.1*max(data));
% for rr = 1:(nMets+nMM)
%     plot(ampl(rr)*A(:,rr));
% end
% title('Preliminary Analysis with full basis set (unregularized)');
% hold;


end 

