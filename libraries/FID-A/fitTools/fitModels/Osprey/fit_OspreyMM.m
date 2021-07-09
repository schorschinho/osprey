% fit_OspreyMM.m
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
% This model is attempting to emulate the rather elusive original LCModel
% algorithm (Provencher, Magn Reson Med 30:672-679 (1993)).
% 
% OUTPUTS:
% fitParams   = Structure containing all necessary fit parameters
%                   - ampl: amplitudes for metabolite/MM/lipids
%                   - beta_j: amplitudes for baseline spline functions
%                   - lineshape: coefficients for lineshape convolution
%                   - ph0: zero-order phase correction
%                   - ph1: first-order phase correction
%                   - gaussLB: Gaussian dampening (common to all basis
%                       functions)
%                   - lorentzLB: Lorentzian dampening (individual to each
%                       basis function)
%                   - freqShift: frequency shift (individual to each basis
%                       function)
%
% INPUTS:
% dataToFit   = FID-A data structure
% basisSet    = FID-A basis set container
% fitOpts     = Structure containing fit options

function [fitParams, resBasisSet] = fit_OspreyMM(dataToFit, basisSet, fitOpts)

%%% 0. PREPARE DATA AND BASIS SET %%%
dataToFit               = op_zeropad(dataToFit, 2);
% Resample basis set to match data resolution and frequency range
resBasisSet             = fit_resampleBasis(dataToFit, basisSet);


%%% 1. EXTRACT OPTIONS AND PREPARE FIT %%%
% Extract ppm fit range
fitRangePPM             = fitOpts.range;
% Initialize the baseline spline parameters
minKnotSpacingPPM       = fitOpts.bLineKnotSpace; % this is the DKNTMN parameter in LCModel


%%% 2. INITIAL REFERENCING %%%
% Determine initial coarse frequency shift from cross-correlation with
% landmark delta functions for NAA, Cr, Cho.
% if ~isfield(dataToFit,'refShift') %NO referenceing so far
    disp('Running initial referencing...');
    [refShift, ~] = fit_OspreyReferencingMM(dataToFit);
    refFWHM = dataToFit.refFWHM;
    % Apply initial referencing shift
    dataToFitRef = op_freqshift(dataToFit, -refShift);
% else %Referencing was performed on another Subspec
%     disp('Initial was performed on another Subspec...');
%     refShift = dataToFit.refShift;
%     refFWHM = dataToFit.refFWHM;
%     % Apply initial referencing shift
%     dataToFitRef = op_freqshift(dataToFit, -refShift);
% end

%%% 3. PRELIMINARY ANALYSIS STEP 1 %%%
% In step 1 of the preliminary analysis, a reduced basis set (Cr, Glu, Ins,
% GPC, NAA) is used to obtain a first estimate for the frequency shift
% and the zero- and first-order phase shift parameters.
% Well-phased spectra should be fine with the starting values for the phase
% corrections, which are zero degrees (zero order) and zero degrees/ppm
% (first order).
%
% The preliminary analysis performed by LCModel allows to cycle those
% starting values, i.e. do the preliminary analysis with pairs of e.g. [30,
% 30] etc. - this is governed by the choice of the DEGZER/SDDEGZ and
% DEGPPM/SDDEGP parameters.
%
% In that case, the preliminary analysis is run multiple times, and the
% best phasing / referencing shift parameters are chosen.
% (Worth exploring in future versions by passing the starting values for
% the phase corrections as arguments.)
disp('Running preliminary analysis with reduced basis set...');
[fitParamsStep1] = fit_Osprey_PrelimReducedMM(dataToFitRef, resBasisSet, minKnotSpacingPPM, fitRangePPM);


%%% 4. FINAL PRELIMINARY ANALYSIS STEP 2 %%%
% In the final step of the preliminary analysis, the full basis set is used
% with the full LCModel (except for baseline regularization) to obtain
% the final optimal starting values.
disp('Running final preliminary analysis step with full basis set...');
[fitParamsStep2] = fit_OspreyPrelimStep2MM(dataToFitRef, resBasisSet, minKnotSpacingPPM, fitRangePPM, fitParamsStep1, refFWHM);

RFWHM               = 2.5; % choose slightly larger than the LCM default
convolutionRange    = RFWHM * refFWHM/2;
PPMINC              = abs(dataToFit.ppm(1) - dataToFit.ppm(2));
% Calculate number of points
N_s                 = round(convolutionRange/PPMINC);
if mod(N_s,2) == 0
    N_s = N_s - 1; % make odd number
end
if N_s < 0
    N_s = - N_s; % make odd number
end

lineShape = zeros(N_s,1);
lineShape(ceil(N_s/2)) = 1;

% Normalize
lineShape = lineShape/sum(lineShape);


%%% 5. CREATE OUTPUT %%%
% Return fit parameters
fitParams = fitParamsStep2;
fitParams.refShift = refShift;
fitParams.refFWHM = refFWHM;
%fitParams.lineShape = lineShape;
fitParams.prelimParams = fitParamsStep1;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
