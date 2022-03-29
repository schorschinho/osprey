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

function [fitParams, resBasisSet] = fit_OspreyMMGaussians(dataToFit, basisSet, fitOpts)

%%% 0. PREPARE DATA AND BASIS SET %%%
dataToFit               = op_zeropad(dataToFit, 2);
% Resample basis set to match data resolution and frequency range
resBasisSet             = fit_resampleBasis(dataToFit, basisSet);


%%% 1. EXTRACT OPTIONS AND PREPARE FIT %%%
% Extract ppm fit range
fitRangePPM             = fitOpts.range;

dataToFitRef = op_freqshift(dataToFit, -fitOpts.Params2.refShift);
%%% 4. FINAL PRELIMINARY ANALYSIS STEP 2 %%%
% In the final step of the preliminary analysis, the full basis set is used
% with the full LCModel (except for baseline regularization) to obtain
% the final optimal starting values.
disp('Running final preliminary analysis step with full basis set...');
[fitParamsStep4] = fit_OspreyPrelimStep2MMGaussians(dataToFitRef, resBasisSet, fitRangePPM, fitOpts.Params2);



%%% 5. CREATE OUTPUT %%%
% Return fit parameters
fitParams = fitParamsStep4;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
