% fit_runFitWater.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [fitParams, resBasisSet] = fit_runFitWater(dataToFit, basisSet, fitModel, fitOpts);
% 
% DESCRIPTION:
% This function performs spectral fitting of the water data in dataToFit, using
% the water basis function in basisSet, using a set of options defined in fitOpts.
% 
% Prior to handing over the data and the basis set to the model function,
% the basis set is resampled to match the spectral resolution of the data.
%
% The function 'fitModel' describes the actual fitting process that is
% being performed. The workflow is designed to encourage development and
% easy testing/implementation of new quantitative models.
%
% OUTPUTS:
% fitParamsWater    = Structure containing all necessary fit parameters
%                       (model-dependent)
% resBasisSet       = resampled basis set (to match the spectral data resolution)
%
% INPUTS:
% dataToFit   = FID-A water data structure
% basisSet    = FID-A water basis set container
% fitModel    = Function containing the water fitting model
% fitOpts     = Structure containing fit options

function [fitParamsWater, resBasisSet] = fit_runFitWater(dataToFit, basisSet, fitModel, fitOpts)

% Parse input arguments
if nargin<4
    fitOpts = fit_defaultFitOpts;
end

%%% 1. PREPARE DATA AND BASIS SETS %%%
% Resample basis set to match data resolution and frequency range
resBasisSet             = fit_resampleBasis(dataToFit, basisSet);


%%% 2. SELECT THE MODEL %%%
switch fitModel
    case 'Osprey'
        [fitParamsWater] = fit_waterOsprey(dataToFit, resBasisSet, fitOpts);
    case 'LCModel'
        [fitParamsWater] = fit_waterOsprey(dataToFit, resBasisSet, fitOpts);
end


end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fitOpts = fit_defaultFitOpts

% Determine fitting range (in ppm) for the water spectra
fitOpts.rangeWater         = [2.0 7.4];        % [ppm] Default: [2.0 7.4]

end
