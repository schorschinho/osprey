% fit_runFit.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [fitParams, resBasisSet] = fit_runFit(dataToFit, basisSet, fitModel, fitOpts);
% 
% DESCRIPTION:
% This function performs spectral fitting of the data in dataToFit, using
% the basis set in basisSet, using a set of options defined in fitOpts.
%
% Prior to handing over the data and the basis set to the model function,
% the basis set is resampled to match the spectral resolution of the data.
%
% The function 'fitModel' describes the actual fitting process that is
% being performed. The workflow is designed to encourage development and
% easy testing/implementation of new quantitative models.
% 
% OUTPUTS:
% fitParams   = Structure containing all necessary fit parameters
%               (model-dependent)
% resBasisSet = resampled basis set (to match the spectral data resolution)
%
% INPUTS:
% dataToFit   = FID-A data structure
% basisSet    = FID-A basis set container
% fitModel    = Function containing the fitting model
% fitOpts     = Structure containing fit options

function [fitParams, resBasisSet] = fit_runFit(dataToFit, basisSet, fitModel, fitOpts)

% Parse input arguments
if nargin<4
    fitOpts = fit_defaultFitOpts(fitModel);
end

%%% 1. SELECT THE MODEL %%%
switch fitModel
    case 'Osprey'
        [fitParams, resBasisSet] = fit_Osprey(dataToFit, basisSet, fitOpts);
end

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fitOpts = fit_defaultFitOpts(fitModel)

% Depending on the model chosen, set some default values here

switch fitModel
    case 'Osprey'
        % Determine fitting range (in ppm) for the metabolite and water spectra
        fitOpts.range              = [0.2 4.2];        % [ppm] Default: [0.2 4.2]
        fitOpts.rangeWater         = [2.0 7.4];        % [ppm] Default: [2.0 7.4]
        
        % Determine the baseline knot spacing (in ppm) for the metabolite spectra
        fitOpts.bLineKnotSpace     = 0.4;              % [ppm] Default: 0.4.
end

end
