function [ModelOutput] = fit_LCModelParamsToModel(fitParams)
%   This function returns the complete fit, the baseline, the residual,
%   individual contributions from each basis function, the data that are
%   being fit, and the frequency axis.
%
%   USAGE:
%       [ModelOutput] = fit_LCModelParamsToModel(fitParams)
%
%   INPUTS:
%       fitParams    = Struct containing the set of LCModel model parameters
%
%   OUTPUTS:
%       ModelOutput = Struct containing:
%                       - Complete fit
%                       - Baseline
%                       - Residual
%                       - Individual basis function contributions
%                       - Original data that are being fit
%                       - ppm axis of the fit range


%%% 1.CREATE OUTPUT %%%
% Return a struct with all the output parameters
ModelOutput.ppm             = fitParams.ppm;
ModelOutput.data            = fitParams.data;
ModelOutput.completeFit     = fitParams.completeFit;
ModelOutput.baseline        = fitParams.baseline;
ModelOutput.residual        = fitParams.residual;
ModelOutput.indivMets       = fitParams.indivMets;


end
