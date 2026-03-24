%%  SecondDerivative
%   This function describes a penalty matrix for parameter regualrization.
%   The second order differencing matrix imposes smoothness between
%   neighbooring parameters.
%   It is specified in the RegFun field in the model procedure json.
%
%   USAGE:
%       Specify this in the parametrizations.(parameter).RegFun field in the model
%       procedure. This is currently used for the p-spline baseline
%       optimization but can in principle be applied as smoother to an
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Handle set up for optimizer
function fh = SecondDerivative
    fh.fun = @ScndDer;
end

%% Generate penalty matrix
function RegMatrix = ScndDer(nPars)
%   This function describes a penalty matrix for parameter regualrization.
%   The second order differencing matrix imposes smoothness between
%   neighbooring parameters. More information to be found in:
%   Eilers et al., Statist. Sci. 11 (2) 89 - 121,(1996)
%
%   USAGE:
%       RegMatrix  = ScndDer(nPars)
%
%   INPUTS:
%       nPars      = number of parameters to regularize

%   OUTPUTS:
%       RegMatrix  = Regularizer matrix
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Calcualte matrix
    e = ones(nPars,1);                                  % Create single column vector 
    RegMatrix = spdiags([e -2*e e],-1:1,nPars,nPars);   % Create second differencing matrix
    RegMatrix = full(RegMatrix);                        % Convert sparse representation to full matrix
end
