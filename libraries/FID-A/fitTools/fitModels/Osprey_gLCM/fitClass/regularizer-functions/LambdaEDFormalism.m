%%  LambdaEDFormalism
%   This function contains all functions to perform a effective dimension
%   optimization of p-splines described in Martin Wilson's ABfit algorithm 
%   It is specified in the regularizer.fun field in the model procedure
%   json. For more information about p-Splines read this:
%   Eilers et al., Statist. Sci. 11 (2) 89 - 121,(1996)
%   ABFit algorithm is descirbed here:
%   Wilson Magn Reson Med. (2021):85:13-29 
%   or here
%   https://github.com/martin3141/spant/blob/master/R/abfit.R
%
%   USAGE:
%       Specify this in the regularizer.fun.
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
function fh = LambdaEDFormalism
    fh.SetLambdaSpace = @set_lambda_space;      % Generate lambda space
    fh.LtoED = @calc_lambda_from_ed;            % Calculate lambda from effective dimension
    fh.EDtoL = @calc_ed_from_lambda;            % Calculate effective dimension from lambda
    fh.optimLambda = @find_optimal_lambda_AIC;  % Find optimal lambda using AIC
end

%% Functions for Lambda to effective dimension formalism

function LambdaSpace = set_lambda_space(ppm,optimFreqFitRange,SplineBasis,lambdaRange,RegMatrix,steps)
% This function generates a regularizer paramter space for a set of ppm
% range, spline basis set, lambda range, regularizer matrix, and steps.
% This follows the description in the ABfit algorithm.
% Wilson Magn Reson Med. (2021):85:13-29 
% or here
% https://github.com/martin3141/spant/blob/master/R/abfit.R
%
%   USAGE:
%       LambdaSpace = set_lambda_space(ppm,optimFreqFitRange,SplineBasis,lambdaRange,RegMatrix,steps)
%
%   INPUTS:
%       ppm                 = ppm axis
%       optimFreqFitRange   = model frequency range
%       SplineBasis         = spline basis struct
%       lambdaRange         = lambda range struct
%       RegMatrix           = function handle to regularizer matrix
%       steps               = number of steps in the regularizer space
%
%   OUTPUTS:
%       LambdaSpace        = space of regularizer parameters
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
%% Generate lambda space
    [indMin, indMax] = ppmToIndex(ppm, optimFreqFitRange);                  % Get fit range ppm indices
    spline_basis =SplineBasis(indMin:indMax,:);                             % Cut out fit range from spline basis
    deriv_mat = RegMatrix.fun(size(spline_basis,2));                        % Generate regularizer matrix
    lambda_min = calc_ed_from_lambda(real(spline_basis), ...            % Calculate minimum lambda from effective dimension with given parameters
                                         real(deriv_mat), ...
                                         lambdaRange.ub*abs(optimFreqFitRange(1)-optimFreqFitRange(2)));
    % Calculate maximum lambda from effective dimension     
    if isempty(lambdaRange.lb)                                              % No lower ED given set to 2.01                                                                              
        lambda_max = calc_lambda_from_ed(real(spline_basis), real(deriv_mat), 2.01);    % Calcualte lambda max
    else
        lambda_max = calc_lambda_from_ed(real(spline_basis), ...            % Calculate maximum lambda from effective dimension with given parameters
                                             real(deriv_mat), ...
                                             lambdaRange.lb*abs(optimFreqFitRange(1)-optimFreqFitRange(2)));
    end
    LambdaSpace = logspace(log10(lambda_min),log10(lambda_max),steps);      % Generate lambda space in log spacing
end

function ed = calc_ed_from_lambda(SplineBasis, RegMatrix, lambda)
% This function calculates the effective dimension for given spline basis
% regularizer matrix and lambda value.
% This follows the description in the ABfit algorithm.
% Wilson Magn Reson Med. (2021):85:13-29 
% or here
% https://github.com/martin3141/spant/blob/master/R/abfit.R
%
%   USAGE:
%       ed  = calc_ed_from_lambda(spline_basis, deriv_mat, lambda)
%
%   INPUTS:
%       SplineBasis     = spline basis struct
%       RegMatrix       = function handle to regularizer matrix
%       lambda          = lambda value
%
%   OUTPUTS:
%       ed        = effective dimension
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
%%  Calculate ED
    inv_mat = pinv([SplineBasis; lambda^0.5 * RegMatrix]);          % Invert splinebasis and regularizer matrix
    zero_mat = zeros(size(RegMatrix, 1), size(RegMatrix, 2));       % Generate zero matrix
    H = inv_mat * [SplineBasis; zero_mat];                          % Calculate H   (equation 4 in Wilson 2021)
    ed = sum(diag(H));                                              % Calculate effective dimension (equation 5 in Wilson 2021)
end

function lambda = calc_lambda_from_ed(SplineBasis, RegMatrix, target_ed, ...
                                      upper_lim, lower_lim, start_val)
% This function calculates lambda for a given spline basis, regularizer
% matrix and lambda value.
% This follows the description in the ABfit algorithm.
% Wilson Magn Reson Med. (2021):85:13-29 
% or here
% https://github.com/martin3141/spant/blob/master/R/abfit.R
%
%   USAGE:
%       lambda = calc_lambda_from_ed(spline_basis, deriv_mat, target_ed, upper_lim, lower_lim, start_val)
%                                      
%
%   INPUTS:
%       SplineBasis     = spline basis struct
%       RegMatrix       = function handle to regularizer matrix
%       target_ed       = target effective dimension
%       upper_lim       = upper limit for lambda (optional default 1e10)
%       lower_lim       = lower limit for lambda (optional default 1e-6)
%       start_val       = start value for lambda (optional default 1)
%
%   OUTPUTS:
%       lambda        = final lambda value
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
%%  Diverge to default options if required 
    if nargin < 4
        upper_lim = 1e10;       % Set upper limit lambda value
    end
    if nargin < 5
        lower_lim = 1e-6;       % Set lower limit lambda value
    end
    if nargin < 6
        start_val = 1.0;        % Set start lambda value
    end
%% Calculate lambda    
    options = optimset('TolX', 1e-6);                   % Set up optimizer
    [lambda,~,flag] = fminbnd(@(x) ed_obj_fn(x, SplineBasis, RegMatrix, target_ed), ... % find lambda
                               lower_lim, upper_lim, options);
    
    if flag ~= 1                                % Optimizer failed
        warning('correct lambda not found');
    end
end

function residual = ed_obj_fn(par, SplineBasis, RegMatrix, target_ed)
% This is the object function to get a target effective dimension
% This follows the description in the ABfit algorithm.
% Wilson Magn Reson Med. (2021):85:13-29 
% or here
% https://github.com/martin3141/spant/blob/master/R/abfit.R
%
%   USAGE:
%       residual = ed_obj_fn(par, spline_basis, deriv_mat, target_ed)
%
%   INPUTS:
%       par          = parameter to optimize
%       SplineBasis  = spline basis struct
%       RegMatrix    = function handle to regularizer matrix
%       target_ed    = target effective dimension
%
%   OUTPUTS:
%       residual     = residual to optimize 
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
    ed = calc_ed_from_lambda(SplineBasis, RegMatrix, par(1));   % Calculate ED from parameters
    residual = (ed - target_ed) ^ 2;                            % Return resiudal
end

function [optimLambda, optimLambdaIndex,ed,AIC] = find_optimal_lambda_AIC(residual,ppm,optimFreqFitRange,SplineBasis,RegMatrix,steps,m,RegParSpace)
% This function finds the optimal lambda indes and value, and also returns
% the ED and AIC space for plotting.
% This follows the description in the ABfit algorithm.
% Wilson Magn Reson Med. (2021):85:13-29 
% or here
% https://github.com/martin3141/spant/blob/master/R/abfit.R
%
%   USAGE:
%       [optimLambda, optimLambdaIndex,ed,AIC] = find_optimal_lambda_AIC(residual,ppm,optimFreqFitRange,SplineBasis,RegMatrix,steps,m,RegParSpace)
%
%   INPUTS:
%       residual          = fit residual matrix
%       ppm               = ppm axis
%       optimFreqFitRange = model frequency range
%       SplineBasis       = spline basis struct
%       RegMatrix         = function handle to regularizer matrix
%       steps             = number of steps in the regularizer space
%       m                 = modification factor for AIC
%       RegParSpace       = space of regularizer parameters
%
%   OUTPUTS:
%       optimLambda       = optimal regularization parameter
%       optimLambdaIndex  = index to optimal regularization parameter
%       ed                = effective dimension space
%       AIC               = AIC space
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
%% Calculate optimal lambda

    [indMin, indMax] = ppmToIndex(ppm, optimFreqFitRange);                  % Get fit range ppm indices
    spline_basis =SplineBasis(indMin:indMax,:);                             % Cut out fit range from spline basis
    deriv_mat = RegMatrix.fun(size(spline_basis,2));                        % Generate regularizer matrix

    ed = zeros(steps,1);                                                    % Setup ED space
    AIC = zeros(steps,1);                                                   % Setup AIC space
    for oo = 1 : steps                                                      % Loop over steps
        ed(oo) = calc_ed_from_lambda(real(spline_basis), real(deriv_mat), RegParSpace(oo));                          % Calculate ED
        AIC(oo) = log(sum((real(residual(indMin:indMax,oo)).^2))) + 2 * m * real(ed(oo))/length(residual(indMin:indMax,oo));    % Calculate AICs
    end                                                                     % End loop over steps
    
    [~,AIC_ind] = min(AIC);             % Find minimal AIC for optimal L

    optimLambdaIndex = AIC_ind;         % Return optimal lambda index
    optimLambda = RegParSpace(AIC_ind); % Return optimal lambda value

end