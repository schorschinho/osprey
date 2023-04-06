function fh = LambdaEDFormalism
    fh.SetLambdaSpace = @set_lambda_space;
    fh.LtoED = @calc_lambda_from_ed;
    fh.EDtoL = @calc_ed_from_lambda;
    fh.optimLambda = @find_optimal_lambda_AIC;
end

function LambdaSpace = set_lambda_space(ppm,optimFreqFitRange,SplineBasis,lambdaRange,RegMatrix,steps)
    [indMin, indMax] = ppmToIndex(ppm, optimFreqFitRange);
    spline_basis =SplineBasis(indMin:indMax,:);
    deriv_mat = RegMatrix.fun(size(spline_basis,2));
    lambda_min = osp_calc_lambda_from_ed(real(spline_basis), real(deriv_mat), lambdaRange.ub*abs(optimFreqFitRange(1)-optimFreqFitRange(2)));
    if isempty(lambdaRange.lb)
        lambda_max = osp_calc_lambda_from_ed(real(spline_basis), real(deriv_mat), 2.01);
    else
        lambda_max = osp_calc_lambda_from_ed(real(spline_basis), real(deriv_mat), lambdaRange.lb*abs(optimFreqFitRange(1)-optimFreqFitRange(2)));
    end
    LambdaSpace = logspace(log10(lambda_min),log10(lambda_max),steps);
end

function [ed] = calc_ed_from_lambda(spline_basis, deriv_mat, lambda)
    inv_mat = pinv([spline_basis; lambda^0.5 * deriv_mat]);
    zero_mat = zeros(size(deriv_mat, 1), size(deriv_mat, 2));
    G = inv_mat * [spline_basis; zero_mat];
    ed = sum(diag(G));
end

function lambda = calc_lambda_from_ed(spline_basis, deriv_mat, target_ed, ...
                                      upper_lim, lower_lim, start_val)
    if nargin < 4
        upper_lim = 1e10;
    end
    if nargin < 5
        lower_lim = 1e-6;
    end
    if nargin < 6
        start_val = 1.0;
    end
    
    options = optimset('TolX', 1e-6);
    [lambda,~,flag] = fminbnd(@(x) ed_obj_fn(x, spline_basis, deriv_mat, target_ed), ...
                               lower_lim, upper_lim, options);
    
    if flag ~= 1
        warning('correct lambda not found');
    end
end

function [output] = ed_obj_fn(par, spline_basis, deriv_mat, target_ed)
    ed = calc_ed_from_lambda(spline_basis, deriv_mat, par(1));
    output = (ed - target_ed) ^ 2;
end

function [optimLambda, optimLambdaIndex,ed,AIC] = find_optimal_lambda_AIC(residual,ppm,optimFreqFitRange,SplineBasis,RegMatrix,steps,m,RegParSpace)
    [indMin, indMax] = ppmToIndex(ppm, optimFreqFitRange);
    spline_basis =SplineBasis(indMin:indMax,:);
    deriv_mat = RegMatrix.fun(size(spline_basis,2));

    ed = zeros(steps,1);
    AIC = zeros(steps,1);
    for oo = 1 : steps
        ed(oo) = osp_calc_ed_from_lambda_stable(real(spline_basis), real(deriv_mat), RegParSpace(oo));
        AIC(oo) = log(sum((real(residual(indMin:indMax,oo)).^2))) + 2 * m * real(ed(oo))/length(residual(indMin:indMax,oo));
    end
    
    [~,AIC_ind] = min(AIC);

    optimLambdaIndex = AIC_ind;
    optimLambda = RegParSpace(AIC_ind);

end