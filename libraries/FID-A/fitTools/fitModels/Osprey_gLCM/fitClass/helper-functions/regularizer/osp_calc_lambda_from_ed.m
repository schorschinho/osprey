function lambda = osp_calc_lambda_from_ed(spline_basis, deriv_mat, target_ed, ...
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
    ed = osp_calc_ed_from_lambda_stable(spline_basis, deriv_mat, par(1));
    output = (ed - target_ed) ^ 2;
end
