function [ed] = osp_calc_ed_from_lambda_stable(spline_basis, deriv_mat, lambda)
    inv_mat = pinv([spline_basis; lambda^0.5 * deriv_mat]);
    zero_mat = zeros(size(deriv_mat, 1), size(deriv_mat, 2));
    G = inv_mat * [spline_basis; zero_mat];
    ed = sum(diag(G));
end
