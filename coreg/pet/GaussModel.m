function F = GaussModel(params,x)
% Function for Gaussian distribution model with linear baseline term

% x(1) = gaussian amplitude
% x(2) = sigma (standard deviation)
% x(3) = mu (mean), ie center of distribution
% x(4) = amplitude of linear baseline
% x(5) = constant amplitude offset

F = params(1)/(params(2)*sqrt(2*pi)) * exp(-(1/2)*((x-params(3))/params(2)).^2) + (params(4)*(x-params(3))+params(5));
end