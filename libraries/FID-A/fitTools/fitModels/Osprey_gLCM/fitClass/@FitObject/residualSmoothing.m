function residualSmoothing(obj,FilterWindow, ExtrapolationPoints, FilterFunction)
%%  residualSmoothing(obj,FilterWindow)
%   This method generates a plot from the fit object.
%
%   USAGE:
%       obj.plotFit(step,secDim, plotRange)
%
%   INPUTS:
%       FilterWindow        = number of points for smoothing filter
%       ExtrapolationPoints = number of points for extrapolation
%       FilterFunction      = Type of filter function
%
%   OUTPUTS:
%       updated fit
%       baseline estimate
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
        FilterFunction = 'Gaussian';
        if nargin < 3
            ExtrapolationPoints = 10;
            if nargin < 2
                FilterWindow = 100;                                    % Gaussian Window
            end
        end
    end


%% Get residual with broad baseline signals
data = fftshift(fft(obj.Data.fids,[],1),1);                         % Get data matrix
fit = obj.Model{end}.fit.fit;                                      % Get fit matrix
residual = data - fit;                                              % Calculate residual

%% Perform the Gaussian smoothing 
% This code is analog/translated from Tarquin  
    wind_fun = zeros(1, 2 * FilterWindow + 1);                             % Initialise window function vector and set elements to equal zero
    
    % Generate Filter window function
    switch FilterFunction
        case 'Gaussian'
            for k = -FilterWindow:FilterWindow 
                wind_fun(k+FilterWindow+1) = exp(-4 * k^2 / FilterWindow^2);                    % Gaussian-shape
            end
        case 'SineBell'
            for k = -FilterWindow:FilterWindow       
                wind_fun(k+FilterWindow+1) = cos(k * pi / (2.0 * FilterWindow + 2.0));            % Sine-bell-shape 
            end
    end
    
    % Normalise window function
    wind_fun_sum = sum(wind_fun);
    wind_fun = wind_fun / wind_fun_sum;
    
    % Create baseline vector
    N = length(residual);
    baseline = zeros(N,1);
    
    % Apply the convolution
    for n = (1 + FilterWindow):(N - FilterWindow)
        prod = 0;
        for m = (n - FilterWindow):(n + FilterWindow)
            prod = prod + wind_fun(m - (n - FilterWindow) + 1) * residual(m);
        end
        baseline(n,1) = prod;
    end
    
    % Go back to the start and extrapolate first GaussianWindow points
    for k = -FilterWindow:-1
        kc = k;
        baseline(1 + FilterWindow + k,1) = baseline(1 + FilterWindow,1) - kc * (baseline(1 + FilterWindow,1) - baseline(1 + FilterWindow + ExtrapolationPoints,1)) / ExtrapolationPoints;
    end
    
    % Go to the end and extrapolate last GaussianWindow points
    for k = 1:FilterWindow
        kc = k;
        baseline(N - FilterWindow + k,1) = baseline(N - FilterWindow,1) + kc * (baseline(N - FilterWindow,1) - baseline(N - FilterWindow - ExtrapolationPoints,1)) / ExtrapolationPoints;
    end
%% Update fit and residual in object
    obj.Model{end}.fit.fit = fit + baseline;             % Calculate fit matrix
    obj.Model{end}.fit.residual = residual - baseline;   % Calculate residual matrix           
    obj.Model{end}.fit.baseline = baseline;              % Store baseline matrix

end