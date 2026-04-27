function plotBaselineBasis(obj, plotRange)
%%  plotBaselineBasis(obj, plotRange)
%   This method generates a plot of the baseline basis 
%
%   USAGE:
%       obj.plotBaselineBasis(plotRange)
%
%   INPUTS:
%       plotRange       = set plot range default is optimFreqFitRange
%
%   OUTPUTS:
%       figure
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

    if nargin < 2
        plotRange = obj.Options{step}.optimFreqFitRange;                    % Set plot range
    end
    
%%  Get fit data from object 

    ppm = obj.Data.ppm;                                                     % Get ppm vector
    splineBasis = obj.BaselineBasis;                                        % Get baseline basis

%%  Generate figure
    figure                                                                  % Initialize figure
    plot(ppm, real(splineBasis), 'k');                                      % Plot baseline basis
    set(gca, 'XDir', 'reverse', 'XLim', plotRange);                         % Clean appearance
    xlabel('chemical shift (ppm');
    
end