function plotBaselineBasis(obj, plotRange)
    
    % default to the provided fit range
    if nargin < 2
        plotRange = obj.Options.optimFreqFitRange;
    end
        
    figure;
    plot(obj.Data.ppm, real(obj.BaselineBasis), 'k');
    set(gca, 'XDir', 'reverse', 'XLim', plotRange);
    xlabel('chemical shift (ppm');
    
end