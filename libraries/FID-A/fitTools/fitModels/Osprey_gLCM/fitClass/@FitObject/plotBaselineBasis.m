function plotBaselineBasis(obj, step, plotRange)
    % default to the provided fit range
    if nargin < 3
        plotRange = obj.Options{step}.optimFreqFitRange;
        if nargin < 2
            step = obj.step;
            plotRange = obj.Options{step}.optimFreqFitRange;
        end
    end
        
    figure;
    plot(obj.Data.ppm, real(obj.BaselineBasis), 'k');
    set(gca, 'XDir', 'reverse', 'XLim', plotRange);
    xlabel('chemical shift (ppm');
    
end