function plotFit(obj,step, plotRange)
    % default to the provided fit range
    if nargin < 3
        plotRange = obj.Options{step}.optimFreqFitRange;
        if nargin < 2
            step = obj.step;
            plotRange = obj.Options{step}.optimFreqFitRange;
        end
    end
        
    % calculate plots
    ppm = obj.Data.ppm;
    data = fftshift(fft(obj.Data.fids,[],1),1);
    fit = obj.Model{step}.fit.fit;
    residual = obj.Model{step}.fit.residual;
    baseline = obj.Model{step}.fit.baseline;
    metabs = obj.Model{step}.fit.metabs;
    figure;
    hold;
    plot(ppm, real(data), 'k');
    plot(ppm, real(fit), 'r');
    plot(ppm, real(residual) + max(real(data)), 'k');
    for rr = 1:size(metabs,2)
        plot(ppm, real(metabs(:,rr) + baseline'), 'g');
    end
    plot(ppm, real(baseline'), 'b');
    hold off;
    
    set(gca, 'XDir', 'reverse', 'XLim', plotRange);
    xlabel('chemical shift (ppm');
    
    set(gca, 'XDir', 'reverse', 'XLim', plotRange);
    xlabel('chemical shift (ppm');
end