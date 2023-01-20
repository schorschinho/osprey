function plotFit(obj,step,secDim, plotRange)
    % default to the provided fit range
    if nargin < 4
        plotRange = obj.Options{step}.optimFreqFitRange;
        if nargin < 3
            secDim = 1;
            plotRange = obj.Options{step}.optimFreqFitRange;
            if nargin < 2
                step = obj.step;
                secDim = 1;
                plotRange = obj.Options{step}.optimFreqFitRange;
            end
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
    plot(ppm, real(data(:,secDim)), 'k');
    plot(ppm, real(residual(:,secDim)) + max(real(data(:,secDim))), 'k');
    for rr = 1:size(metabs,2)
        plot(ppm, real(metabs(:,rr,secDim) + baseline(:,secDim)), 'g');
    end
    plot(ppm, real(baseline(:,secDim)), 'b');
    plot(ppm, real(fit(:,secDim)), 'r', 'LineWidth', 0.1);
    hold off;
    
    set(gca, 'XDir', 'reverse', 'XLim', plotRange);
    xlabel('chemical shift (ppm');
    
    legend('data', 'residual', 'metabolites', 'baseline', 'fit');
end