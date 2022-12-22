function plotInitFit(obj, plotRange)
    % default to the provided fit range
    if nargin < 2
        plotRange = obj.Options{1}.optimFreqFitRange;
    end
        
    % calculate plots
    data        = fftshift(fft(obj.Data.fids, [], 1),1);
    ppm         = obj.Data.ppm;
    fit         = obj.Model{1}.fit.fit;
    residual    = obj.Model{1}.fit.residual;
    baseline    = obj.Model{1}.fit.baseline;
    metabs      = obj.Model{1}.fit.metabs;
    
    figure;
    hold;
    plot(ppm, real(data), 'k');
    plot(ppm, real(fit), 'r');
    plot(ppm, real(residual) + max(real(data)), 'k');
    for rr = 1:size(metabs,2)
        plot(ppm, real(metabs(:,rr) + baseline), 'g');
    end
    plot(ppm, real(baseline), 'b');
    hold off;
    
    set(gca, 'XDir', 'reverse', 'XLim', plotRange);
    xlabel('chemical shift (ppm');
    
end