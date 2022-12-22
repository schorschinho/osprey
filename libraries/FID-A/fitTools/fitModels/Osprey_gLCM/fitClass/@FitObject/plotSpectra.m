function plotSpectra(obj, plotRange)
    % default to the provided fit range
    if nargin < 2
        plotRange = obj.Options.optimFreqFitRange;
    end
        
    figure;
    plot(obj.Data.ppm, real(fftshift(fft(obj.Data.fids,[],1),1)), 'k');
    set(gca, 'XDir', 'reverse', 'XLim', plotRange);
    
end