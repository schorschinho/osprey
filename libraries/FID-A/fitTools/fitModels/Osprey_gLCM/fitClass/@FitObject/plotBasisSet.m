function plotBasisSet(obj, step, plotRange)
% default to the provided fit range
if nargin < 3
    plotRange = obj.Options{step}.optimFreqFitRange;
    if nargin < 2
        step = obj.step;
        plotRange = obj.Options{step}.optimFreqFitRange;
    end
end
            
        % fft time-domain basis set
        fdBasisSpecs = real(fftshift(fft(obj.BasisSets.fids,[],1),1));
        % calculate stagger
        stag = median(max(abs(real(fdBasisSpecs))));
        % calculate text position
        ppm = obj.Data.ppm;
        xText = round((max(plotRange) - min(plotRange))*0.1 + min(plotRange));
                    
        figure;
        hold;

        for rr = 1:length(obj.BasisSets.names)
            % plot basis functions not included in the fit in red
            if obj.BasisSets.includeInFit(rr) == 1
                colorLine = 'k';
            else
                colorLine = 'r';
            end
            plot(ppm, fdBasisSpecs(:,rr) - stag*(rr-1), colorLine);
            text(xText, -stag*(rr-1)+0.5*stag, obj.BasisSets.names{rr}, 'Color', colorLine, 'Clipping', 'on');
            
            
        end
        hold off;
        set(gca, 'XDir', 'reverse', 'XLim', plotRange);
        xlabel('chemical shift (ppm');
            
end