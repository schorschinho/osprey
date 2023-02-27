function plotFit(obj,step,secDim, plotRange)

    % default to the provided fit range
    if nargin < 4
        plotRange = obj.Options{step}.optimFreqFitRange;
        if nargin < 3
            secDim = [];
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
    names   = obj.BasisSets.names(logical(obj.BasisSets.includeInFit(obj.step,:)));
    MMind = cat(2,find(contains(names,'MM')),find(contains(names,'Lip')));
    if isempty(MMind)
       MMind = size(metabs,2) + 1;
    end

    figure
    if size(data,2) > 1 && isempty(secDim)
        plotMat = round(sqrt(size(data,2)));
        tiledlayout(plotMat,plotMat)
        shift = max(max(real(data(ppm>plotRange(1) & ppm<plotRange(2),:))) + min(real((residual(ppm>plotRange(1) & ppm<plotRange(2),:)))));
        YAxLim = [min(min(real(data(ppm>plotRange(1) & ppm<plotRange(2),:)))) max(max(real(data(ppm>plotRange(1) & ppm<plotRange(2),:))) + max(abs(residual(ppm>plotRange(1) & ppm<plotRange(2),:))))];
        for secDim = 1 : size(data,2)
            nexttile 
    
            hold;
            
            for rr = 1:size(metabs,2)
                if rr < MMind(1)
                    plot(ppm, real(metabs(:,rr,secDim)) - shift/3,'Color', [110/255 136/255 164/255]);
                else
                    plot(ppm, real(metabs(:,rr,secDim)) - shift*2/3,'Color', [110/255 136/255 164/255]);
                end
            end
            plot(ppm, real(data(:,secDim)),'Color', [11/255 71/255 111/255]);
            plot(ppm, real(residual(:,secDim)) + shift,'Color', [11/255 71/255 111/255]);
            plot(ppm, real(fit(:,secDim)), 'Color',[254/255 186/255 47/255], 'LineWidth', 0.1);
            plot(ppm, real(baseline(:,secDim)),'Color', [11/255 71/255 111/255], 'LineWidth', 0.1);
            hold off;
            
            set(gca, 'XDir', 'reverse', 'XLim', plotRange, 'YLim', YAxLim);
            xlabel('chemical shift (ppm)');
            
        %     legend('data', 'residual', 'metabolites', 'baseline', 'fit');
        end
    else
        secDim = 1;
        hold;
        shift = max(max(real(data(ppm>plotRange(1) & ppm<plotRange(2),:))) + min(real((residual(ppm>plotRange(1) & ppm<plotRange(2),:)))));    
        for rr = 1:size(metabs,2)
            if rr < MMind(1)
                    plot(ppm, real(metabs(:,rr,secDim)) - shift/3,'Color', [110/255 136/255 164/255]);
            else
                plot(ppm, real(metabs(:,rr,secDim)) - shift*2/3,'Color', [110/255 136/255 164/255]);
            end
        end
        
        plot(ppm, real(data(:,secDim)),'Color', [11/255 71/255 111/255]);
        plot(ppm, real(residual(:,secDim)) + shift,'Color', [11/255 71/255 111/255]);
        plot(ppm, real(fit(:,secDim)), 'Color',[254/255 186/255 47/255], 'LineWidth', 0.1);
        plot(ppm, real(baseline(:,secDim)),'Color', [11/255 71/255 111/255], 'LineWidth', 0.1);
        hold off;
        
        set(gca, 'XDir', 'reverse', 'XLim', plotRange);
        xlabel('chemical shift (ppm)');
            
        %     legend('data', 'residual', 'metabolites', 'baseline', 'fit');
    end
end