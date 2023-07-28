function plotFit(obj,step,secDim, plotRange)
%%  plotFit(obj,step,secDim, plotRange)
%   This method generates a plot from the fit object.
%
%   USAGE:
%       obj.plotFit(step,secDim, plotRange)
%
%   INPUTS:
%       step            = step to plot      
%       secDim          = spectrum along indirect dimension to plot % OPTIONS:   - [] default (plot all) 
%                                                                                - n (plot spectrum with index n)
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

    if nargin < 4
        plotRange = obj.Options{step}.optimFreqFitRange;            % Set plot range
        if nargin < 3
            secDim = 1;                                             % Set second dimensions to plot
            if nargin < 2
                step = obj.step;                                    % Set to last step
            end
        end
    end

%%  Get fit data from object 

    ppm = obj.Data.ppm;                                                                 % Get ppm vector
    data = fftshift(fft(obj.Data.fids,[],1),1);                                         % Get data matrix
    fit = obj.Model{step}.fit.fit;                                                      % Get fit matrix
    residual = obj.Model{step}.fit.residual;                                            % Get residual matrix           
    baseline = obj.Model{step}.fit.baseline;                                            % Get baseline matrix
    metabs = obj.Model{step}.fit.metabs;                                                % Get metabolite matrix
    names   = obj.BasisSets.names(logical(obj.BasisSets.includeInFit(obj.step,:)));     % Get names cell with included metabolites
    MMind = cat(2,find(contains(names,'MM')),find(contains(names,'Lip')));              % Find the first index with macromolecules 
    if isempty(MMind)                                                                   % If no MMs are included 
       MMind = size(metabs,2) + 1;                                                      % Set the MMind to nMM + 1
    end

%%  Generate figure

    figure                                                                              % Initialize figure
    if size(data,2) > 1 && isempty(secDim)                                              % 2D fit and no user provided indirect dimension index = generate a tiled plot
        plotMat = round(sqrt(size(data,2)));                                            % Set dimensions of tiledlayout                                            
        tiledlayout(plotMat,plotMat)                                                    % Initialize tiledlayout
        shift = max(max(real(data(ppm>plotRange(1) & ppm<plotRange(2),:))) + ...        % Calculate shift for residual and individual basis functions
                min(real((residual(ppm>plotRange(1) & ppm<plotRange(2),:)))));
        YAxLim = [min(min(real(data(ppm>plotRange(1) & ppm<plotRange(2),:)))) ...       % Calculate the y-axis limits to make them constant across the tiledlayout
                  max(max(real(data(ppm>plotRange(1) & ppm<plotRange(2),:))) + ...
                  max(abs(residual(ppm>plotRange(1) & ppm<plotRange(2),:))))];
        for secDim = 1 : size(data,2)                                                   % Loop over second dimension 
            nexttile                                                                    % Initialize new tile
            hold;                                                                       % Hold plot becuase we want to see all results
            for rr = 1:size(metabs,2)                                                   % Loop over basis functions
                if rr < MMind(1)                                                        % Plot metabolite basis functions
                    plot(ppm, real(metabs(:,rr,secDim)) - shift/3,...
                        'Color', [110/255 136/255 164/255]);
                else                                                                    % Plot macromolecular basis functions
                    plot(ppm, real(metabs(:,rr,secDim)) - shift*2/3,...
                        'Color', [110/255 136/255 164/255]);
                end
            end                                                                         % End loop over basis functions
            plot(ppm, real(data(:,secDim)), ...                                         % Plot data                              
                'Color', [11/255 71/255 111/255]);
            plot(ppm, real(residual(:,secDim)) + shift, ...                             % Plot residual
                'Color', [11/255 71/255 111/255]);
            plot(ppm, real(fit(:,secDim)), ...                                          % Plot fit
                'Color',[254/255 186/255 47/255], 'LineWidth', 0.1);
            plot(ppm, real(baseline(:,secDim)), ...                                     % Plot baseline
                'Color', [11/255 71/255 111/255], 'LineWidth', 0.1);
            hold off;            
            set(gca, 'XDir', 'reverse', 'XLim', plotRange, 'YLim', YAxLim);             % Clean appearance
            xlabel('chemical shift (ppm)');
        end
    else                                                                                % 1D fit or a user provided indirect dimension 
        hold;                                                                           % Hold plot becuase we want to see all results
        shift = max(max(real(data(ppm>plotRange(1) & ppm<plotRange(2),:))) + ...        % Calculate shift for residual and individual basis functions
            min(real((residual(ppm>plotRange(1) & ppm<plotRange(2),:)))));    
        for rr = 1:size(metabs,2)                                                       % Loop over basis functions
            if rr < MMind(1)                                                            % Plot metabolite basis functions
                    plot(ppm, real(metabs(:,rr,secDim)) - shift/3, ...
                        'Color', [110/255 136/255 164/255]);
            else                                                                        % Plot macromolecular basis functions
                plot(ppm, real(metabs(:,rr,secDim)) - shift*2/3, ...
                    'Color', [110/255 136/255 164/255]);
            end
        end                                                                             % End loop over basis functions
        plot(ppm, real(data(:,secDim)), ...                                             % Plot data 
            'Color', [11/255 71/255 111/255]);
        plot(ppm, real(residual(:,secDim)) + shift, ...                                 % Plot residual
            'Color', [11/255 71/255 111/255]);
        plot(ppm, real(fit(:,secDim)), ...                                              % Plot fit
            'Color',[254/255 186/255 47/255], 'LineWidth', 0.1);
        plot(ppm, real(baseline(:,secDim)), ...                                         % Plot baseline
            'Color', [11/255 71/255 111/255], 'LineWidth', 0.1);
        hold off;
        set(gca, 'XDir', 'reverse', 'XLim', plotRange);                                 % Clean appearance
        xlabel('chemical shift (ppm)');
    end
end