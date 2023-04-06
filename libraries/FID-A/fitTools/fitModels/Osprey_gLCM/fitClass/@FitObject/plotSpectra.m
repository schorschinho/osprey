function plotSpectra(obj, secDim, plotRange)
%%  plotSpectra(obj,secDim, plotRange)
%   This method generates a plot from the fit object.
%
%   USAGE:
%       obj.plotFit(step,secDim, plotRange)
%
%   INPUTS:
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

    if nargin < 3
        plotRange = obj.Options{step}.optimFreqFitRange;            % Set plot range
        if nargin < 2
            secDim = 1;                                             % Set second dimensions to plot
        end
    end
%%  Get fit data from object 

    ppm = obj.Data.ppm;                                                     % Get ppm vector
    data = fftshift(fft(obj.Data.fids,[],1),1);                             % Get data matrix    
%%  Generate figure
    figure                                                                  % Initialize figure
    if size(data,2) > 1 && isempty(secDim)                                  % 2D data and no user provided indirect dimension index = generate a tiled plot
        plotMat = round(sqrt(size(data,2)));                                % Set dimensions of tiledlayout                                            
        tiledlayout(plotMat,plotMat)                                        % Initialize tiledlayout
        for secDim = 1 : size(data,2)                                       % Loop over second dimension 
            nexttile                                                        % Initialize new tile
            hold;                                                           % Hold plot becuase we want to see all results
            plot(obj.Data.ppm, real(data(:,secDim)), ... % Plot data
                'Color', [11/255 71/255 111/255]);
            set(gca, 'XDir', 'reverse', 'XLim', plotRange);                 % Clean appearance
        end
    else                                                                    % 1D fit or a user provided indirect dimension 
        secDim = 1;                                                         % Set second dimension index
        hold;                                                               % Hold plot becuase we want to see all results
        plot(obj.Data.ppm, real(data(:,secDim)), ...                        % Plot data
            'Color', [11/255 71/255 111/255]);
        set(gca, 'XDir', 'reverse', 'XLim', plotRange);                     % Clean appearance
    end
    
end