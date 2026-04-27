function plotResidual3D(obj,newFigure,step,secDim, plotRange)
%%  plotResidual3D(obj,step,secDim, plotRange)
%   This method generates a 3D plot of the residual from the fit object.
%
%   USAGE:
%       obj.plotFit3D(step,secDim, plotRange)
%
%   INPUTS:
%       newFigure       = create a new figure 
%       step            = step to plot      
%       secDim          = spectrum along indirect dimension to plot % OPTIONS:   - default (plot all) 
%                                                                                - n (plot the first n spectra)
%                                                                                - plot spectra in a range e.g. [7, 17]
%       plotRange       = set plot range default is optimFreqFitRange
%
%   OUTPUTS:
%       figure
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2025-07-02)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%%  Diverge to default options if required 

    if nargin < 5
        if nargin < 3
            step = obj.step;                                    % Set to last step
        end
        plotRange = obj.Options{step}.optimFreqFitRange;            % Set plot range
        if nargin < 4
            secDim = 0;                                             % Set second dimensions to plot (if not defined, it will set it o 0 and plot all)
            if nargin < 3
                step = obj.step;                                    % Set to last step
                if nargin < 2
                    newFigure = 1;                                  % create new figure
                end
            end
        end
    end

%%  Get fit data from object 

    ppm = obj.Data.ppm;                                                                 % Get ppm vector
    residual = obj.Model{step}.fit.residual;                                            % Get residual matrix           
    data = fftshift(fft(obj.Data.fids,[],1),1);                                         % Get data matrix
%%  Generate figure
    if newFigure
        figure;                                                                     % Initialize figure
    end                                                                            % Initialize figure
    
    if isVariableRange(secDim)                                                          % Check if the input variable is a Range
        startRange= floor(secDim(1));
        endRange=ceil(secDim(2));
        if startRange>=1                                                                % Check that the range has correct boundaries
            if endRange<=size(data,2) 
                 dim_r = 1: size (data,2);
                 dim = repmat(dim_r, [size(data,1) 1]);
                 for ss = startRange : endRange
                       plot3(ppm,dim(:,ss),real(residual(:,ss)),'k', ...
                           'Linewidth',0.4,'Color', [11/255 71/255 111/255])            % plot residual
                         
                         hold on
                 end

                 set(gca, 'XLim', plotRange);                        % Clean appearance
                 xlabel('chemical shift (ppm)');
                 view(191,32)
                 hold off
            end
        end
    
    elseif size(data,2) > 1 && secDim==0                                                % no defined second dimention = plot all
        dim_r = 1: size (data,2);
        dim = repmat(dim_r, [size(data,1) 1]);
    
        for ss = 1:size(data,2)

            plot3(ppm,dim(:,ss),real(residual(:,ss)),'k','Linewidth',0.4, ...               % plot residual
                'Color', [11/255 71/255 111/255])
            hold on;
        end


    elseif secDim > 0 && secDim <=size(data,2)                                          % if you defined the dimension, it will plot the no. of spectra you defined
        dim_r = 1: secDim;
        dim = repmat(dim_r, [size(data,1) 1]);

        for ss = 1:secDim

            plot3(ppm,dim(:,ss),real(residual(:,ss)),'k','Linewidth',0.4, ...               % plot residual
                'Color', [11/255 71/255 111/255])
            hold on;
        end
        
  
        end
        set(gca, 'XLim', plotRange,...
            'LineWidth', 1, 'TickDir', 'out',...
            'YTickLabel',{},'YTick',{},...
            'ZTickLabel',{},'ZTick',{},...
            'XDir','reverse',...
            'XColor', [11/255 71/255 111/255], ...
            'YColor', [11/255 71/255 111/255], ...
            'ZColor', [11/255 71/255 111/255]);                                  % Clean appearance
        xlabel('chemical shift (ppm)');
        % view(191,32)
        hold off

end
    function isRange = isVariableRange(variable)
    % checks if the variable is a range
          if isnumeric(variable) && numel(variable)==2 && isvector(variable)
            isRange=true;
          else
            isRange=false;
          end
    end
   