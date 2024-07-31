function plotFit3D(obj,step,secDim, plotRange)
%%  plotFit3D(obj,step,secDim, plotRange)
%   This method generates a 3D plot from the fit object.
%
%   USAGE:
%       obj.plotFit3D(step,secDim, plotRange)
%
%   INPUTS:
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
%       Dr. Dunja Simicic (Johns Hopkins University, 2023-08-04)
%       dsimici1@jh.edu
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
            secDim = 0;                                             % Set second dimensions to plot (if not defined, it will set it o 0 and plot all)
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
    
    if isVariableRange(secDim)                                                          % Check if the input variable is a Range
        startRange= floor(secDim(1));
        endRange=ceil(secDim(2));
        if startRange>=1                                                                % Check that the range has correct boundaries
            if endRange<=size(data,2) 
                 dim_r = 1: size (data,2);
                 dim = repmat(dim_r, [size(data,1) 1]);
                 for ss = startRange : endRange
                       plot3(ppm,dim(:,ss),real(data(:,ss)),'k', ...
                           'Linewidth',0.4,'Color', [11/255 71/255 111/255])            % plot data
                         
                         hold on
                 end

                 hold on
                 dim=dim+0.05;                                                          % for some reason matlab has a problem when you put the fit exactly on top of the data, so this adds a small offset
                 for ss = startRange : endRange
                    
                       plot3(ppm,dim(:,ss),real(fit(:,ss)),'k', ...
                           'Linewidth' ,1.5, 'Color', [254/255 186/255 47/255])         % plot fit
                       
    
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

            plot3(ppm,dim(:,ss),real(data(:,ss)),'k','Linewidth',0.4, ...               % plot data
                'Color', [11/255 71/255 111/255])
            hold on
        end
        hold on
        dim=dim+0.05;                                                      
        for ss = 1:size(fit,2)

            plot3(ppm,dim(:,ss),real(fit(:,ss)),'k','Linewidth',1.5, ...                % plot fit
                'Color', [254/255 186/255 47/255])
    
        end
    
        set(gca, 'XLim', plotRange);                                 % Clean appearance
        xlabel('chemical shift (ppm)');
        view(191,32)
        hold off

    elseif secDim > 0 && secDim <=size(data,2)                                          % if you defined the dimension, it will plot the no. of spectra you defined
        dim_r = 1: secDim;
        dim = repmat(dim_r, [size(data,1) 1]);

        for ss = 1:secDim

            plot3(ppm,dim(:,ss),real(data(:,ss)),'k','Linewidth',0.4, ...               % plot data
                'Color', [11/255 71/255 111/255])
            hold on
        end
        hold on
        dim=dim+0.05;
        for ss = 1:secDim

            plot3(ppm,dim(:,ss),real(fit(:,ss)),'k','Linewidth',1.5, ...                 % plot fit
                'Color', [254/255 186/255 47/255])
    
        end
  
        set(gca, 'XLim', plotRange);                                  % Clean appearance
        %set(gca, 'XDir', 'reverse');
        xlabel('chemical shift (ppm)');
        view(191,32)
        hold off
        end
    

    % check if the variable is a range

end
    function isRange = isVariableRange(variable)
    % checks if the variable is a range
          if isnumeric(variable) && numel(variable)==2 && isvector(variable)
            isRange=true;
          else
            isRange=false;
          end
    end
   