function plotBasisSet(obj, step, plotRange)
%%  plotBasisSet(obj, step, plotRange)
%   This method generates a plot of the basis set 
%
%   USAGE:
%       obj.plotBaselineBasis(step, plotRange)
%
%   INPUTS:
%       step            = step to plot    
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
        plotRange = obj.Options{step}.optimFreqFitRange;                    % Set plot range
        if nargin < 2
            step = obj.step;                                                % Set to last step
        end
    end
  
%%  Get fit data from object 
    if ~isempty(obj.BasisSets.fids)                                         % basis set has been removed to reduce file size
        fdBasisSpecs = real(fftshift(fft(obj.BasisSets.fids,[],1),1));      % fft time-domain basis set
        stag = median(max(abs(real(fdBasisSpecs))));                        % calculate stagger 
    end
    ppm = obj.Data.ppm;                                                     % Get ppm vector
    xText = round((max(plotRange) - min(plotRange))*0.1 + min(plotRange));  % calculate text position
         
    %%  Generate figure
    if ~isempty(obj.BasisSets.fids)                                         % basis set has been removed to reduce file size
        figure;                                                             % Initialize figure
        hold;                                                               % Hold plot becuase we want to see all results  
        for rr = 1:length(obj.BasisSets.names)                              % Loop over basis functions
            % plot basis functions not included in the fit in red
            if obj.BasisSets.includeInFit(rr) == 1
                colorLine = 'k';                                            % Basis function is included in fit
            else
                colorLine = 'r';                                            % Basis function is not included in fit
            end
            plot(ppm, fdBasisSpecs(:,rr) - stag*(rr-1), colorLine);         % Plot basis function
            text(xText, -stag*(rr-1)+0.5*stag, obj.BasisSets.names{rr}, ... % Add basis function name
                'Color', colorLine, 'Clipping', 'on');                
        end                                                                 % End loop over basis functions
        hold off;
        set(gca, 'XDir', 'reverse', 'XLim', plotRange);                     % Clean appearance
        xlabel('chemical shift (ppm');
    else
        fprintf('Basis set has been removed to reduce file size. \n');      % Write warning that basis set has been removed
    end
            
end