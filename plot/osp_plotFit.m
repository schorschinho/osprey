function out = osp_plotFit(MRSCont, kk, which, stagFlag, xlab, ylab, figTitle)
%% out = osp_plotFit(MRSCont, kk, which, stagFlag, xlab, ylab, figTitle)
%   Creates a figure showing data stored in an Osprey data container, as
%   well as the fit to it, the baseline, the residual, and contributions
%   from the individual metabolites.
%
%   USAGE:
%       out = osp_plotFit(MRSCont, kk, which, stagFlag, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       which    = String for the spectrum to fit (optional)
%                   OPTIONS:    'off' (default)
%                               'diff1'
%                               'diff2'
%                               'sum'
%                               'ref'
%                               'w'
%       stagFlag  = flag to decide whether basis functions should be plotted
%                   vertically staggered or simply over one another
%                   (optional)
%                    OPTIONS:   1 = staggered (default)
%                               0 = not staggered
%       xlab      = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab      = label for the y-axis (optional.  Default = '');
%       figTitle  = label for the title of the plot (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-10-02)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-10-02: First version of the code.

% Check that OspreyFit has been run before
if ~MRSCont.flags.didFit
    error('Trying to plot fitted data, but fit has not been performed yet. Run OspreyFit first.')
end


%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
if nargin<7
    figTitle = sprintf('Fit plot') ;
    if nargin<6
        ylab='';
        if nargin<5
            xlab='Frequency (ppm)';
            if nargin<4
                stagFlag = 1;
                if nargin < 3
                    which = 'off';
                    if nargin < 2
                        kk = 1;
                        if nargin<1
                            error('ERROR: no input Osprey container specified.  Aborting!!');
                        end
                    end
                end
            end
        end
    end
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract processed spectra and fit parameters
if strcmp(which, 'off')
    dataToPlot  = MRSCont.processed.A{kk};
else
    dataToPlot  = MRSCont.processed.(which){kk};
end
if strcmp(which, 'ref') || strcmp(which, 'w')
    fitRangePPM = MRSCont.opts.fit.rangeWater;
    basisSet    = MRSCont.fit.resBasisSet.water;
else
    fitRangePPM = MRSCont.opts.fit.range;
    basisSet    = MRSCont.fit.resBasisSet.(which);
end
specToPlot  = op_freqrange(dataToPlot, fitRangePPM(1), fitRangePPM(end));
fitParams   = MRSCont.fit.results.(which).fitParams{kk};


%%% 3. PREPARE LINES TO DISPLAY %%%
% Extract raw and processed spectra in the plot range
% Concatenate parameters together into one large vector.
switch MRSCont.opts.fit.method
    case 'Osprey'
        ampl        = fitParams.ampl;
        zeroPhase   = fitParams.ph0;
        firstPhase  = fitParams.ph1;
        gaussLB     = fitParams.gaussLB;
        lorentzLB   = fitParams.lorentzLB;
        freqShift   = fitParams.freqShift;
        if strcmp(which, 'ref') || strcmp(which, 'w')
            % If water, extract and apply nonlinear parameters
            x       = [zeroPhase; firstPhase; gaussLB; lorentzLB; freqShift];
            [appliedBasisSet] = fit_applynlinwaterOsprey(basisSet, x, fitRangePPM);
            fitted      = (appliedBasisSet.specs*ampl);
            fitted      = fitted .* MRSCont.fit.scale{kk};
        else
            % If metabolites, extract and apply nonlinear parameters
            beta_j      = fitParams.beta_j;
            spl_pos     = fitParams.spl_pos;
            x          = [zeroPhase; firstPhase; gaussLB; lorentzLB; freqShift; beta_j];
            [appliedBasisSet, B] = fit_applynlinOsprey(basisSet, x, spl_pos, fitRangePPM);
            % Convert to complex baseline and scaled basis functions
            temp_B      = reshape(B, [length(B)/2, 2]);
            complex_B   = (temp_B(:,1) + 1i*temp_B(:,2));
            fitted      = (appliedBasisSet.specs*ampl + complex_B);
            % Scale
            complex_B   = complex_B .* MRSCont.fit.scale{kk};
            fitted      = fitted .* MRSCont.fit.scale{kk};
        end
end


%%% 4. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
out = figure;
% Prepare a couple of useful variables
nBasisFct = basisSet.nMets;
if isfield(basisSet, 'nMM')
    nBasisFct = nBasisFct + basisSet.nMM;
end


%%% 5. PLOT DATA, FIT, RESIDUAL, BASELINE %%%
% positive stagger to offset data, fit, residual, and baseline from the
% individual metabolite contributions
stagData = 0.2*(max(abs(min(real(dataToPlot.specs))), abs(max(real(dataToPlot.specs)))));
% Add the data and plot
hold on;
plot(specToPlot.ppm, real(specToPlot.specs) + stagData, 'k');
plot(specToPlot.ppm, real(fitted) + stagData, 'r', 'LineWidth', 1.5);
plot(specToPlot.ppm, real(specToPlot.specs) - real(fitted) + 1.2*max(real(specToPlot.specs)) + stagData, 'k', 'LineWidth', 1);
if ~(strcmp(which, 'ref') || strcmp(which, 'w'))
    plot(specToPlot.ppm, real(complex_B) + stagData, 'k', 'LineWidth', 1, 'Color', 'b');
end


%%% 6. PLOT BASIS FUNCTIONS %%%
% Plot separate metabolite basis functions only if not water
if ~(strcmp(which, 'ref') || strcmp(which, 'w'))
    if stagFlag
        % Staggered plots will be in all black and separated by the mean of the
        % maximum across all spectra
        stag = max(abs(mean(max(real(appliedBasisSet.specs)))), abs(mean(min(real(appliedBasisSet.specs))))) * MRSCont.fit.scale{kk};
        % Loop over all basis functions
        
        for rr = 1:nBasisFct
            plot(appliedBasisSet.ppm, ampl(rr).*appliedBasisSet.specs(:,rr).*MRSCont.fit.scale{kk} - rr*stag, 'k');
            % Instead of a MATLAB legend, annotate each line separately with the
            % name of the metabolite
            text(fitRangePPM(1), - rr*stag, appliedBasisSet.name{rr}, 'FontSize', 14);
        end
        
        % Preliminary formatting; might need some more stability here, or
        % differentiation based on sequence type
        set(gca, 'YLim', [-nBasisFct*stag-0.05*abs(min(real(dataToPlot.specs)))  Inf]);
        hold off;
        
    else
        % If not staggered, plots will simply be made with distinguishable
        % colors
        colours = distinguishable_colors(nBasisFct);
        
        % Loop over all basis functions
        for rr = 1:nBasisFct
            plot(appliedBasisSet.ppm, ampl(rr).*appliedBasisSet.specs(:,rr).*MRSCont.fit.scale{kk}, 'Color', colours(rr,:), 'LineWidth', 1);
        end
        legend([newline, appliedBasisSet.name], 'Orientation', 'horizontal');
        hold off
        
    end
    
else
    
    % If water is being shown, show a simple legend
    legend('Data', 'Fit', 'Residual');
    
end


%%% 7. DESIGN FINETUNING %%%
% Adapt common style for all axes
set(gca, 'XDir', 'reverse', 'XLim', [fitRangePPM(1), fitRangePPM(end)]);
set(gca, 'LineWidth', 1, 'TickDir', 'out');
set(gca, 'FontSize', 16);
% If no y caption, remove y axis
if isempty(ylab)
    set(gca, 'YColor', 'w');
else
    set(gca, 'YColor', 'k');
end
% Black axes, white background
set(gca, 'XColor', 'k');
set(gca, 'Color', 'w');
set(gcf, 'Color', 'w');
box off;
title(figTitle);
xlabel(xlab, 'FontSize', 16);
ylabel(ylab, 'FontSize', 16);


%%% 8. ADD OSPREY LOGO %%%
[I, map] = imread('osprey.gif','gif');
axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
imshow(I, map);
axis off;


end

   