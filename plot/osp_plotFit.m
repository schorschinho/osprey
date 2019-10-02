% fit_plotBasis.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% out = fit_plotBasis(basisSet, stagFlag, ppmmin, ppmmax, xlab, ylab, figTitle)
% 
% DESCRIPTION:
% Creates a figure showing all a basis functions in FID-A basis set struct
% (generated with fit_makeBasis).
% 
% INPUTS:
% MRSCont  = Osprey data container.
% kk       = Index for the kk-th dataset
% which    = String for the spectrum to fit
%            OPTIONS:   'off' (default)
%                       'diff1'
%                       'diff2'
%                       'sum'
% stagFlag  = flag to decide whether basis functions should be plotted
%               vertically staggered or simply over one another (optional.
%               Default: 1 = staggered; 0 = not staggered)
% xlab      = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
% ylab      = label for the y-axis (optional.  Default = '');
% figTitle  = label for the title of the plot (optional.  Default = '');


function out = osp_plotFit(MRSCont, kk, which, stagFlag, xlab, ylab, figTitle)

%% out = osp_plotFit(MRSCont, kk, which, stagFlag, xlab, ylab, figTitle)
%   Creates a figure showing data stored in an Osprey data container, as
%   well as the fit to it, the baseline, the residual, and contributions
%   from the individual metabolites.
%
%   USAGE:
%       out = fit_plotBasis(MRSCont, kk, which, stagFlag, xlab, ylab, figTitle)
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

% Parse input arguments and fall back to defaults if not provided
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

% Extract parameters
if strcmp(which, 'off')
    dataToPlot  = MRSCont.processed.A{kk};
else
    dataToPlot  = MRSCont.processed.(which){kk};
end
fitRangePPM = MRSCont.opts.fit.range;
specToPlot  = op_freqrange(dataToPlot, fitRangePPM(1), fitRangePPM(end));
basisSet    = MRSCont.fit.resBasisSet.(which);
ampl        = MRSCont.fit.results.(which).ampl{kk};
zeroPhase   = MRSCont.fit.results.(which).ph0{kk};
firstPhase  = MRSCont.fit.results.(which).ph1{kk};
gaussLB     = MRSCont.fit.results.(which).gaussLB{kk};
lorentzLB   = MRSCont.fit.results.(which).lorentzLB{kk};
freqShift   = MRSCont.fit.results.(which).freqShift{kk};
beta_j      = MRSCont.fit.results.(which).beta_j{kk};
spl_pos     = MRSCont.fit.results.(which).spl_pos{kk};

% Concatenate parameters together into one large vector.
x          = [zeroPhase; firstPhase; gaussLB; lorentzLB; freqShift; beta_j];
% Apply nonlinear parameters
[basisSet, B] = fit_applyMetabModel(basisSet, x, spl_pos, fitRangePPM);

% Convert to complex baseline and scaled basis functions
temp_B      = reshape(B, [length(B)/2, 2]);
complex_B   = (temp_B(:,1) + 1i*temp_B(:,2));
fitted      = (basisSet.specs*ampl + complex_B);
% Scale
complex_B   = complex_B .* MRSCont.fit.scale{kk};
fitted      = fitted .* MRSCont.fit.scale{kk};

% Generate a new figure and keep the handle memorized
out = figure;
% Prepare a couple of useful variables
nBasisFct = basisSet.nMets;
if isfield(basisSet, 'nMM')
    nBasisFct = nBasisFct + basisSet.nMM;
end

% Plot the data, fit, residual, and baseline
% positive stagger
stagData = 0.2*(max(abs(min(real(dataToPlot.specs))), abs(max(real(dataToPlot.specs)))));
% Add the data and plot
hold on;
plot(specToPlot.ppm, real(specToPlot.specs) + stagData, 'k');
plot(specToPlot.ppm, real(fitted) + stagData, 'r', 'LineWidth', 1.5);
plot(specToPlot.ppm, real(complex_B) + stagData, 'k', 'LineWidth', 1);
plot(specToPlot.ppm, real(specToPlot.specs) - real(fitted) + 1.2*max(real(specToPlot.specs)) + stagData, 'k', 'LineWidth', 1);
    
% Plot the basis functions
if stagFlag
    % Staggered plots will be in all black and separated by the mean of the
    % maximum across all spectra
    stag = max(abs(mean(max(real(basisSet.specs)))), abs(mean(min(real(basisSet.specs))))) * MRSCont.fit.scale{kk};
    % Loop over all basis functions

    for rr = 1:nBasisFct
        plot(basisSet.ppm, ampl(rr).*basisSet.specs(:,rr).*MRSCont.fit.scale{kk} - rr*stag, 'k');
        % Instead of a MATLAB legend, annotate each line separately with the
        % name of the metabolite
        text(fitRangePPM(1), - rr*stag, basisSet.name{rr}, 'FontSize', 14);
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
        plot(basisSet.ppm, ampl(rr).*basisSet.specs(:,rr).*MRSCont.fit.scale{kk}, 'Color', colours(rr,:), 'LineWidth', 1);
    end
    legend(['Data', 'Fit', 'Baseline', 'Residual', newline, basisSet.name], 'Orientation', 'horizontal' );
    hold off
    
end

% Common style for all outputs
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
xlabel(xlab, 'FontSize', 20);
ylabel(ylab, 'FontSize', 20);
% Set linewidth coherently
Fig1Ax1 = get(out, 'Children');
Fig1Ax1Line1 = get(Fig1Ax1, 'Children');
if iscell(Fig1Ax1Line1)
    Fig1Ax1Line1 = Fig1Ax1Line1(~cellfun('isempty', Fig1Ax1Line1));
    Fig1Ax1Line1 = Fig1Ax1Line1{1};
end
%set(Fig1Ax1Line1, 'LineWidth', 1);


end

   