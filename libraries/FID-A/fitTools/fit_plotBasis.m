% fit_plotBasis.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% out = fit_plotBasis(basisSet, dim, stagFlag, ppmmin, ppmmax, xlab, ylab, figTitle)
% 
% DESCRIPTION:
% Creates a figure showing all a basis functions in FID-A basis set struct
% (generated with fit_makeBasis).
% 
% INPUTS:
% basisSet  = Simulated basis set in FID-A structure format.
% dim       = Dimension of the basis set
% stagFlag  = flag to decide whether basis functions should be plotted
%               vertically staggered or simply over one another (optional.
%               Default: 1 = staggered; 0 = not staggered)
% ppmmin    = lower limit of ppm scale to plot (optional.  Default = 0.2 ppm).
% ppmmax    = upper limit of ppm scale to plot (optional.  Default = 5.2 ppm).
% xlab      = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
% ylab      = label for the y-axis (optional.  Default = '');
% figTitle  = label for the title of the plot (optional.  Default = '');


function out = fit_plotBasis(basisSet, dim, stagFlag, ppmmin, ppmmax, xlab, ylab, figTitle)

% Parse input arguments
if nargin<8
    figTitle = sprintf('Basis set for %s at TE = %i ms', basisSet.seq{1}, basisSet.te) ;
    if nargin<7
        ylab='';
        if nargin<6
            xlab='Frequency (ppm)';
            if nargin<5
                ppmmax=5.2;
                if nargin<4
                    ppmmin=0.2;
                    if nargin<3
                        stagFlag = 1;
                        if nargin < 2
                            if length(size(basisSet.fids)) > 2
                                disp('More than one basis set dimension, but no number to display specified. Displaying first.');
                            end
                            dim = 1;
                            if nargin<1
                                error('ERROR: no input basis set specified.  Aborting!!');
                            end
                        end
                    end
                end
            end
        end
    end
end

% Generate a new figure and keep the handle memorized
out = figure;
% Prepare a couple of useful variables
nBasisFct = basisSet.nMets;
if isfield(basisSet, 'nMM')
    nBasisFct = nBasisFct + basisSet.nMM;
end

% Plot the basis functions
if stagFlag
    % Staggered plots will be in all black and separated by the mean of the
    % maximum across all spectra
    stag = mean(max(real(basisSet.specs(:,:,dim))));
    
    % Loop over all basis functions
    hold on
    for kk = 1:nBasisFct
        plot(basisSet.ppm, squeeze(basisSet.specs(:,kk,dim)) - kk*stag, 'k');
        % Instead of a MATLAB legend, annotate each line separately with the
        % name of the metabolite
        text(ppmmin, - kk*stag, basisSet.name{kk}, 'FontSize', 14);
    end
    hold off
    
    
    
else
    % If not staggered, plots will simply be made with distinguishable
    % colors
    colours = distinguishable_colors(nBasisFct);
    
    % Loop over all basis functions
    hold on
    for kk = 1:nBasisFct
        plot(basisSet.ppm, squeeze(basisSet.specs(:,kk,dim)), 'Color', colours(kk,:));
    end
    legend(basisSet.name);
    hold off
end

% Common style for all outputs
set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax]);
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
set(Fig1Ax1Line1, 'LineWidth', 1);


end

   