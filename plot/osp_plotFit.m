function out = osp_plotFit(MRSCont, kk, which, GUI, stagFlag, xlab, ylab, figTitle)
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
%       GUI       = flag to decide whether plot is used in GUI
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
if nargin<8    
    [~,filen,ext] = fileparts(MRSCont.files{kk});
    figTitle = sprintf([MRSCont.opts.fit.method ' ' MRSCont.opts.fit.style ' ' which ' fit plot:\n' filen ext]);
    if nargin<7
        ylab='';
        if nargin<6
            xlab='Frequency (ppm)';
            if nargin<5
                stagFlag = 1;
                if nargin<4
                    GUI = 0;    
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
    basisSet    = MRSCont.fit.resBasisSet.(which){kk};
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
stagData = 0.1*(max(abs(min(real(dataToPlot.specs))), abs(max(real(dataToPlot.specs)))));
maxPlot = max(real(specToPlot.specs) +  abs(min(real(specToPlot.specs) - real(fitted)))) + abs(max(real(specToPlot.specs) - real(fitted))) + stagData;
% Add the data and plot
hold on;
if ~GUI
    plot(specToPlot.ppm, (zeros(1,length(specToPlot.ppm)) + stagData)/maxPlot, 'k'); % Zeroline
    plot(specToPlot.ppm, (real(specToPlot.specs) + stagData)/maxPlot, 'k'); % Data
    plot(specToPlot.ppm, (real(fitted) + stagData)/maxPlot, 'r', 'LineWidth', 1.5); % Fit
    plot(specToPlot.ppm, (zeros(1,length(specToPlot.ppm)) + max(real(specToPlot.specs)) + stagData)/maxPlot, 'k', 'LineWidth', 1); % Maximum Data

    plot(specToPlot.ppm, (real(specToPlot.specs) - real(fitted) + max(real(specToPlot.specs) +  abs(min(real(specToPlot.specs) - real(fitted)))) + stagData)/maxPlot, 'k', 'LineWidth', 1); % Residue
    plot(specToPlot.ppm, (zeros(1,length(specToPlot.ppm)) + max(real(specToPlot.specs) +  abs(min(real(specToPlot.specs) - real(fitted)))) + stagData)/maxPlot, '--k', 'LineWidth', 0.5); % Zeroline Residue
    plot(specToPlot.ppm, (zeros(1,length(specToPlot.ppm)) + max(real(specToPlot.specs) +  abs(min(real(specToPlot.specs) - real(fitted)))) + abs(max(real(specToPlot.specs) - real(fitted))) + stagData)/maxPlot, 'k', 'LineWidth', 1); % Max Residue
else
    plot(specToPlot.ppm, (zeros(1,length(specToPlot.ppm)) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground); % Zeroline
    plot(specToPlot.ppm, (real(specToPlot.specs) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground); % Data
    plot(specToPlot.ppm, (real(fitted) + stagData)/maxPlot, 'Color', MRSCont.colormap.Accent, 'LineWidth', 1.6); % Fit
    plot(specToPlot.ppm, (zeros(1,length(specToPlot.ppm)) + max(real(specToPlot.specs)) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Maximum Data

    plot(specToPlot.ppm, (real(specToPlot.specs) - real(fitted) + max(real(specToPlot.specs) +  abs(min(real(specToPlot.specs) - real(fitted)))) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Residue
    plot(specToPlot.ppm, (zeros(1,length(specToPlot.ppm)) + max(real(specToPlot.specs) +  abs(min(real(specToPlot.specs) - real(fitted)))) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineStyle','--', 'LineWidth', 0.5); % Zeroline Residue
    plot(specToPlot.ppm, (zeros(1,length(specToPlot.ppm)) + max(real(specToPlot.specs) +  abs(min(real(specToPlot.specs) - real(fitted)))) + abs(max(real(specToPlot.specs) - real(fitted))) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Max Residue 
end

if ~GUI
    text(fitRangePPM(1), (0 + stagData)/maxPlot, '0', 'FontSize', 14); %Zeroline text
    text(fitRangePPM(1), (0 + max(real(specToPlot.specs)) + stagData)/maxPlot -0.05, num2str(max(real(specToPlot.specs)),'%1.2e'), 'FontSize', 14); % Maximum Data Text
    text(fitRangePPM(1), (0 + max(real(specToPlot.specs) +  abs(min(real(specToPlot.specs) - real(fitted)))) + stagData)/maxPlot, '0', 'FontSize', 14); %Zeroline Residue text
    text(fitRangePPM(1), (0 + max(real(specToPlot.specs) +  abs(min(real(specToPlot.specs) - real(fitted)))) + abs(max(real(specToPlot.specs) - real(fitted))) + stagData)/maxPlot+0.05, num2str(abs(max(real(specToPlot.specs) - real(fitted))),'%1.2e'), 'FontSize', 14); %Max Residue text
else
    text(fitRangePPM(1), (0 + stagData)/maxPlot, '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline text
    text(fitRangePPM(1), (0 + max(real(specToPlot.specs)) + stagData)/maxPlot-0.05, num2str(max(real(specToPlot.specs)),'%1.2e'), 'FontSize', 10,'Color',MRSCont.colormap.Foreground); % Maximum Data Text
    text(fitRangePPM(1), (0 + max(real(specToPlot.specs) +  abs(min(real(specToPlot.specs) - real(fitted)))) + stagData)/maxPlot, '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline Residue text
    text(fitRangePPM(1), (0 + max(real(specToPlot.specs) +  abs(min(real(specToPlot.specs) - real(fitted)))) + abs(max(real(specToPlot.specs) - real(fitted))) + stagData)/maxPlot +0.05, num2str(abs(max(real(specToPlot.specs) - real(fitted))),'%1.2e'), 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Max Residue text
end

if ~(strcmp(which, 'ref') || strcmp(which, 'w'))
    if ~GUI
        plot(specToPlot.ppm, (real(complex_B) + stagData)/maxPlot, 'k', 'LineWidth', 1, 'Color', 'b');
    else
        plot(specToPlot.ppm, (real(complex_B) + stagData)/maxPlot, 'k', 'LineWidth', 1, 'Color', MRSCont.colormap.LightAccent);
    end
end


%%% 6. PLOT BASIS FUNCTIONS %%%
% Plot separate metabolite basis functions only if not water
if ~(strcmp(which, 'ref') || strcmp(which, 'w'))
    if stagFlag
        % Staggered plots will be in all black and separated by the mean of the
        % maximum across all spectra
%         stag = max(abs(mean(max(real(appliedBasisSet.specs)))), abs(mean(min(real(appliedBasisSet.specs))))) * MRSCont.fit.scale{kk};
        stag = maxPlot *  2.5 / nBasisFct;
        % Loop over all basis functions
        
        for rr = 1:nBasisFct
            % Instead of a MATLAB legend, annotate each line separately with the
            % name of the metabolite
            if ~GUI
                plot(appliedBasisSet.ppm, (ampl(rr).*appliedBasisSet.specs(:,rr).*MRSCont.fit.scale{kk} - rr*stag)/maxPlot, 'k');
                text(fitRangePPM(1), (- rr*stag)/maxPlot, appliedBasisSet.name{rr}, 'FontSize', 14);
            else
                plot(appliedBasisSet.ppm, (ampl(rr).*appliedBasisSet.specs(:,rr).*MRSCont.fit.scale{kk} - rr*stag)/maxPlot, 'Color',MRSCont.colormap.Foreground);
                text(fitRangePPM(1), (- rr*stag)/maxPlot, appliedBasisSet.name{rr}, 'FontSize', 10,'Color',MRSCont.colormap.Foreground);
            end
        end
        
        % Preliminary formatting; might need some more stability here, or
        % differentiation based on sequence type
        set(gca, 'YLim', [(-nBasisFct-1)*stag/maxPlot  1]);
        hold off;
        
    else
        % If not staggered, plots will simply be made with distinguishable
        % colors
        colours = distinguishable_colors(nBasisFct);
        
        % Loop over all basis functions
        for rr = 1:nBasisFct
            plot(appliedBasisSet.ppm, (ampl(rr).*appliedBasisSet.specs(:,rr).*MRSCont.fit.scale{kk})/maxPlot, 'Color', colours(rr,:), 'LineWidth', 1);
        end
        legend([newline, appliedBasisSet.name], 'Orientation', 'horizontal');
        hold off
        
    end
    
else
    % Preliminary formatting; might need some more stability here, or
    % differentiation based on sequence type
    set(gca, 'YLim', [0  1]);
    hold off;
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
    if ~GUI
        set(gca, 'YColor', 'w');
        % Black axes, white background
        set(gca, 'XColor', 'k');
        set(gca, 'Color', 'w');
        set(gcf, 'Color', 'w');
        title(figTitle, 'Interpreter', 'none');
    else
        set(gca, 'YColor', MRSCont.colormap.Background);
        set(gca,'YTickLabel',{});
        set(gca,'YTick',{});
        % Dirtywhite axes, light gray background
        set(gca, 'XColor', MRSCont.colormap.Foreground);
        set(gca, 'Color', MRSCont.colormap.Background);
        set(gcf, 'Color', MRSCont.colormap.Background);
        title(figTitle, 'Interpreter', 'none', 'Color', MRSCont.colormap.Foreground);
    end
else
    set(gca, 'YColor', 'k');
end

box off;
xlabel(xlab, 'FontSize', 16);
ylabel(ylab, 'FontSize', 16);


%%% 8. ADD OSPREY LOGO %%%
if ~GUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end

end

   