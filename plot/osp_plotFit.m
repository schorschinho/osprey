function out = osp_plotFit(MRSCont, kk, which, GUI, conc, stagFlag, xlab, ylab, figTitle)
%% out = osp_plotFit(MRSCont, kk, which, stagFlag, xlab, ylab, figTitle)
%   Creates a figure showing data stored in an Osprey data container, as
%   well as the fit to it, the baseline, the residual, and contributions
%   from the individual metabolites.
%
%   USAGE:
%       out = osp_plotFit(MRSCont, kk, which, GUI, conc, stagFlag, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       which    = String for the spectrum to plot (optional)
%                   OPTIONS:    'off' (default)
%                               'diff1'
%                               'diff2'
%                               'sum'
%                               'ref'
%                               'w'
%       GUI       = flag to decide whether plot is used in GUI
%       conc      = 
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
% Get the fit method and style
fitMethod   = MRSCont.opts.fit.method;
fitStyle    = MRSCont.opts.fit.style;
% Fall back to defaults if not provided
if nargin<9    
    [~,filen,ext] = fileparts(MRSCont.files{kk});
    if strcmp(which, 'conc')
        figTitle = sprintf([fitMethod ' ' fitStyle ' ' conc ' fit plot:\n' filen ext]);
    else
        figTitle = sprintf([fitMethod ' ' fitStyle ' ' which ' fit plot:\n' filen ext]);
    end
    if nargin<8
        ylab='';
        if nargin<7
            xlab='Frequency (ppm)';
            if nargin<6
                stagFlag = 1;
                if nargin<5
                    conc = 'diff1';
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
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract processed spectra and fit parameters
if  strcmp(which, 'conc')
    dataToPlot  = MRSCont.processed.(conc){kk};
else
    if strcmp(which, 'off')
        dataToPlot  = MRSCont.processed.A{kk};
    else
        dataToPlot  = MRSCont.processed.(which){kk};
    end
end
if strcmp(which, 'ref') || strcmp(which, 'w')
    fitRangePPM = MRSCont.opts.fit.rangeWater;
    basisSet    = MRSCont.fit.resBasisSet.(which).water{kk};
else if strcmp(which, 'conc')
        fitRangePPM = MRSCont.opts.fit.range;
        basisSet    = MRSCont.fit.resBasisSet.(which){kk};
    else
        fitRangePPM = MRSCont.opts.fit.range;
        basisSet    = MRSCont.fit.resBasisSet.(which){kk};
    end
end

% Get the fit parameters
fitParams   = MRSCont.fit.results.(which).fitParams{kk};
% Pack up into structs to feed into the reconstruction functions
inputData.dataToFit                 = dataToPlot;
inputData.basisSet                  = basisSet;
inputSettings.scale                 = MRSCont.fit.scale{kk};
inputSettings.fitRangePPM           = fitRangePPM;
inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;


%%% 3. PREPARE LINES TO DISPLAY %%%
% Extract data, ppm axes, fit, residual, baseline, and individual
% metabolite contributions.
% Obviously, depending on the fit method that has been used, we need to
% reconstruct the plots with fit method functions
switch fitMethod
    % Depending on whether metabolite or water data are to be
    % displayed, create the plots via different models
    case 'Osprey'
        if strcmp(which, 'ref') || strcmp(which, 'w')
            % if water, use the water model
            [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
        else
            % if metabolites, use the metabolite model
            [ModelOutput] = fit_OspreyParamsToModel(inputData, inputSettings, fitParams);
        end
    case 'LCModel'
        if strcmp(which, 'ref') || strcmp(which, 'w')
            % if water, use the water model
            [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
        else
            % if metabolites, use the metabolite model
            [ModelOutput] = fit_LCModelParamsToModel(inputData, inputSettings, fitParams);
        end
        % leave the concatenated case out for now?
        
%         else if strcmp(which, 'conc')
%                 % If concatened metabolites, extract and apply nonlinear
%                 % parameters and shifts
%                 beta_j      = fitParams.beta_j;
%                 spl_pos     = fitParams.spl_pos;
%                 addFreqShift = fitParams.addFreqShift;
%                 x          = [zeroPhase; firstPhase; gaussLB; lorentzLB; freqShift; beta_j; addFreqShift];
%                 [appliedBasisSet, B] = fit_applynlinOspreyMultiplex(basisSet, x, spl_pos, fitRangePPM);
%                 % Convert to complex baseline and scaled basis functions
%                 switch conc
%                     case 'diff1'
%                         B_conc = B(:,1);
%                         appliedBasisSet.specs = appliedBasisSet.specs(:,:,1);
%                     case 'sum'
%                         B_conc = B(:,2);
%                         appliedBasisSet.specs = appliedBasisSet.specs(:,:,2);
%                     end
%                 temp_B      = reshape(B_conc, [length(B_conc)/2, 2]);
%                 baseline   = (temp_B(:,1) + 1i*temp_B(:,2));
%                 fitted      = (appliedBasisSet.specs*ampl + baseline);
%                 % Scale
%                 baseline   = baseline .* MRSCont.fit.scale{kk};
%                 fitted      = fitted .* MRSCont.fit.scale{kk};
%             end
%         end
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
% Unpack the ModelOutput struct
ppm         = ModelOutput.ppm;
dataToPlot  = ModelOutput.data;
fit         = ModelOutput.completeFit;
residual    = ModelOutput.residual;
% If water, don't get baseline and individual fits
if ~(strcmp(which, 'ref') || strcmp(which, 'w'))
    baseline    = ModelOutput.baseline;
    indivPlots  = ModelOutput.indivMets;
end

% Determine a positive stagger to offset data, fit, residual, and 
% baseline from the individual metabolite contributions
stagData = 0.1*(max(abs(min(dataToPlot)), abs(max(dataToPlot))));
maxPlot = max(dataToPlot + abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData;
% Add the data and plot
hold on;
if ~GUI
    plot(ppm, (zeros(1,length(ppm)) + stagData)/maxPlot, 'k'); % Zeroline
    plot(ppm, (dataToPlot + stagData)/maxPlot, 'k'); % Data
    plot(ppm, (fit + stagData)/maxPlot, 'r', 'LineWidth', 1.5); % Fit
    plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot) + stagData)/maxPlot, 'k', 'LineWidth', 1); % Maximum Data

    plot(ppm, (residual + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, 'k', 'LineWidth', 1); % Residual
    plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, '--k', 'LineWidth', 0.5); % Zeroline Residue
    plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot +  abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData)/maxPlot, 'k', 'LineWidth', 1); % Max Residue
else
    plot(ppm, (zeros(1,length(ppm)) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground); % Zeroline
    plot(ppm, (dataToPlot + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground); % Data
    plot(ppm, (fit + stagData)/maxPlot, 'Color', MRSCont.colormap.Accent, 'LineWidth', 1.6); % Fit
    plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Maximum Data

    plot(ppm, (residual + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Residual
    plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineStyle','--', 'LineWidth', 0.5); % Zeroline Residue
    plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot +  abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Max Residue 
end

if ~GUI
    text(fitRangePPM(1), (0 + stagData)/maxPlot, '0', 'FontSize', 14); %Zeroline text
    text(fitRangePPM(1), (0 + max(dataToPlot) + stagData)/maxPlot -0.05, num2str(max(dataToPlot),'%1.2e'), 'FontSize', 14); % Maximum Data Text
    text(fitRangePPM(1), (0 + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, '0', 'FontSize', 14); %Zeroline Residua; text
    text(fitRangePPM(1), (0 + max(dataToPlot +  abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData)/maxPlot+0.05, num2str(abs(max(dataToPlot - fit)),'%1.2e'), 'FontSize', 14); %Max Residue text
else
    text(fitRangePPM(1), (0 + stagData)/maxPlot, '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline text
    text(fitRangePPM(1), (0 + max(dataToPlot) + stagData)/maxPlot-0.05, num2str(max(dataToPlot),'%1.2e'), 'FontSize', 10,'Color',MRSCont.colormap.Foreground); % Maximum Data Text
    text(fitRangePPM(1), (0 + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline Residual text
    text(fitRangePPM(1), (0 + max(dataToPlot +  abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData)/maxPlot +0.05, num2str(abs(max(dataToPlot - fit)),'%1.2e'), 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Max Residue text
end

if ~(strcmp(which, 'ref') || strcmp(which, 'w'))
    if ~GUI
        plot(ppm, (real(baseline) + stagData)/maxPlot, 'k', 'LineWidth', 1, 'Color', 'b');
    else
        plot(ppm, (real(baseline) + stagData)/maxPlot, 'k', 'LineWidth', 1, 'Color', MRSCont.colormap.LightAccent);
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
                plot(ppm, (indivPlots(:,rr) - rr*stag)/maxPlot, 'k');
                text(fitRangePPM(1), (- rr*stag)/maxPlot, basisSet.name{rr}, 'FontSize', 14);
            else
                plot(ppm, (indivPlots(:,rr) - rr*stag)/maxPlot, 'Color',MRSCont.colormap.Foreground);
                text(fitRangePPM(1), (- rr*stag)/maxPlot, basisSet.name{rr}, 'FontSize', 10,'Color',MRSCont.colormap.Foreground);
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
            plot(ppm, indivPlots(:,rr)/maxPlot, 'Color', colours(rr,:), 'LineWidth', 1);
        end
        legend([newline, basisSet.name], 'Orientation', 'horizontal');
        hold off
        
    end
    
else
    % Preliminary formatting; might need some more stability here, or
    % differentiation based on sequence type
    set(gca, 'YLim', [0  1]);
    hold off;
    % If water is being shown, show a simple legend
    %legend('Data', 'Fit', 'Residual');
    
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

   