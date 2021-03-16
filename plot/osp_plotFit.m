function out = osp_plotFit(MRSCont, kk, which_spec,VoxelIndex, conc, stagFlag, xlab, ylab, figTitle)
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
%                                 'mm' re_mm
%       conc      = flag if concatenate fitting was used
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
    if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
        if strcmp(which_spec, 'conc')
            figTitle = sprintf([fitMethod ' ' fitStyle ' ' conc ' fit plot:\n' filen ext]);
        else
            figTitle = sprintf([fitMethod ' ' fitStyle ' ' which_spec ' fit plot:\n' filen ext]);
        end
    elseif  (MRSCont.flags.isPRIAM && isfield(MRSCont.flags,'isPRIAM')) 
        if nargin<4
            VoxelIndex = 1; 
         end
        if strcmp(which_spec, 'conc')
            figTitle = sprintf([fitMethod ' ' fitStyle ' ' conc ' fit plot:\n' filen ext '\n Voxel ' num2str(VoxelIndex)]);
        else
            figTitle = sprintf([fitMethod ' ' fitStyle ' ' which_spec ' fit plot:\n' filen ext  '\n Voxel ' num2str(VoxelIndex)]);
        end 
    else
            if nargin<4
                VoxelIndex = [1 1]; 
             end
            if strcmp(which_spec, 'conc')
                figTitle = sprintf([fitMethod ' ' fitStyle ' ' conc ' fit plot:\n' filen ext '\n Voxel ' num2str(VoxelIndex(1)) ' ' num2str(VoxelIndex(2))]);
            else
                figTitle = sprintf([fitMethod ' ' fitStyle ' ' which_spec ' fit plot:\n' filen ext  '\n Voxel ' num2str(VoxelIndex(1)) ' ' num2str(VoxelIndex(2))]);
            end 
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
                        VoxelIndex = 1; 
                        if nargin < 3
                            which_spec = 'off';
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
if (MRSCont.flags.isPRIAM == 1) || (MRSCont.flags.isMRSI == 1)
    if ~exist('VoxelIndex') && (MRSCont.flags.isPRIAM == 1)
            VoxelIndex = 1;
        elseif ~exist('VoxelIndex') && (MRSCont.flags.isMRSI == 1)
            VoxelIndex = [1 1];  
        end
    if  strcmp(which_spec, 'conc')
        dataToPlot=op_takeVoxel(MRSCont.processed.(conc){kk},VoxelIndex);
    else
        if strcmp(which_spec, 'off')
            dataToPlot=op_takeVoxel(MRSCont.processed.A{kk},VoxelIndex);
        else
            dataToPlot=op_takeVoxel(MRSCont.processed.(which_spec){kk},VoxelIndex);
        end
    end

    if (MRSCont.flags.isPRIAM == 1)
        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            fitRangePPM = MRSCont.opts.fit.rangeWater;
            basisSet    = MRSCont.fit.resBasisSet{VoxelIndex}.(which_spec).water{MRSCont.info.(which_spec).unique_ndatapoint_indsort(kk)};
        else if strcmp(which_spec, 'conc')
                fitRangePPM = MRSCont.opts.fit.range;
                basisSet    = MRSCont.fit.resBasisSet{VoxelIndex}.(which_spec){MRSCont.info.diff1.unique_ndatapoint_indsort(kk)};
            else if strcmp(which_spec, 'off')
                    fitRangePPM = MRSCont.opts.fit.range;
                    basisSet    = MRSCont.fit.resBasisSet{VoxelIndex}.(which_spec){kk};
                else
                    fitRangePPM = MRSCont.opts.fit.range;
                    basisSet    = MRSCont.fit.resBasisSet{VoxelIndex}.(which_spec){kk};
                end
            end
        end
    else
        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            fitRangePPM = MRSCont.opts.fit.rangeWater;
            basisSet    = MRSCont.fit.resBasisSet{VoxelIndex(1), VoxelIndex(2)}.(which_spec).water{MRSCont.info.(which_spec).unique_ndatapoint_indsort(kk)};
        else if strcmp(which_spec, 'conc')
                fitRangePPM = MRSCont.opts.fit.range;
                basisSet    = MRSCont.fit.resBasisSet{VoxelIndex(1), VoxelIndex(2)}.(which_spec){MRSCont.info.diff1.unique_ndatapoint_indsort(kk)};
            else if strcmp(which_spec, 'off')
                    fitRangePPM = MRSCont.opts.fit.range;
                    basisSet    = MRSCont.fit.resBasisSet{VoxelIndex(1), VoxelIndex(2)}.(which_spec){kk};
                else
                    fitRangePPM = MRSCont.opts.fit.range;
                    basisSet    = MRSCont.fit.resBasisSet{VoxelIndex(1), VoxelIndex(2)}.(which_spec){kk};
                end
            end
        end        
    end
else
    if  strcmp(which_spec, 'conc')
        dataToPlot  = MRSCont.processed.(conc){kk};
    else
        if strcmp(which_spec, 'off')
            dataToPlot  = MRSCont.processed.A{kk};
        else
            dataToPlot  = MRSCont.processed.(which_spec){kk};
        end
    end


    if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
        fitRangePPM = MRSCont.opts.fit.rangeWater;
        basisSet    = MRSCont.fit.resBasisSet.(which_spec).water{MRSCont.info.(which_spec).unique_ndatapoint_indsort(kk)};
    else if strcmp(which_spec, 'conc')
            fitRangePPM = MRSCont.opts.fit.range;
            basisSet    = MRSCont.fit.resBasisSet.(which_spec){MRSCont.info.diff1.unique_ndatapoint_indsort(kk)};
        else if strcmp(which_spec, 'off')
                fitRangePPM = MRSCont.opts.fit.range;
                basisSet    = MRSCont.fit.resBasisSet.(which_spec){kk};
            else
                fitRangePPM = MRSCont.opts.fit.range;
                basisSet    = MRSCont.fit.resBasisSet.(which_spec){kk};
            end
        end
    end
end



% Get the fit parameters

if (MRSCont.flags.isPRIAM == 1)
    fitParams   = MRSCont.fit.results{VoxelIndex}.(which_spec).fitParams{kk};
elseif (MRSCont.flags.isMRSI == 1)
    fitParams   = MRSCont.fit.results{VoxelIndex(1), VoxelIndex(2)}.(which_spec).fitParams{kk};
else
    fitParams   = MRSCont.fit.results.(which_spec).fitParams{kk};
end
% Pack up into structs to feed into the reconstruction functions
inputData.dataToFit                 = dataToPlot;
inputData.basisSet                  = basisSet;
if (length(fitParams.ampl) == 3) 
    inputData.basisSet_mm                  = MRSCont.fit.basisSet_mm;
end
if (MRSCont.flags.isPRIAM == 1)
    inputSettings.scale                 = MRSCont.fit.scale{kk};
else
    inputSettings.scale                 = MRSCont.fit.scale{kk};
end
inputSettings.fitRangePPM           = fitRangePPM;
inputSettings.minKnotSpacingPPM     = MRSCont.opts.fit.bLineKnotSpace;
inputSettings.fitStyle              = MRSCont.opts.fit.style;
inputSettings.flags.isMEGA          = MRSCont.flags.isMEGA;
inputSettings.flags.isHERMES        = MRSCont.flags.isHERMES;
inputSettings.flags.isHERCULES      = MRSCont.flags.isHERCULES;
inputSettings.flags.isPRIAM         = MRSCont.flags.isPRIAM;
inputSettings.concatenated.Subspec  = conc;

%%% 3. PREPARE LINES TO DISPLAY %%%
% Extract data, ppm axes, fit, residual, baseline, and individual
% metabolite contributions.
% Obviously, depending on the fit method that has been used, we need to
% reconstruct the plots with fit method functions
switch fitMethod
    % Depending on whether metabolite or water data are to be
    % displayed, create the plots via different models
    case 'Osprey'
        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            % if water, use the water model
            [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
        else
            % if metabolites, use the metabolite model
            if strcmp(inputSettings.fitStyle,'Concatenated')
                [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams);
            else
                [ModelOutput] = fit_OspreyParamsToModel(inputData, inputSettings, fitParams);
            end
        end
    case 'OspreyAsym'
        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            % if water, use the water model
            [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
        else
            % if metabolites, use the metabolite model
            if strcmp(inputSettings.fitStyle,'Concatenated')
                [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams);
            else
                [ModelOutput] = fit_OspreyAsymParamsToModel(inputData, inputSettings, fitParams);
            end
        end        
    case 'OspreyNoLS'
        if strcmp(which_spec, 'ref') || strcmp(which_spec, 'w')
            % if water, use the water model
            [ModelOutput] = fit_waterOspreyParamsToModel(inputData, inputSettings, fitParams);
        else
            % if metabolites, use the metabolite model
            if strcmp(inputSettings.fitStyle,'Concatenated')
                [ModelOutput] = fit_OspreyParamsToConcModel(inputData, inputSettings, fitParams);
            else
                [ModelOutput] = fit_OspreyNoLSParamsToModel(inputData, inputSettings, fitParams);
            end
        end        
end

%re_mm
%For MM, prepare a 'clean' MM spectrum that has metabolite signals pulled
%out.
if (strcmp(which_spec, 'mm'))
   Met_corr_spectrum  = sum(ModelOutput.indivMets(:,1:4),2);
end


%%% 4. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
canvasSize  = get(0,'defaultfigureposition');
if stagFlag && ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
    canvasSize(4) = canvasSize(4) * 1.8;
end
out = figure('Position', canvasSize);
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
if ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
    baseline    = ModelOutput.baseline;
    indivPlots  = ModelOutput.indivMets;
end

if isfield(MRSCont.plot,'fit') && MRSCont.plot.fit.match
    if strcmp(which_spec, 'conc')
        stagData = MRSCont.plot.fit.(conc).stagData(kk);
        maxPlot = MRSCont.plot.fit.(conc).maxPlot(kk);        
    else
        stagData = MRSCont.plot.fit.(which_spec).stagData(kk);
        maxPlot = MRSCont.plot.fit.(which_spec).maxPlot(kk);
    end
else
    % Determine a positive stagger to offset data, fit, residual, and 
    % baseline from the individual metabolite contributions
    stagData = 0.1*(max(abs(min(dataToPlot)), abs(max(dataToPlot))));
    maxPlot = max(dataToPlot + abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData;
end
% Add the data and plot
hold on;
plot(ppm, (zeros(1,length(ppm)) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground); % Zeroline
plot(ppm, (dataToPlot + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground); % Data
plot(ppm, (fit + stagData)/maxPlot, 'Color', MRSCont.colormap.Accent, 'LineWidth', 1.6); % Fit
plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Maximum Data

plot(ppm, (residual + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Residual
plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineStyle','--', 'LineWidth', 0.5); % Zeroline Residue
plot(ppm, (zeros(1,length(ppm)) + max(dataToPlot +  abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData)/maxPlot, 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Max Residue 

if (strcmp(which_spec, 'mm'))
   plot(ppm, (dataToPlot + stagData-Met_corr_spectrum)/maxPlot, 'Color',[1 0 0.1]); % Data
end

text(fitRangePPM(1), (0 + stagData)/maxPlot, '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline text
text(fitRangePPM(1), (0 + max(dataToPlot) + stagData)/maxPlot-0.05, num2str(max(dataToPlot),'%1.2e'), 'FontSize', 10,'Color',MRSCont.colormap.Foreground); % Maximum Data Text
text(fitRangePPM(1), (0 + max(dataToPlot +  abs(min(dataToPlot - fit))) + stagData)/maxPlot, '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline Residual text
text(fitRangePPM(1), (0 + max(dataToPlot +  abs(min(dataToPlot - fit))) + abs(max(dataToPlot - fit)) + stagData)/maxPlot +0.05, num2str(abs(max(dataToPlot - fit)),'%1.2e'), 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Max Residue text

if ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
    plot(ppm, (real(baseline) + stagData)/maxPlot, 'k', 'LineWidth', 1, 'Color', MRSCont.colormap.LightAccent);
end


%%% 6. PLOT BASIS FUNCTIONS %%%
% Plot separate metabolite basis functions only if not water
if ~(strcmp(which_spec, 'ref') || strcmp(which_spec, 'w'))
    if stagFlag
        % Staggered plots will be in all black and separated by the mean of the
        % maximum across all spectra
%         stag = max(abs(mean(max(real(appliedBasisSet.specs)))), abs(mean(min(real(appliedBasisSet.specs))))) * MRSCont.fit.scale{kk};
        stag = maxPlot *  2.5 / nBasisFct;
        % Loop over all basis functions
        
        for rr = 1:nBasisFct
            % Instead of a MATLAB legend, annotate each line separately with the
            % name of the metabolite
            plot(ppm, (indivPlots(:,rr) - rr*stag)/maxPlot, 'Color',MRSCont.colormap.Foreground);
            text(fitRangePPM(1), (- rr*stag)/maxPlot, basisSet.name{rr}, 'FontSize', 10,'Color',MRSCont.colormap.Foreground);
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
    set(gca, 'YLim', [0  1.2]);
    hold off;
    % If water is being shown, show a simple legend
    %legend('Data', 'Fit', 'Residual');
    
end


%%% 7. DESIGN FINETUNING %%%
% Adapt common style for all axes
set(gca, 'XDir', 'reverse', 'XLim', [fitRangePPM(1), fitRangePPM(end)], 'XMinorTick', 'On');
set(gca, 'LineWidth', 1, 'TickDir', 'out');
set(gca, 'FontSize', 16);
% If no y caption, remove y axis
if isempty(ylab)
    if ~MRSCont.flags.isGUI
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


%%% 8. ADD OSPREY LOGO AND TIGHTEN FIGURE %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
else
    out = tightfig(out);
end

end

   