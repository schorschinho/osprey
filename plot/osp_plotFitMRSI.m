function out = osp_plotFitMRSI(MRSCont, kk, which_spec,conc,mask)
%% out = osp_plotFitMRSI(MRSCont, kk, which, ppmmin, ppmmax)
%   Creates a figure showing processed data stored in an Osprey data container,
%   ie in the raw fields. This function will display the *processed and
%   averaged* data, i.e. after spectral alignment, averaging, water removal,
%   and other processing steps carried out in OspreyProcess.
%
%   USAGE:
%       out = osp_plotProcess(MRSCont, kk, which, ppmmin, ppmmax, xlab, ylab)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       which    = String for the spectrum to fit (optional)
%                   OPTIONS:    'A' (default)
%                               'B' (for MEGA, HERMES, HERCULES)
%                               'C' (for HERMES, HERCULES)
%                               'D' (for HERMES, HERCULES)
%                               'diff1' (for MEGA, HERMES, HERCULES)
%                               'diff2' (for HERMES, HERCULES)
%                               'sum' (for MEGA, HERMES, HERCULES)
%                               'ref'
%                               'w'
%       VoxelIndex = Index for the Voxel
%       xlab     = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab     = label for the y-axis (optional.  Default = '');
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

if nargin<5
    mask = 1;
    if nargin<4
        conc = 'diff1'; 
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



figTitle = ['MRSI fits: ' which_spec ];
% Set up colormaps
if isfield(MRSCont,'colormap')
    colormap = MRSCont.colormap;
    tintFactor = 0.75;
    colormap.ForegroundTint = [colormap.Foreground(1)+(1-colormap.Foreground(1))*tintFactor...
                               colormap.Foreground(2)+(1-colormap.Foreground(2))*tintFactor...
                               colormap.Foreground(3)+(1-colormap.Foreground(3))*tintFactor ];
else
    colormap.Background     = [1 1 1];
    colormap.LightAccent    = [110/255 136/255 164/255];
    colormap.Foreground     = [0 0 0];
    colormap.Accent         = [11/255 71/255 111/255];
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract raw and processed spectra in the plot range

if  strcmp(which_spec, 'conc')
    dataToPlot=op_takeVoxel(MRSCont.processed.(conc){kk},[1 1]);
else
    if strcmp(which_spec, 'off')
        dataToPlot=op_takeVoxel(MRSCont.processed.A{kk},[1 1]);
    else
        dataToPlot=op_takeVoxel(MRSCont.processed.(which_spec){kk},[1 1]);
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




% Get the fit parameters


fitParams   = MRSCont.fit.results{1, 1}.(which_spec).fitParams{kk};

% Pack up into structs to feed into the reconstruction functions
inputData.dataToFit                 = dataToPlot;
inputData.basisSet                  = basisSet;

inputSettings.scale                 = MRSCont.fit.scale{kk};

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

XVox = MRSCont.raw{kk}.nXvoxels;
YVox = MRSCont.raw{kk}.nYvoxels;  
npoints = length(ModelOutput.data);
procDataMarixToPlot = zeros((npoints+50) * XVox,YVox);
FitMarixToPlot = zeros((npoints+50) * XVox,YVox);
ResMarixToPlot = zeros((npoints+50) * XVox,YVox);

ppmLineToPlot = [];
for x = 1 : XVox
    ppmLineToPlot = horzcat(ppmLineToPlot, ModelOutput.ppm,ones(1,50)*nan);
end


for y = 1 : YVox
    procDataLineToPlot = [];
    FitLineToPlot = [];
    ResLineToPlot = [];
    for x = 1 : XVox
        if  strcmp(which_spec, 'conc')
            dataToPlot=op_takeVoxel(MRSCont.processed.(conc){kk},[x y]);
        else
            if strcmp(which_spec, 'off')
                dataToPlot=op_takeVoxel(MRSCont.processed.A{kk},[x y]);
            else
                dataToPlot=op_takeVoxel(MRSCont.processed.(which_spec){kk},[x y]);
            end
        end
         inputData.dataToFit                 = dataToPlot;
         fitParams   = MRSCont.fit.results{x, y}.(which_spec).fitParams{kk};
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
         procDataLineToPlot = vertcat(procDataLineToPlot,ModelOutput.data,ones(50,1)*nan);
         FitLineToPlot = vertcat(FitLineToPlot,ModelOutput.completeFit,ones(50,1)*nan);
         ResLineToPlot = vertcat(ResLineToPlot,ModelOutput.residual,ones(50,1)*nan);
    end
    procDataMarixToPlot(:,y) = procDataLineToPlot;
    FitMarixToPlot(:,y) = FitLineToPlot;
    ResMarixToPlot(:,y) = ResLineToPlot;
end



%%% 3. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
if ~MRSCont.flags.isGUI
    out = figure;
else
    out = figure('Visible','off');
end

shift = 0.075;
for y = 1 : YVox
    for x = 1 : XVox
        if mask
            if MRSCont.mask{kk}(x,y)
                plot((x-1)*(npoints+50)+1:x*(npoints+50),procDataMarixToPlot((x-1)*(npoints+50)+1:x*(npoints+50),y)+shift*y,'Color',colormap.Foreground);
                plot((x-1)*(npoints+50)+1:x*(npoints+50),FitMarixToPlot((x-1)*(npoints+50)+1:x*(npoints+50),y)+shift*y,'Color',colormap.Accent);
                plot((x-1)*(npoints+50)+1:x*(npoints+50),ResMarixToPlot((x-1)*(npoints+50)+1:x*(npoints+50),y)+shift*y+shift/2,'Color',colormap.Foreground);
                hold on
                text((x-1)*(npoints+50)+1, procDataMarixToPlot((x-1)*(npoints+50)+1,y)+shift*y-shift/10, num2str(fitRangePPM(1)), 'Color', colormap.ForegroundTint);
                text(x*(npoints+50)-50, procDataMarixToPlot((x-1)*(npoints+50)+1,y)+shift*y-shift/10, num2str(fitRangePPM(2)), 'Color', colormap.ForegroundTint);
            end
        else
            plot((x-1)*(npoints+50)+1:x*(npoints+50),procDataMarixToPlot((x-1)*(npoints+50)+1:x*(npoints+50),y)+shift*y,'Color',colormap.Foreground);
            plot((x-1)*(npoints+50)+1:x*(npoints+50),procDataMarixToPlot((x-1)*(npoints+50)+1:x*(npoints+50),y)+shift*y,'Color',colormap.Foreground);
            plot((x-1)*(npoints+50)+1:x*(npoints+50),ResMarixToPlot((x-1)*(npoints+50)+1:x*(npoints+50),y)+shift*y+shift/10,'Color',colormap.Foreground);
            text((x-1)*(npoints+50)+1, procDataMarixToPlot((x-1)*(npoints+50)+1,y)+shift*y-shift/10, num2str(fitRangePPM(1)), 'Color', colormap.ForegroundTint);
                text(x*(npoints+50)-50, procDataMarixToPlot((x-1)*(npoints+50)+1,y)+shift*y-shift/10, num2str(fitRangePPM(2)), 'Color', colormap.ForegroundTint);
            hold on
        end
    end
end

hold off



%%% 8. DESIGN FINETUNING %%%
% Adapt common style for all axes

    

gcf = out;
set(gcf, 'Color', MRSCont.colormap.Background);        
box off;
set(gca, 'XDir', 'reverse');


if ~MRSCont.flags.isGUI
    set(gca, 'YColor', 'w');
    % Black axes, white background
    set(gca, 'XColor', 'w');
    set(gca, 'Color', 'w');
    set(gcf, 'Color', 'w');
    title(figTitle, 'Interpreter', 'none');
else
    set(gca, 'YColor', MRSCont.colormap.Background);
    set(gca,'YTickLabel',{});
    set(gca,'YTick',{});
    set(gca,'XTickLabel',{});
    set(gca,'XTick',{});
    set(gca, 'XColor', MRSCont.colormap.Background);
    set(gca, 'Color', MRSCont.colormap.Background);
    set(gcf, 'Color', MRSCont.colormap.Background);
    title(figTitle, 'Interpreter', 'none', 'Color', MRSCont.colormap.Foreground);
end

%%% 9. ADD OSPREY LOGO %%%
% Add to the printout, but not if displayed in the GUI.
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    text(gca, 0, -0.1, [MRSCont.ver.Osp],'Color', colormap.Foreground);
    imshow(I, map);
    axis off;
end
end

   