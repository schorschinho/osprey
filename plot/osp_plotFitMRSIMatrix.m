function out = osp_plotFitMRSIMatrix(MRSCont, kk, which_spec,conc,coord, ppmmin, ppmmax,which_slice,yLim,mask)
%% out = osp_plotProcessMRSI(MRSCont, kk, which, ppmmin, ppmmax)
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

% Check that OspreyProcess has been run before
if ~MRSCont.flags.didProcess
    error('Trying to plot processed data, but data has not been processed yet. Run OspreyProcess first.')
end
out= figure;
%%% 1. PARSE INPUT ARGUMENTS %%%
% Get the fit method and style
fitMethod   = MRSCont.opts.fit.method;
fitStyle    = MRSCont.opts.fit.style;
% Fall back to defaults if not provided
% if nargin <10
%     yLim = [];
%     if nargin < 9
%         which_slice =1;
%         if nargin < 8
%         add_No_MoCo = 0;
%             if nargin < 7
%             lb = 0;
%                 if nargin<6
%                 switch which_spec{1}
%                 case {'A', 'B', 'C', 'D', 'diff1', 'diff2','diff3', 'sum','mm'}
%                 ppmmax = 5;
%                 case {'ref', 'w'}
%                 ppmmax = 2*4.68;
%                 otherwise
%                 error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
%                 end
%                     if nargin<5
%                     switch which_spec{1}
%                     case {'A', 'B', 'C', 'D', 'diff1', 'diff2','diff3', 'sum'}
%                     ppmmin = 0.2;
%                     case {'ref', 'w','mm'}
%                     ppmmin = 0;
%                     otherwise
%                     error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
%                     end
%                         if nargin < 4
%                         coord = [2 8; 2 8];
%                             if nargin < 3
%                             which_spec = 'A';
%                                 if nargin < 2
%                                 kk = 1;
%                                     if nargin<1
%                                         error('ERROR: no input Osprey container specified.  Aborting!!');
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end   
%         end
%     end
% end
% Set up colormaps
if isfield(MRSCont,'colormap')
    colormap = MRSCont.colormap;
    tintFactor = 0.50;
    colormap.ForegroundTint{1} = [colormap.Foreground(1)+(1-colormap.Foreground(1))*tintFactor...
                               colormap.Foreground(2)+(1-colormap.Foreground(2))*tintFactor...
                               colormap.Foreground(3)+(1-colormap.Foreground(3))*tintFactor ];
                           
    tintFactor = 0.30;
    colormap.ForegroundTint{2} = [colormap.Foreground(1)+(1-colormap.Foreground(1))*tintFactor...
                               colormap.Foreground(2)+(1-colormap.Foreground(2))*tintFactor...
                               colormap.Foreground(3)+(1-colormap.Foreground(3))*tintFactor ];
   
else
    colormap.Background     = [1 1 1];
    colormap.LightAccent    = [110/255 136/255 164/255];
    colormap.Foreground     = [0 0 0];
    colormap.Accent         = [11/255 71/255 111/255];
end
mask = ((mask(:,:,which_slice)));

%%% 2. EXTRACT DATA TO PLOT %%%
% Extract raw and processed spectra in the plot range


% if MRSCont.flags.hasRef || MRSCont.flags.hasWater
%     if MRSCont.flags.hasRef
%         Norm = MRSCont.plot.processed.ref.ContMax;
%     else
%         Norm = MRSCont.plot.processed.w.ContMax;
%     end
%     if isempty(yLim)
%         yLim = [min(MRSCont.plot.processed.(which_spec{1}).min) max(MRSCont.plot.processed.(which_spec{1}).max*2)]/Norm;
%         if strcmp(which_spec{1},'diff1') || strcmp(which_spec{1},'diff2')
%             yLim = [min(MRSCont.plot.processed.(which_spec{1}).min)/2 max(MRSCont.plot.processed.(which_spec{1}).max)]/Norm;
%         end
%     end
% else
%     Norm = MRSCont.plot.processed.A.ContMax;
%     if isempty(yLim)
%         yLim = [min(MRSCont.plot.processed.A.min) max(MRSCont.plot.processed.A.max)]/Norm;
%     end
% end
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
    basisSet    = MRSCont.fit.resBasisSet.(which_spec).water{1};
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



%%% 3. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized

plot_index = 1;

t=tiledlayout(coord(2,2)-coord(2,1)+1,coord(1,2)-coord(1,1)+1,'TileSpacing','compact');


for y = coord(2,1) : coord(2,2)
    for x = coord(1,1) : coord(1,2)
%         out = subaxis(coord(2,2)-coord(2,1)+1,coord(1,2)-coord(1,1)+1,plot_index, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0, 'MarginLeft', 0.075, 'MarginBottom', 0.1);
        nexttile
        hold on
        
        if  strcmp(which_spec, 'conc')
            dataToPlot=op_takeVoxel(MRSCont.processed.(conc){kk},[x y which_slice]);
        else
            if strcmp(which_spec, 'off')
                dataToPlot=op_takeVoxel(MRSCont.processed.A{kk},[x y which_slice]);
            else
                dataToPlot=op_takeVoxel(MRSCont.processed.(which_spec){kk},[x y which_slice]);
            end
        end
         inputData.dataToFit                 = dataToPlot;
         fitParams   = MRSCont.fit.results{x, y, which_slice}.(which_spec).fitParams{kk};
%          fitParams.LM_out.iteration
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
         
        plot(ModelOutput.ppm,ModelOutput.data,'Color',colormap.Foreground);    
        plot(ModelOutput.ppm,ModelOutput.completeFit,'Color',colormap.Accent);
        plot(ModelOutput.ppm,ModelOutput.residual + yLim(2)*1.1,'Color',colormap.ForegroundTint{1});
%          text(ppmmax,yLim(2),[num2str(x) ' ' num2str(y)]);
        plot_index = plot_index +1;    
        box off;
        set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', [yLim(1) yLim(2) * 1.2]);
        set(gca, 'LineWidth', 1, 'TickDir', 'in', 'XMinorTick', 'Off','TickLength', [0.03 0.025]);
        set(gca, 'XColor', colormap.Foreground, 'YColor', colormap.Foreground);

        if ~mask(x,y)
            set(gca, 'XColor', colormap.Background, 'YColor', colormap.Background);
        end
        
        if x ~= coord(1,1)
            set(gca,'YTickLabel',{})
            else
            text(ppmmax,yLim(2),[num2str(y)]);
        end
        if y == coord(2,1)
            text(ppmmax,yLim(2),[num2str(x)]);
        end
        if y ~= coord(2,2)
            set(gca,'XTickLabel',{})
        else
            xlabel('Frequency (ppm)', 'Color', colormap.Foreground)
        end
    end
end

set(gcf, 'Color', MRSCont.colormap.Background);   

end

   