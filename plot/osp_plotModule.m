function out = osp_plotModule(MRSCont, Module, kk,Index, which, metab, corr)
%% osp_plotModule
%   Callback function on print figure button click.
%
%
%   USAGE:
%       osp_plotModule(MRSCont, Module, kk, which, metab, corr)
%
%   INPUT:     MRSCont  = Osprey data container.
%              Module       = String for the Module
%              OPTIONS:    - 'OspreyLoad' (default)
%                          - 'OspreyProcess'
%                          - 'OspreyFit'
%                          - 'OspreyCoreg'
%                          - 'OspreySeg'
%                          - 'OspreySpecOverview'
%                          - 'OspreyMeanOverview'
%                          - 'OspreyRaincloudOverview'
%                          - 'OspreyScatterOverview'
%              kk       = Index for the kk-th dataset (optional. Default = 1)
%              which    = String for the plot this can either indicate the
%                           spectrum or the fit style/model
%              OPTIONS: 'mets' (OspreyLoad)
%                       'A' (OspreyProcess,OspreySpecOverview,OspreyMeanOverview)
%                       'B' (OspreyProcess,OspreySpecOverview,OspreyMeanOverview)
%                       'C' (OspreyProcess,OspreySpecOverview,OspreyMeanOverview)
%                       'D' (OspreyProcess,OspreySpecOverview,OspreyMeanOverview)
%                       'off' (OspreyFit)
%                       'diff1' (OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'diff2' (OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'sum' (OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'ref' (OspreyLoad,OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'w' (OspreyLoad,OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       %%%%%%%%%%%%%%% For Unedited %%%%%%%%%%%%%%%%%%%%%%%%
%                       'off-tCr' 'off-rawWaterScaled' 'off-CSFWaterScaled' 'off-TissCorrWaterScaled' (OspreyRaincloudOverview,OspreyScatterOverview)
%                        %%%%%%%%%%%%%%% Seperate MEGA-PRESS %%%%%%%%%%%%%%%%%%%%%%%%
%                       'A-tCr' 'A-rawWaterScaled' 'A-CSFWaterScaled' 'A-TissCorrWaterScaled' (OspreyRaincloudOverview,OspreyScatterOverview)
%                        %%%%%%%%%%%%%%% Seperate MEGA-PRESS/HERMES/HERCULES %%%%%%%%%%%%%%%%%%%%%%%%
%                       'diff1-tCr' 'diff1-rawWaterScaled' 'diff1-CSFWaterScaled' 'diff1-TissCorrWaterScaled' (OspreyRaincloudOverview,OspreyScatterOverview)
%                        %%%%%%%%%%%%%%% Seperate HERMES/HERCULES %%%%%%%%%%%%%%%%%%%%%%%%
%                       'diff2-tCr' 'diff2-rawWaterScaled' 'diff2-CSFWaterScaled' 'diff2-TissCorrWaterScaled' (OspreyRaincloudOverview,OspreyScatterOverview)
%                        %%%%%%%%%%%%%%% Seperate HERMES/HERCULES %%%%%%%%%%%%%%%%%%%%%%%%
%                       'sum-tCr' 'sum-rawWaterScaled' 'sum-CSFWaterScaled' 'sum-TissCorrWaterScaled' (OspreyRaincloudOverview,OspreyScatterOverview)
%                        %%%%%%%%%%%%%%% Concatenated MEGA-PRESS/HERMES/HERCULES %%%%%%%%%%%%%%%%%%%%%%%%
%                       'conc-tCr' 'conc-rawWaterScaled' 'conc-CSFWaterScaled' 'conc-TissCorrWaterScaled' (OspreyRaincloudOverview,OspreyScatterOverview)
%              metab    = String for the metab e.g. (GABA)
%              corr     = String for the correlation
%              OPTIONS: 'corr-X' (OspreyScatterOverview) with X indicating the correlation measure X
%                       'metab-X' (OspreyScatterOverview) with X indicating the metab X
%                       'SNR' (OspreyScatterOverview) correlation with SNR
%                       'FWHM' (OspreyScatterOverview) correlation with FWHM
%             Index     Index to basis set or subspectrum.
%             OPTIONS:
%
%
%
%   OUTPUT:     figure handle
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-02-03)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-01-16: First version of the code.


%%% ----- 1. CREATE MAIN FIGURE & SETTINGS ----- %%%
% Set color theme
colormapfig.Background = [255/255 254/255 254/255];
colormapfig.LightAccent = [110/255 136/255 164/255];
colormapfig.Foreground = [11/255 71/255 111/255];
colormapfig.Accent = [11/255 71/255 111/255];
% Determine font accorind to os
font=osp_platform('fonts');
font = font.helvetica;


% Create main window with fixed size
screenSize      = get(0,'ScreenSize');
canvasSize      = screenSize;
canvasSize(4)   = screenSize(4) * 0.7;
canvasSize(3)   = canvasSize(4) * (11/8.5);
canvasSize(2)   = (screenSize(4) - canvasSize(4))/2;
canvasSize(1)   = (screenSize(3) - canvasSize(3))/2;
out = figure('NumberTitle', 'off', 'Visible', 'on', 'Menu', 'none', ...
             'Position', canvasSize, 'ToolBar', 'none', ...
             'HandleVisibility', 'off', 'Renderer', 'painters', ...
             'Color', colormapfig.Background, 'Tag','MainFigure');
MRSCont.flags.isGUI = 1;
Title = MRSCont.ver.Osp;

Frame = uix.Panel('Parent', out, 'Padding', 1, 'Title', Title, ...
                  'FontName', font, 'BackgroundColor', colormapfig.Background, ...
                  'ForegroundColor', colormapfig.Foreground, ...
                  'HighlightColor', colormapfig.Background, ...
                  'ShadowColor', colormapfig.Background);
input_figure = uix.VBox('Parent', Frame,  'BackgroundColor',colormapfig.Background, 'Spacing', 5);
box = uix.HBox('Parent', input_figure,'BackgroundColor',colormapfig.Background, 'Spacing',6);
Info = uix.Panel('Parent',box, 'Padding', 5, 'Title', MRSCont.files{1,kk},...
                 'FontName', font, 'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground,...
                 'HighlightColor', colormapfig.Foreground, 'ShadowColor', colormapfig.Foreground);
LogoFig = figure('Visible','off');
[I, map] = imread('osprey.gif','gif');
axes(LogoFig, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
imshow(I, map);
axis off;
ViewAxes = gca();
set(ViewAxes,'Box','off','XColor', 'none','YColor','none', 'Parent', box);
set(box, 'Width', [-0.9 -0.1]);

% Clean up sequence name if needed
if ~strcmp('',MRSCont.raw{1,kk}.seq)
    if strcmp(sprintf('\n'),MRSCont.raw{1,kk}.seq(end))
        Seq = MRSCont.raw{1,kk}.seq(1:end-1);
    else
        Seq = MRSCont.raw{1,kk}.seq;
    end
else
    Seq = MRSCont.raw{1,kk}.seq;
end

Geom = fieldnames(MRSCont.raw{1,1}.geometry.size);


%%% ----- 2. CREATE INDIVIDUAL FIGURES FOR EACH ANALYSIS STAGE ----- %%%
% Depending on the Module input, the layout for the different analysis
% steps is defined below.

switch Module

    %%% --- 2a. OspreyLoad --- %%%
    case 'OspreyLoad'
        outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyLoad');
        [~,filename,~]  = fileparts(MRSCont.files{kk});

        % Grid for Plot and Data control sliders
        if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) % HBox for HERMES/HERCULES
            Plot = uix.HBox('Parent', input_figure, 'BackgroundColor',colormapfig.Background, 'Units', 'normalized');
        else
            Plot = uix.VBox('Parent', input_figure, 'BackgroundColor',colormapfig.Background, 'Units', 'normalized');
        end

        % Information text panel
        InfoText  = uicontrol('Parent',Info,'style','text',...
            'FontSize', 12, 'FontName', font,'ForegroundColor', colormapfig.Foreground,...
            'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',colormapfig.Background);

        % Get some information about the data from MRSCont to fill the info panel
        switch which
                case 'metabolites'
                StatText = ['Metabolite Data -> Sequence: ' Seq '; B0: ' num2str(MRSCont.raw{1,kk}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,kk}.te) ' / ' num2str(MRSCont.raw{1,kk}.tr) '\naverages: ' num2str(MRSCont.raw{1,kk}.averages)...
                             '; Sz: ' num2str(MRSCont.raw{1,kk}.sz) '; dimensions: ' num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{1})) ' x ' num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{2})) ' x ' num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{3})) ' mm = '...
                             num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{1}) * MRSCont.raw{1,kk}.geometry.size.(Geom{2}) * MRSCont.raw{1,kk}.geometry.size.(Geom{3})/1000) ' ml'];
                case 'MM'
                        StatText = ['MM Data -> Sequence: ' Seq '; B0: ' num2str(MRSCont.raw_mm{1,kk}.Bo) '; TE / TR: ' num2str(MRSCont.raw_mm{1,kk}.te) ' / ' num2str(MRSCont.raw_mm{1,kk}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_mm{1,kk}.spectralwidth) ' Hz'...   %re_mm
                             '\nraw subspecs: ' num2str(MRSCont.raw_mm{1,kk}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_mm{1,kk}.rawAverages) '; averages: ' num2str(MRSCont.raw_mm{1,kk}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_mm{1,kk}.sz) '; dimensions: ' num2str(MRSCont.raw_mm{1,kk}.geometry.size.(Geom{1})) ' x ' num2str(MRSCont.raw_mm{1,kk}.geometry.size.(Geom{2})) ' x ' num2str(MRSCont.raw_mm{1,kk}.geometry.size.(Geom{3})) ' mm = '...   %re_mm
                             num2str(MRSCont.raw_mm{1,kk}.geometry.size.(Geom{1}) * MRSCont.raw_mm{1,kk}.geometry.size.(Geom{2}) * MRSCont.raw_mm{1,kk}.geometry.size.(Geom{3})/1000) ' ml'];   %re_mm
                case 'ref'
                StatText = ['Reference Data -> Sequence: ' Seq '; B0: ' num2str(MRSCont.raw_ref{1,kk}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,kk}.te) ' / ' num2str(MRSCont.raw_ref{1,kk}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,kk}.spectralwidth) ' Hz'...   %re_mm
                             '\nraw subspecs: ' num2str(MRSCont.raw_ref{1,kk}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,kk}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,kk}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_ref{1,kk}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,kk}.geometry.size.(Geom{1})) ' x ' num2str(MRSCont.raw_ref{1,kk}.geometry.size.(Geom{2})) ' x ' num2str(MRSCont.raw_ref{1,kk}.geometry.size.(Geom{3})) ' mm = '...   %re_mm
                             num2str(MRSCont.raw_ref{1,kk}.geometry.size.(Geom{1}) * MRSCont.raw_ref{1,kk}.geometry.size.(Geom{2}) * MRSCont.raw_ref{1,kk}.geometry.size.(Geom{3})/1000) ' ml'];   %re_mm
               case 'MM reference'
                        StatText = ['MM reference Data -> Sequence: ' Seq '; B0: ' num2str(MRSCont.raw_mm_ref{1,kk}.Bo) '; TE / TR: ' num2str(MRSCont.raw_mm_ref{1,kk}.te) ' / ' num2str(MRSCont.raw_mm_ref{1,kk}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_mm_ref{1,kk}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw_mm_ref{1,kk}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_mm_ref{1,kk}.rawAverages) '; averages: ' num2str(MRSCont.raw_mm_ref{1,kk}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_mm_ref{1,kk}.sz) '; dimensions: ' num2str(MRSCont.raw_mm_ref{1,kk}.geometry.size.(Geom{1})) ' x ' num2str(MRSCont.raw_mm_ref{1,kk}.geometry.size.(Geom{2})) ' x ' num2str(MRSCont.raw_mm_ref{1,kk}.geometry.size.(Geom{3})) ' mm = '...
                             num2str(MRSCont.raw_mm_ref{1,kk}.geometry.size.(Geom{1}) * MRSCont.raw_mm_ref{1,kk}.geometry.size.(Geom{2}) * MRSCont.raw_mm_ref{1,kk}.geometry.size.(Geom{3})/1000) ' ml'];
                case 'w'
                    StatText = ['Water Data -> Sequence: ' Seq '; B0: ' num2str(MRSCont.raw_w{1,kk}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,kk}.te) ' / ' num2str(MRSCont.raw_w{1,kk}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,kk}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw_w{1,kk}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,kk}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,kk}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_w{1,kk}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,kk}.geometry.size.(Geom{1})) ' x ' num2str(MRSCont.raw_w{1,kk}.geometry.size.(Geom{2})) ' x ' num2str(MRSCont.raw_w{1,kk}.geometry.size.(Geom{3})) ' mm = '...
                             num2str(MRSCont.raw_w{1,kk}.geometry.size.(Geom{1}) * MRSCont.raw_w{1,kk}.geometry.size.(Geom{2}) * MRSCont.raw_w{1,kk}.geometry.size.(Geom{3})/1000) ' ml'];
        end
        set(InfoText, 'String', sprintf(StatText));

        %%% 2aa. VISUALIZATION PART OF OSPREYLOAD %%%
        % osp_plotLoad is used to visualize the raw data. Number of subplots
        % depends on the number of subspectra of the sequence
        Exp = Index(1);
        VoxelIndex = Index(2);
        switch which
            case 'metabolites'
            temp = osp_plotLoad(MRSCont, kk,'mets',Exp);
            outputFile      = [filename '_Voxel_' num2str(VoxelIndex) '_Exp_' num2str(Exp) '_OspreyLoad_metabolites.pdf'];

            if MRSCont.flags.isUnEdited % One window for UnEdited
                drawnow;
                set( temp.Children(1), 'Parent', Plot );
            end
            if MRSCont.flags.isMEGA %Two windows for MEGA
                set( temp.Children(2), 'Parent', Plot );
                set( temp.Children(1), 'Parent', Plot );
                set(Plot,'Heights', [-0.49 -0.49]);
                set(Plot.Children(2), 'Units', 'normalized')
                set(Plot.Children(2), 'OuterPosition', [0,0.5,1,0.5])
                set(Plot.Children(1), 'Units', 'normalized')
                set(Plot.Children(1), 'OuterPosition', [0,0,1,0.5])
            end
            if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) %Four windows for HERMES/HERCULES
                multiACload = uix.VBox('Parent', Plot, 'Padding', 5, 'BackgroundColor',colormapfig.Background);
                multiAload = uix.VBox('Parent', multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);
                multiCload = uix.VBox('Parent', multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);
                multiBDload = uix.VBox('Parent', Plot,'Padding', 5, 'BackgroundColor',colormapfig.Background);
                multiBload = uix.VBox('Parent', multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);
                multiDload = uix.VBox('Parent', multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);
                set( temp.Children(1), 'Parent', multiDload );
                set( temp.Children(1), 'Parent', multiCload );
                set( temp.Children(1), 'Parent', multiBload );
                set( temp.Children(1), 'Parent', multiAload );
                set(Plot,'Width', [-0.49 -0.49]);
                set(multiDload.Children(1), 'Units', 'normalized')
                set(multiDload.Children(1), 'OuterPosition', [0,0,1,1])
                set(multiCload.Children(1), 'Units', 'normalized')
                set(multiCload.Children(1), 'OuterPosition', [0,0,1,1])
                set(multiBload.Children(1), 'Units', 'normalized')
                set(multiBload.Children(1), 'OuterPosition', [0,0,1,1])
                set(multiAload.Children(1), 'Units', 'normalized')
                set(multiAload.Children(1), 'OuterPosition', [0,0,1,1])
            end

        case 'MM'

                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, kk,'mm',Exp);
                         VoxelIndex = 1;
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, kk,'mm',Exp,Index(2));
                         VoxelIndex = Index(1);
                          else
                            temp = osp_plotLoad(MRSCont, kk,'mm',Exp,[Index(2) Index(3) Index(4)]);
                end
                outputFile      = [filename '_Voxel_' num2str(VoxelIndex) '_Exp_' num2str(Exp) '_OspreyLoad_MM.pdf'];
                if MRSCont.flags.isUnEdited % One window for UnEdited
                    ViewAxes = gca();
                    set( ViewAxes, 'Parent', Plot );
                end
                if MRSCont.flags.isMEGA %Two windows for MEGA
                    set( temp.Children(2), 'Parent', Plot );
                    set( temp.Children(1), 'Parent', Plot );
                    set(Plot,'Heights', [-0.49 -0.49]);
                    set(Plot.Children(2), 'Units', 'normalized')
                    set(Plot.Children(2), 'OuterPosition', [0,0.5,1,0.5])
                    set(Plot.Children(1), 'Units', 'normalized')
                    set(Plot.Children(1), 'OuterPosition', [0,0,1,0.5])
                end
                if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) %Four windows for HERMES/HERCULES
                    multiACload = uix.VBox('Parent', Plot, 'Padding', 5, 'BackgroundColor',colormapfig.Background);
                        multiAload = uix.VBox('Parent', multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);
                        multiCload = uix.VBox('Parent', multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);
                    multiBDload = uix.VBox('Parent', Plot,'Padding', 5, 'BackgroundColor',colormapfig.Background);
                        multiBload = uix.VBox('Parent', multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);
                        multiDload = uix.VBox('Parent', multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);
                    set( temp.Children(1), 'Parent', multiDload );
                    set( temp.Children(1), 'Parent', multiCload );
                    set( temp.Children(1), 'Parent', multiBload );
                    set( temp.Children(1), 'Parent', multiAload );
                    set(Plot,'Width', [-0.49 -0.49]);
                    set(multiDload.Children(1), 'Units', 'normalized')
                    set(multiDload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(multiCload.Children(1), 'Units', 'normalized')
                    set(multiCload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(multiBload.Children(1), 'Units', 'normalized')
                    set(multiBload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(multiAload.Children(1), 'Units', 'normalized')
                    set(multiAload.Children(1), 'OuterPosition', [0,0,1,1])

                end
            case 'ref'

                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, kk,'ref',Exp);
                        VoxelIndex = 1;
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, kk,'ref',Exp,Index(2));
                        VoxelIndex = Index(1);
                      else
                        temp = osp_plotLoad(MRSCont, kk,'ref',Exp,[Index(2) Index(3) Index(4)]);
                end
                drawnow;
                ViewAxes = gca(); %re_mm
                set( ViewAxes, 'Parent', Plot );
                outputFile      = [filename '_Voxel_' num2str(VoxelIndex) '_Exp_' num2str(Exp) '_OspreyLoad_ref.pdf'];


             case 'w'
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                    temp = osp_plotLoad(MRSCont, kk,'w',Exp);
                    VoxelIndex = 1;
                elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    temp = osp_plotLoad(MRSCont, kk,'w',Exp,Index(2));
                    VoxelIndex = Index(1);
              else
                temp = osp_plotLoad(MRSCont, kk,'w',Exp,[Index(2) Index(3) Index(4)]);
                end
                drawnow;
                ViewAxes = gca(); %re_mm
                set( ViewAxes, 'Parent', Plot );
                outputFile      = [filename '_Voxel_' num2str(VoxelIndex)  '_Exp_' num2str(Exp) '_OspreyLoad_w.pdf'];
            end
            set(input_figure, 'Heights', [-0.1 -0.9]);
            % Get rid of the Load figure
            close( temp );


    %%% --- 2b. OspreyProcess --- %%%
    case 'OspreyProcess'
        outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyProcess');
        [~,filename,~]  = fileparts(MRSCont.files{kk});
        Exp = Index(1);
        SubSpec = Index(2);
        VoxelIndex = 1;
        Plot = uix.HBox('Parent', input_figure, ...
            'Padding', 5,'BackgroundColor', colormapfig.Background);
        set(input_figure, 'Heights', [-0.11 -0.89]);
         % Get parameter from file to fill the info panel
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                if strcmp(which,'metab')
                     StatText = ['SNR(' MRSCont.processed.metab{1,kk}.QC_names{SubSpec}  '): '  num2str(MRSCont.QM.SNR.(which)(Exp,kk,SubSpec)) '; FWHM (' MRSCont.processed.metab{1,kk}.QC_names{SubSpec} '): '...
                                num2str(MRSCont.QM.FWHM.(which)(Exp,kk,SubSpec)) ' / ' (num2str(MRSCont.QM.FWHM.(which)(Exp,kk,SubSpec)/MRSCont.processed.(which){kk}.txfrq(Exp)*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM.freqShift.(which)(Exp,kk,SubSpec)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.A(kk)*MRSCont.processed.(which){kk}.txfrq(Exp)/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.A(kk)*MRSCont.processed.(which){kk}.txfrq(Exp)/1e6) ' Hz'];
                else
                    StatText = ['SNR(' MRSCont.processed.(which){1,kk}.QC_names{SubSpec} '): ' num2str(MRSCont.QM.SNR.(which)(Exp,kk)) '; FWHM (' MRSCont.processed.(which){1,kk}.QC_names{SubSpec} '): '...
                            num2str(MRSCont.QM.FWHM.(which)(Exp,kk)) ' / ' (num2str(MRSCont.QM.FWHM.(which)(Exp,kk)/MRSCont.processed.(which){kk}.txfrq(Exp)*1e6))...
                            ' Hz / ppm'];
                end
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                if strcmp(which,'metab')
                     StatText = ['SNR(' MRSCont.processed.metab{1,kk}.QC_names{SubSpec}  '): '  num2str(MRSCont.QM.SNR.(which)(Exp,kk,SubSpec)) '; FWHM (' MRSCont.processed.metab{1,kk}.QC_names{SubSpec} '): '...
                                num2str(MRSCont.QM.FWHM.(which)(Exp,kk,SubSpec)) ' / ' (num2str(MRSCont.QM.FWHM.(which)(Exp,kk,SubSpec)/MRSCont.processed.(which){kk}.txfrq(Exp)*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM.freqShift.(which)(Exp,kk,SubSpec)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.A(kk)*MRSCont.processed.(which){kk}.txfrq(Exp)/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.A(kk)*MRSCont.processed.(which){kk}.txfrq(Exp)/1e6) ' Hz'];
                else
                    StatText = ['SNR(' MRSCont.processed.metab{1,kk}.QC_names{SubSpec} '): ' num2str(MRSCont.QM.SNR.(which)(Exp,kk)) '; FWHM (' MRSCont.processed.metab{1,kk}.QC_names{SubSpec} '): '...
                            num2str(MRSCont.QM.FWHM.(which)(Exp,kk)) ' / ' (num2str(MRSCont.QM.FWHM.(which)(Exp,kk)/MRSCont.processed.(which){kk}.txfrq(Exp)*1e6))...
                            ' Hz / ppm'];
                end
            end

            if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                StatText = ['Voxel ' num2str(Index(1)) ': ' StatText];
            end
            InfoText  = uicontrol('Parent',Info,'style','text','FontSize', 12, 'FontName', font,...
                                         'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                         'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);

 %%% 4. VISUALIZATION PART OF THIS TAB %%%
 %osp_plotProcess is used to visualize the processed spectra
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                temp = osp_plotProcess(MRSCont, kk,which,SubSpec,Exp); % Create figure
                VoxelIndex = 1;
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                temp = osp_plotProcess(MRSCont, kk,which,SubSpec,Exp,Index(1)); %Create figure
                VoxelIndex = Index(1);
            end
            %Subplots are distributed here
                proSpecs = uix.VBox('Parent', Plot, 'Padding', 5, 'BackgroundColor',colormapfig.Background);
                    proPre = uix.VBox('Parent', proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);
                    proPost = uix.VBox('Parent', proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);
                proOut = uix.VBox('Parent', Plot,'Padding', 5, 'BackgroundColor',colormapfig.Background);
                    proDrift = uix.VBox('Parent', proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);
                    proAlgn = uix.VBox('Parent', proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',colormapfig.Background);

            set( temp.Children(1), 'Parent', proDrift );
            set( temp.Children(1), 'Parent', proAlgn );
            set( temp.Children(1), 'Parent', proPost );
            set( temp.Children(1), 'Parent', proPre );
            close( temp );
%%% 5. DATA CONTROLS FOR THIS TAB %%%
            set(Plot,'Widths', [-0.49 -0.49]);
            set(proPre.Children(1), 'Units', 'normalized')
            set(proPre.Children(1), 'OuterPosition', [0,0,1,1])
            set(proPost.Children(1), 'Units', 'normalized')
            set(proPost.Children(1), 'OuterPosition', [0,0,1,1])
            set(proDrift.Children(1), 'Units', 'normalized')
            set(proDrift.Children(1), 'OuterPosition', [0,0,1,1])
            set(proAlgn.Children(1), 'Units', 'normalized')
            set(proAlgn.Children(1), 'OuterPosition', [0,0,1,1])

            outputFile      = [filename '_Voxel_' num2str(VoxelIndex) '_OspreyProcess_' which '_' MRSCont.processed.(which){kk}.names{SubSpec} '.pdf'];


     %%% --- 2c. OspreyFit --- %%%
    case 'OspreyFit'
        basis = Index(1);
        subspectrum = Index(2);
        % For this visualization, we will have to make a few
        % distinctions upfront since the modeling algorithms (LCModel
        % vs. Osprey) do not always return the same kinds of data, or they
        % return them in different formats.
        switch MRSCont.opts.fit.method
            case 'LCModel'
                % Number of metabolites and lipid/MM basis functions
                basisNames = MRSCont.fit.results.metab.fitParams{kk}.name;
                nLip    = sum(~cellfun(@isempty, strfind(basisNames, 'Lip')));
                nMM     = sum(~cellfun(@isempty, strfind(basisNames, 'MM')));
                nMMLip  = nLip + nMM;
                nMets   = length(basisNames) - nMMLip;
                nComb   = sum(~cellfun(@isempty, strfind(basisNames, '_')));
                % No info panel string for the water fit range
                waterFitRangeString = '';
                % Where are the metabolite names stored?
                basisSetNames = MRSCont.fit.results.(which).fitParams{kk}.name;
                subSpecName = 'A';
                % Smaller fonts for the results
                resultsFontSize = 6;
            case 'Osprey'
                % Additional info panel string for the water fit range
                waterFitRangeString = ['Fitting range: ' num2str(MRSCont.opts.fit.rangeWater(1)) ' to ' num2str(MRSCont.opts.fit.rangeWater(2)) ' ppm'];
                % Where are the metabolite names stored?
                if strcmp(which, 'ref') || strcmp(which, 'w')
                    basisSet = MRSCont.fit.resBasisSet.(which).(['np_sw_' num2str(round(MRSCont.processed.metab{kk}.sz(1))) '_' num2str(round(MRSCont.processed.metab{kk}.spectralwidth))]){1};
                    basisSetNames = basisSet.name;
                    subSpecName = which;
                else if strcmp(which, 'conc')
                        basisSet = MRSCont.fit.resBasisSet.(which).(['np_sw_' num2str(round(MRSCont.processed.metab{kk}.sz(1))) '_' num2str(round(MRSCont.processed.metab{kk}.spectralwidth))]){basis,1};
                        basisSetNames = basisSet.name;
                        subSpecName = basisSet.names{1};
                    else
                        basisSet = MRSCont.fit.resBasisSet.(which).(['np_sw_' num2str(round(MRSCont.processed.metab{kk}.sz(1))) '_' num2str(round(MRSCont.processed.metab{kk}.spectralwidth))]){basis,1,subspectrum};
                        basisSetNames = basisSet.name;
                        subSpecName = basisSet.names{1};
                    end
                end
                 % Number of metabolites and lipid/MM basis functions
                nMets   = MRSCont.fit.basisSet.nMets;
                nMMLip  = MRSCont.fit.basisSet.nMM;
                % Larger fonts for the results
                resultsFontSize = 11;
        end

        % Build output folder filename
        outputFolder    = fullfile(MRSCont.outputFolder, 'Figures', 'OspreyFit');
        [~,filename,~]  = fileparts(MRSCont.files{kk});

        Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',colormapfig.Background);
        set(input_figure, 'Heights', [-0.12 -0.88]);
        %%% 2.cc FILLING FITTED AMPLITUDE PANEL %%%
        % Creates the panel on the right side with the fitted amplitudes
        if  ~strcmp (MRSCont.opts.fit.style, 'Concatenated') ||  strcmp(which, 'ref') || strcmp(which, 'w') %Is not concateneted or is reference/water fit
                switch MRSCont.opts.fit.method
                    case 'LCModel'
                        if strcmp(which, 'ref') || strcmp(which, 'w')
                            RawAmpl = MRSCont.fit.results.(which).fitParams{1,kk}.h2oarea .* MRSCont.fit.scale{kk};
                        else
                            RawAmpl = MRSCont.fit.results.(which).fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk};
                            CRLB    = MRSCont.fit.results.(which).fitParams{1,kk}.CRLB;
                        end
                    case 'Osprey'
                        RawAmpl = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.ampl .* MRSCont.fit.scale{kk};
                end
            else %Is concatenated and not water/reference
                which = 'conc';
            end
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                ph0 = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.ph0;
                ph1 = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.ph1;
                if ~strcmp(which, 'ref') && ~strcmp(which, 'w')
                    refShift = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.refShift;
                    refFWHM = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.refFWHM;
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(which).fitParams{1,kk}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(which).fitParams{1,kk}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                end
                RawAmpl = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.ampl .* MRSCont.fit.scale{kk};
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                if ~strcmp(which, 'ref') && ~strcmp(which, 'w')
                    refShift = MRSCont.fit.results{1,Index(1)}.(which).fitParams{1,kk}.refShift;
                    refFWHM = MRSCont.fit.results{1,Index(1)}.(which).fitParams{1,kk}.refFWHM;
                    ph0 = MRSCont.fit.results{1,Index(1)}.(which).fitParams{1,kk}.ph0;
                    ph1 = MRSCont.fit.results{1,Index(1)}.(which).fitParams{1,kk}.ph1;
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(which).fitParams{1,kk}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(which).fitParams{1,kk}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                end
                RawAmpl = MRSCont.fit.results{1,Index(1)}.(which).fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk};
            else
               if ~strcmp(which, 'ref') && ~strcmp(which, 'w')
                    refShift = MRSCont.fit.results{Index(1),Index(2)}.(which).fitParams{1,kk}.refShift;
                    refFWHM = MRSCont.fit.results{Index(1),Index(2)}.(which).fitParams{1,kk}.refFWHM;
                    ph0 = MRSCont.fit.results{Index(1),Index(2)}.(which).fitParams{1,kk}.ph0;
                    ph1 = MRSCont.fit.results{Index(1),Index(2)}.(which).fitParams{1,kk}.ph1;
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(which).fitParams{1,kk}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(which).fitParams{1,kk}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                end
                RawAmpl = MRSCont.fit.results{1,Index(1)}.(which).fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk};
            end
            % Get parameter from file to fill the info panel
            if  ~strcmp (which, 'ref') && ~strcmp (which, 'w') %Metabolite data
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                    StatText = ['Metabolite Data -> Sequence: ' Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' which,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.range(1)) ' to ' num2str(MRSCont.opts.fit.range(2)) ' ppm; Baseline knot spacing: ' num2str(MRSCont.opts.fit.bLineKnotSpace) ' ppm; ph0: ' num2str(ph0,'%1.2f'),...
                        'deg; ph1: ' num2str(ph1,'%1.2f') 'deg; refShift: ' num2str(refShift,'%1.2f') ' Hz; refFWHM: ' num2str(refFWHM,'%1.2f')...
                        ' ppm\nNumber of metabolites: ' num2str(nMets) '; Number of MM/lipids: ' num2str(nMMLip) ...
                        ' scale: '  num2str(MRSCont.fit.scale{kk})];
                elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                    StatText = ['Metabolite Data -> Sequence: ' Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' which,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.range(1)) ' to ' num2str(MRSCont.opts.fit.range(2)) ' ppm; Baseline knot spacing: ' num2str(MRSCont.opts.fit.bLineKnotSpace) ' ppm; ph0: ' num2str(ph0,'%1.2f'),...
                        'deg; ph1: ' num2str(ph1,'%1.2f') 'deg; refShift: ' num2str(refShift,'%1.2f') ' Hz; refFWHM: ' num2str(refFWHM,'%1.2f')...
                        ' ppm\nNumber of metabolites: ' num2str(nMets) '; Number of MM/lipids: ' num2str(nMMLip) ...
                        ' scale: '  num2str(MRSCont.fit.scale{kk})];
                else
                    switch MRSCont.opts.fit.method
                    case 'Osprey'
                        iniph0 = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.prelimParams.ph0;
                        iniph1 = MRSCont.fit.results.(which).fitParams{basis,kk,subspectrum}.prelimParams.ph1;
                    case 'LCModel'
                        iniph0 = nan;
                        iniph1 = nan;
                    end
                    StatText = ['Metabolite Data -> Sequence: ' Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' which,...
                            '\nFitting range: ' num2str(MRSCont.opts.fit.range(1)) ' to ' num2str(MRSCont.opts.fit.range(2)) ' ppm; Baseline knot spacing: ' num2str(MRSCont.opts.fit.bLineKnotSpace) ' ppm; ph0: ' num2str(ph0,'%1.2f'),...
                            'deg; ph1: ' num2str(ph1,'%1.2f') 'deg; refShift: ' num2str(refShift,'%1.2f') ' Hz; refFWHM: ' num2str(refFWHM,'%1.2f')...
                            ' ppm\nNumber of metabolites: ' num2str(nMets) '; Number of MM/lipids: ' num2str(nMMLip) ...
                            ' scale: '  num2str(MRSCont.fit.scale{kk})];
                end

            else if strcmp (which, 'ref') %Reference data?
            StatText = ['Reference Data -> Sequence: ' Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' which,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.rangeWater(1)) ' to ' num2str(MRSCont.opts.fit.rangeWater(2)) ' ppm'];
                else %Is water data
                    StatText = ['Water Data -> Sequence: ' Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' which,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.rangeWater(1)) ' to ' num2str(MRSCont.opts.fit.rangeWater(2)) ' ppm'];
                end
            end
    %%% 4. FILLING FITTED AMPLITUDE PANEL %%%
    % Creates the panel on the right side with the fitted ammplitudes
            InfoText  = uicontrol('Parent',Info,'style','text',...
                                        'FontSize', 12, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
            Results = uix.Panel('Parent', Plot,...
                                       'Title', ['Raw Amplitudes'],'FontName', font,'HighlightColor', colormapfig.Foreground,...
                                       'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground, 'ShadowColor', colormapfig.Foreground);
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                    if ~(strcmp(which, 'ref') || strcmp(which, 'w')) %Metabolite fit
                        NameText = [''];
                        RawAmplText = [''];
                        CRLBText    = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [basisSetNames{m} ' \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                            if strcmp(MRSCont.opts.fit.method, 'LCModel')
                                CRLBText = [CRLBText, [num2str(CRLB(m), '%i') '%%\n']];
                            end
                        end
                    else %Water/reference fit but this should never happen in this loop
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                    end
                    set(Results, 'Title', ['Raw Amplitudes']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                    if ~(strcmp(which, 'ref') || strcmp(which, 'w')) %Metabolite fit
                        switch MRSCont.opts.fit.method
                            case 'Osprey'
                                if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                                    RawAmpl = RawAmpl ./ (MRSCont.fit.results.ref.fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk});
                                else
                                    RawAmpl = RawAmpl ./ (MRSCont.fit.results.water.fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk});
                                end
                            case 'LCModel'
                        end
                        NameText = [''];
                        RawAmplText = [''];
                        CRLBText    = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [basisSetNames{m} ' \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                            if strcmp(MRSCont.opts.fit.method, 'LCModel')
                                CRLBText = [CRLBText, [num2str(CRLB(m), '%i') '%%\n']];
                            end
                        end
                        set(Results, 'Title', ['Raw Water Ratio']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                    else %Water/reference fit
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                       set(Results, 'Title', ['Raw Amplitudes']);
                       FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                       FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                       'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                       'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                       FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                       'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                       'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                    end
                end
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                    if ~(strcmp(which, 'ref') || strcmp(which, 'w')) %Metabolite fit
                        NameText = [''];
                        RawAmplText = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [MRSCont.fit.resBasisSet{Index(1),Index(2)}.(which).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        end
                    else %Water/reference fit but this should never happen in this loop
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                    end
                    set(Results, 'Title', ['Raw Amplitudes']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                    if ~(strcmp(which, 'ref') || strcmp(which, 'w')) %Metabolite fit
                        if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results{1,Index(1)}.ref.fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk});
                        else
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results{1,Index(1)}.w.fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk});
                        end
                        NameText = [''];
                        RawAmplText = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [MRSCont.fit.resBasisSet{1,Index(1)}.(which).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        end
                        set(Results, 'Title', ['Raw Water Ratio']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                    else %Water/reference fit
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                       set(Results, 'Title', ['Raw Amplitudes']);
                       FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                       FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                       'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                       'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                       FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                       'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                       'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                    end
                end
            else
                if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                    if ~(strcmp(which, 'ref') || strcmp(which, 'w')) %Metabolite fit
                        NameText = [''];
                        RawAmplText = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [MRSCont.fit.resBasisSet{Index(1),Index(2)}.(which).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        end
                    else %Water/reference fit but this should never happen in this loop
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                    end
                    set(Results, 'Title', ['Raw Amplitudes']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                    if ~(strcmp(which, 'ref') || strcmp(which, 'w')) %Metabolite fit
                        if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results{1,Index(1)}.ref.fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk});
                        else
                            RawAmpl = RawAmpl ./ (MRSCont.fit.results{1,Index(1)}.w.fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk});
                        end
                        NameText = [''];
                        RawAmplText = [''];
                        for m = 1 : length(RawAmpl) %Names and Amplitudes
                            NameText = [NameText, [MRSCont.fit.resBasisSet{1,Index(1)}.(which).(MRSCont.info.A.unique_ndatapoint_spectralwidth{1}).name{m} ': \n']];
                            RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                        end
                        set(Results, 'Title', ['Raw Water Ratio']);
                        FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                        FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                        FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                        'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                    else %Water/reference fit
                       NameText = ['Water: ' ];
                       RawAmplText = [num2str(RawAmpl,'%1.2e')];
                       set(Results, 'Title', ['Raw Amplitudes']);
                       FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                       FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                       'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                       'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                       FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                       'FontSize', 11, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                       'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                    end
                end
            end

        %%%  5. VISUALIZATION PART OF THIS TAB %%%
        %osp_plotFit is used to visualize the fits (off,diff1,diff2,sum,ref,water)
        temp = figure( 'Visible', 'off' );
        if  ~strcmp (MRSCont.opts.fit.style, 'Concatenated') ||  strcmp(which, 'ref') || strcmp(which, 'w') %Is not concateneted or is reference/water fit
            temp = osp_plotFit(MRSCont, kk, which,[1 subspectrum basis],which);
        end
        ViewAxes = gca();
        set(ViewAxes, 'Parent', Plot );
        close( temp );

        set(Plot,'Widths', [-0.16 -0.84]);
        set(Plot.Children(2), 'Units', 'normalized');
        set(Plot.Children(2), 'OuterPosition', [0.17,0.02,0.75,0.98])
        outputFile      = [filename '_OspreyFit_' which '_' subSpecName '_basis_' num2str(basis) '.pdf'];


    case {'OspreyCoreg','OspreySeg'} %Coreg/Seg
        outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyCoregSeg');
        [~,filename,~]  = fileparts(MRSCont.files{kk});

        % Creates layout for plotting and data control
        Plot = uix.HBox('Parent', input_figure,'BackgroundColor',colormapfig.Background);
        set(input_figure, 'Heights', [-0.1 -0.9]);
        % Get parameter from file to fill the info panel

        StatText = ['Metabolite Data -> Sequence: ' Seq '; B0: ' num2str(MRSCont.raw{1,kk}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,kk}.te) ' / ' num2str(MRSCont.raw{1,kk}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,kk}.spectralwidth) ' Hz'...
            '\nraw subspecs: ' num2str(MRSCont.raw{1,kk}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,kk}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,kk}.averages)...
            '; Sz: ' num2str(MRSCont.raw{1,kk}.sz) '; dimensions: ' num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{1})) ' x ' num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{2})) ' x ' num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{3})) ' mm = '...
            num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{1}) * MRSCont.raw{1,kk}.geometry.size.(Geom{2}) * MRSCont.raw{1,kk}.geometry.size.(Geom{3})/1000) ' ml'];

        InfoText  = uicontrol('Parent',Info,'style','text',...
            'FontSize', 12, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);

        %%% 2.VISUALIZATION PART OF THIS TAB %%%
        % In this case osp_plotCoreg or osp_plotSegment is used to visualize the
        % coregistration or the segmentation
        Results = uix.VBox('Parent', Plot,'BackgroundColor',colormapfig.Background);
        temp = figure( 'Visible', 'off' );
        if MRSCont.flags.didSeg %Did segment. In this case coreg has already been performed. Visualize both
            osp_plotCoreg(MRSCont, kk);
            ViewAxes = gca();
            set(ViewAxes, 'Parent', Results );
            colormap(Results.Children,'gray');
            close( temp );
            temp = figure( 'Visible', 'off' );
            osp_plotSegment(MRSCont, kk);
            ViewAxes = gca();
            set(ViewAxes, 'Parent', Results );
            colormap(Results.Children(1),'gray');
            close( temp );
        else % Only coreg has been run
            osp_plotCoreg(MRSCont, kk);
            ViewAxes = gca();
            set(ViewAxes, 'Parent', Results );
            colormap(Results.Children,'gray');
            close( temp );
        end
        outputFile      = [filename '_OspreyCoregSeg.pdf'];
    case 'OspreySpecOverview' %SpecOverview
        set(Info,'Title', 'Descriptive Information');
        groupString = '';
        for g = 1 : MRSCont.overview.NoGroups
            groupString = [groupString MRSCont.overview.groupNames{g} ' with ' num2str(sum(MRSCont.overview.groups(:) == g)) ' subjects; '];
        end
        StatText = ['Sequence: ' Seq '; Number of subjects: ' num2str(MRSCont.nDatasets) '; Number of Groups: ' num2str(MRSCont.overview.NoGroups) '\n'...
            'Distribution: ' groupString '\n'];

        InfoText  = uicontrol('Parent',Info,'style','text',...
            'FontSize', 12, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
        Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',colormapfig.Background);
        set(input_figure, 'Heights', [-0.1 -0.9]);
        outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview','Individual');
        outputFile  = [which '.pdf'];
        for g = 1 :  MRSCont.overview.NoGroups %Loop over groups
            temp = osp_plotOverviewSpec(MRSCont, which, g, 0.1,'Frequency (ppm)','','',Index);
            if g == 1
                fig_hold = temp;
                set(fig_hold,'Tag','fig_hold');
            else
                ViewAxes=get(temp,'Children');
                copyobj(ViewAxes.Children, fig_hold.Children(1));
                h = findall(groot,'Type','figure');
                for ff = 1 : length(h)
                    if ~(strcmp(h(ff).Tag, 'Osprey') ||  strcmp(h(ff).Tag, 'TMWWaitbar') ||  strcmp(h(ff).Tag, 'fig_hold') ||  strcmp(h(ff).Tag, 'MainFigure'))
                        close(h(ff))
                    end
                end
            end
        end
        drawnow
        set(fig_hold.Children, 'Parent', Plot );
    case 'OspreyMeanOverview' %MeanOverview
        set(Info,'Title', 'Descriptive Information');
        groupString = '';
        for g = 1 : MRSCont.overview.NoGroups
            groupString = [groupString MRSCont.overview.groupNames{g} ' with ' num2str(sum(MRSCont.overview.groups(:) == g)) ' subjects; '];
        end
        StatText = ['Sequence: ' Seq '; Number of subjects: ' num2str(MRSCont.nDatasets) '; Number of Groups: ' num2str(MRSCont.overview.NoGroups) '\n'...
            'Distribution: ' groupString '\n'];

        InfoText  = uicontrol('Parent',Info,'style','text',...
            'FontSize', 12, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
        Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',colormapfig.Background);
        set(input_figure, 'Heights', [-0.1 -0.9]);
        outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Mean');
        outputFile  = [which '.pdf'];
        for g = 1 :  MRSCont.overview.NoGroups %Loop over groups
            if MRSCont.overview.NoGroups > 1
                temp = osp_plotMeanSpec(MRSCont, which,g,1,1/MRSCont.overview.NoGroups);
                if g == 1
                    fig_hold = temp;
                else
                    ax=get(temp,'Children');
                    lines=get(ax,'Children');
                    copyobj(lines, fig_hold.Children(1));
                    close(temp);
                end
            else
                fig_hold = osp_plotMeanSpec(MRSCont, which,g);
            end
        end
        drawnow
        set(fig_hold.Children,'Children',flipud(fig_hold.Children.Children));
        set(fig_hold.Children, 'Parent', Plot );
        close(fig_hold);


    case 'OspreyRaincloudOverview' %Raincloud plot
        set(Info,'Title', 'Descriptive Information');
        groupString = '';
        for g = 1 : MRSCont.overview.NoGroups
            groupString = [groupString MRSCont.overview.groupNames{g} ' with ' num2str(sum(MRSCont.overview.groups(:) == g)) ' subjects; '];
        end
        StatText = ['Sequence: ' Seq '; Number of subjects: ' num2str(MRSCont.nDatasets) '; Number of Groups: ' num2str(MRSCont.overview.NoGroups) '\n'...
            'Distribution: ' groupString '\n'];

        InfoText  = uicontrol('Parent',Info,'style','text',...
            'FontSize', 12, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
        Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',colormapfig.Background);
        set(input_figure, 'Heights', [-0.1 -0.9]);
        outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Raincloud');
        split_which = strsplit(which,'-');
        if strcmp(split_which{3},'AlphaCorrWaterScaled') || strcmp(split_which{3},'AlphaCorrWaterScaledGroupNormed')
            metab = 'GABA';
        end
        fig_hold = osp_plotRaincloud(MRSCont,split_which{2},split_which{3},metab,'Raincloud plot',0,Index(1),Index(2));
        delete( fig_hold.Children(1));
        set(fig_hold.Children,'Children',flipud(fig_hold.Children.Children));
        set( fig_hold.Children, 'Parent', Plot );
        close(fig_hold);
        outputFile  = [metab '_' split_which{1} '_' split_which{2} '_' split_which{3} '.pdf'];
    case 'OspreyScatterOverview' %Correlation plot
        set(Info,'Title', 'Descriptive Information');
        groupString = '';
        for g = 1 : MRSCont.overview.NoGroups
            groupString = [groupString MRSCont.overview.groupNames{g} ' with ' num2str(sum(MRSCont.overview.groups(:) == g)) ' subjects; '];
        end
        StatText = ['Sequence: ' Seq '; Number of subjects: ' num2str(MRSCont.nDatasets) '; Number of Groups: ' num2str(MRSCont.overview.NoGroups) '\n'...
            'Distribution: ' groupString '\n'];

        InfoText  = uicontrol('Parent',Info,'style','text',...
            'FontSize', 12, 'FontName', font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
        Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',colormapfig.Background);
        set(input_figure, 'Heights', [-0.1 -0.9]);
        outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Correlation');
        split_which = strsplit(which,'-');
        if strcmp(split_which{3},'AlphaCorrWaterScaled') || strcmp(split_which{3},'AlphaCorrWaterScaledGroupNormed')
            metab = 'GABA';
        end
        MRSCont.flags.isGUI = 0;
        split_corr = strsplit(corr,'-');
        if strcmp(split_corr{1},'corr')
            fig_hold = osp_plotScatter(MRSCont,split_which{2},split_which{3},metab,MRSCont.overview.corr.Meas{str2num(split_corr{2})},MRSCont.overview.corr.Names{str2num(split_corr{2})},1,Index(1));
            outputFile  = [metab '_' split_which{2} '_' split_which{3} '_basis' num2str(Index(1)) '_'  MRSCont.overview.corr.Names{str2num(split_corr{2})} '.pdf'];
        else if strcmp(split_corr{1},'metab')
                fig_hold = osp_plotScatter(MRSCont,split_which{1},split_which{2},metab,MRSCont.overview.corr.Names{str2num(split_corr{1})},1,Index(1),1,Index(1));
                outputFile  = [metab '_' split_which{2} '_' split_which{3} '_basis' num2str(Index(1)) '_'  metab '.pdf'];
            else if strcmp(split_corr{1},'SNR')
                    fig_hold = osp_plotScatter(MRSCont,split_which{2},split_which{3},metab,MRSCont.QM.SNR.metab(1,:,Index(1))','SNR Subspectrum A',1,Index(1));
                    outputFile  = [metab '_' split_which{2} '_' split_which{3} '_basis' num2str(Index(1)) '_SNR.pdf'];
                else if  strcmp(split_corr{1},'FWHM')
                        fig_hold = osp_plotScatter(MRSCont,split_which{2},split_which{3},metab,MRSCont.QM.FWHM.metab(1,:,Index(1))','FWHM (Hz)',1,Index(1));
                        outputFile  = [metab '_' split_which{2} '_' split_which{3} '_basis' num2str(Index(1)) '_FWHM.pdf'];
                    end
                end
            end
        end

        delete( fig_hold.Children(1));
        delete( fig_hold.Children(1));
        set(fig_hold.Children,'Children',flipud(fig_hold.Children.Children));
        set(fig_hold.Children, 'Parent', Plot );
        close(fig_hold);

end
set(out,'Renderer','painters','Menu','none','Toolbar','none');

%% Clean up and save
% Save the figure to the output folder
% Determine output folder
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
fig_pos = out.PaperPosition;
out.PaperSize = [fig_pos(3) fig_pos(4)];

% print(fig,'-dpdf','-painters','-r600','-bestfit',strcat(plot_path,plot_name));

% print(out,fullfile(outputFolder,outputFile),'-dpdf') % then print it
saveas(out,fullfile(outputFolder,outputFile),'pdf');
h = findall(groot,'Type','figure');
for ff = 1 : length(h)
    if ~(strcmp(h(ff).Tag, 'Osprey') ||  strcmp(h(ff).Tag, 'TMWWaitbar'))
        close(h(ff))
    end
end
end
