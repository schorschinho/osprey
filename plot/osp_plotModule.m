function out = osp_plotModule(MRSCont, Module, kk, which, metab, corr)
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
%              OPTIONS: 'A' (OspreyLoad,OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'B' (OspreyLoad,OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'C' (OspreyLoad,OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'D' (OspreyLoad,OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'diff1' (OspreyLoad,OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'diff2' (OspreyLoad,OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'sum' (OspreyLoad,OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'ref' (OspreyLoad,OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                       'w' (OspreyLoad,OspreyProcess,OspreyFit,OspreySpecOverview,OspreyMeanOverview)
%                        %%%%%%%%%%%%%%% For Unedited  or Seperate MEGA-PRESS %%%%%%%%%%%%%%%%%%%%%%%%        
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
%%% 1. GET DATA %%%
    colormapfig.Background = [255/255 254/255 254/255];
    colormapfig.LightAccent = [110/255 136/255 164/255];
    colormapfig.Foreground = [11/255 71/255 111/255];
    colormapfig.Accent = [11/255 71/255 111/255];
    
    screenSize      = get(0,'ScreenSize');
    canvasSize      = screenSize;
    canvasSize(4)   = screenSize(4) * 0.6;
    canvasSize(3)   = canvasSize(4) * (11/8.5);
    canvasSize(2)   = (screenSize(4) - canvasSize(4))/2;
    canvasSize(1)   = (screenSize(3) - canvasSize(3))/2;
    out = figure('NumberTitle', 'off', 'Visible', 'on', 'Menu', 'none','Position', canvasSize,...
                    'ToolBar', 'none', 'HandleVisibility', 'off', 'Renderer', 'painters', 'Color', colormapfig.Background);
    input_figure = uix.VBox('Parent', out,  'BackgroundColor',colormapfig.Background, 'Spacing', 5);                
    box = uix.HBox('Parent', input_figure,'BackgroundColor',colormapfig.Background, 'Spacing',6);
    Info = uix.Panel('Parent',box, 'Padding', 5, 'Title', MRSCont.files{kk},...
                                 'FontName', 'Arial', 'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground,...
                                 'HighlightColor', colormapfig.Foreground, 'ShadowColor', colormapfig.Foreground);
    LogoFig = figure('Visible','off'); 
    [I, map] = imread('osprey.gif','gif');
    axes(LogoFig, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
    ViewAxes = gca();
    set( ViewAxes,'Box','off','XColor', 'none','YColor','none', 'Parent', box );

    set(box, 'Width', [-0.9 -0.1]);
    if strcmp(sprintf('\n'),MRSCont.raw{1,kk}.seq(end)) %Clean up Sequence Name if needed
        Seq = MRSCont.raw{1,kk}.seq(1:end-1);
    else
        Seq = MRSCont.raw{1,kk}.seq;
    end
    Geom = fieldnames(MRSCont.raw{1,1}.geometry.size);    
    switch Module
        case 'OspreyLoad' %Load tab
            outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyLoad');
            [~,filename,~]  = fileparts(MRSCont.files{kk});
                % Grid for Plot and Data control sliders
                if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) % HBox for HERMES/HERCULES
                    Plot = uix.HBox('Parent', input_figure, 'BackgroundColor',colormapfig.Background, 'Units', 'normalized');
                else
                    Plot= uix.VBox('Parent', input_figure, 'BackgroundColor',colormapfig.Background, 'Units', 'normalized');
                end
                InfoText  = uicontrol('Parent',Info,'style','text',...
                                              'FontSize', 12, 'FontName', 'Arial','ForegroundColor', colormapfig.Foreground,...
                                              'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',colormapfig.Background);
            % Get parameter from file to fill the info panel
            if strcmp(which,'mets') %Is metabolite data?
                StatText = ['Metabolite Data -> Sequence: ' Seq '; B0: ' num2str(MRSCont.raw{1,kk}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,kk}.te) ' / ' num2str(MRSCont.raw{1,kk}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,kk}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw{1,kk}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,kk}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,kk}.averages)...
                             '; Sz: ' num2str(MRSCont.raw{1,kk}.sz) '; dimensions: ' num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{1})) ' x ' num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{2})) ' x ' num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{3})) ' mm = '...
                             num2str(MRSCont.raw{1,kk}.geometry.size.(Geom{1}) * MRSCont.raw{1,kk}.geometry.size.(Geom{2}) * MRSCont.raw{1,kk}.geometry.size.(Geom{3})/1000) ' ml'];
            else if strcmp(which, 'ref') %Is water or ref data?
            StatText = ['Reference Data -> Sequence: ' Seq '; B0: ' num2str(MRSCont.raw_ref{1,kk}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,kk}.te) ' / ' num2str(MRSCont.raw_ref{1,kk}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,kk}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw_ref{1,kk}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,kk}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,kk}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_ref{1,kk}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,kk}.geometry.size.(Geom{1})) ' x ' num2str(MRSCont.raw_ref{1,kk}.geometry.size.(Geom{2})) ' x ' num2str(MRSCont.raw_ref{1,kk}.geometry.size.(Geom{3})) ' mm = '...
                             num2str(MRSCont.raw_ref{1,kk}.geometry.size.(Geom{1}) * MRSCont.raw_ref{1,kk}.geometry.size.(Geom{2}) * MRSCont.raw_ref{1,kk}.geometry.size.(Geom{3})/1000) ' ml'];
                else
                    StatText = ['Water Data -> Sequence: ' Seq '; B0: ' num2str(MRSCont.raw_w{1,kk}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,kk}.te) ' / ' num2str(MRSCont.raw_w{1,kk}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,kk}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw_w{1,kk}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,kk}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,kk}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_w{1,kk}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,kk}.geometry.size.(Geom{1})) ' x ' num2str(MRSCont.raw_w{1,kk}.geometry.size.(Geom{2})) ' x ' num2str(MRSCont.raw_w{1,kk}.geometry.size.(Geom{3})) ' mm = '...
                             num2str(MRSCont.raw_w{1,kk}.geometry.size.(Geom{1}) * MRSCont.raw_w{1,kk}.geometry.size.(Geom{2}) * MRSCont.raw_w{1,kk}.geometry.size.(Geom{3})/1000) ' ml'];
                end
            end
            set(InfoText, 'String', sprintf(StatText));
     %%% 4. VISUALIZATION PART OF THIS TAB %%%
     %osp_plotLoad is used to visualize the raw data. Number of subplots
     %depends on the number of subspectra of the seuqence
            if strcmp(which,'mets') %Metabolite data/tab
                outputFile      = [filename '_OspreyLoad_mets.eps'];
                temp = osp_plotLoad(MRSCont, kk,'mets',1 );
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
            else if strcmp(which,'ref') %ref data/tab
                    temp = osp_plotLoad(MRSCont, kk,'ref',1 );
                    ViewAxes = gca();
                    set( ViewAxes, 'Parent', Plot );
                    outputFile      = [filename '_OspreyLoad_ref.eps'];
                else %water data/tab has only one window all the time
                    temp = osp_plotLoad(MRSCont, kk,'w',1 );
                    ViewAxes = gca();
                    set(ViewAxes, 'Parent', Plot );
                    outputFile      = [filename '_OspreyLoad_w.eps'];
                end
            end
            set(input_figure, 'Heights', [-0.1 -0.9]);
            % Get rid of the Load figure
            close( temp );
        case 'OspreyProcess' %Process Tab
            outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyProcess');
            [~,filename,~]  = fileparts(MRSCont.files{kk});
            Plot = uix.HBox('Parent', input_figure, ...
                'Padding', 5,'BackgroundColor', colormapfig.Background);
            set(input_figure, 'Heights', [-0.11 -0.89]);
            if MRSCont.flags.isUnEdited %Is UnEdited?
                if strcmp(which, 'ref') || strcmp(which, 'w')
                    SNR = 'water';
                else                    
                    SNR = 'tNAA';
                end
            end
            if MRSCont.flags.isMEGA %Is MEGA?
                if strcmp(which, 'ref') || strcmp(which, 'w')
                    SNR = 'water';
                else 
                    switch which
                        case 'A'
                            SNR = 'tNAA';
                        case 'B'
                            SNR = 'tCr';
                        case 'diff1'
                            SNR = MRSCont.processed.diff1{1,kk}.target;
                        case 'sum'
                            SNR = 'tNAA';                            
                    end                        
                end
            end
        if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) %Is HERMES\HERCULES?
           if strcmp(which, 'ref') || strcmp(which, 'w')
                    SNR = 'water';
                else 
                    switch which
                        case 'A'
                            SNR = 'tNAA';
                        case 'B'
                            SNR = 'tCr';
                        case 'C'
                            SNR = 'tNAA';
                        case 'D'
                            SNR = 'tCr';                            
                        case 'diff1'
                            SNR = MRSCont.processed.diff1{1,kk}.target;
                        case 'diff2'
                            SNR = MRSCont.processed.diff1{1,kk}.target;
                        case 'sum'
                            SNR = 'tNAA';
                    end                        
           end
        end                        
            
            % Get parameter from file to fill the info panel
            if (strcmp(which,'A') || strcmp(which,'B') || strcmp(which,'C') || strcmp(which,'D') || strcmp(which,'diff1') || strcmp(which,'diff2') || strcmp(which,'sum'))
                StatText = ['Metabolite Data -> SNR(' SNR '): '  num2str(MRSCont.QM.SNR.(which)(kk)) '; FWHM: '...
                            num2str(MRSCont.QM.FWHM.(which)(kk)) ' / ' (num2str(MRSCont.QM.FWHM.(which)(kk)*MRSCont.processed.(which){kk}.txfrq/1e6))...
                            ' ppm / Hz \nReference shift: ' num2str(MRSCont.QM.freqShift.(which)(kk)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.(which)(kk)*MRSCont.processed.(which){kk}.txfrq/1e6)...
                            ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.(which)(kk)*MRSCont.processed.(which){kk}.txfrq/1e6) ' Hz'];
            else if strcmp(which,'ref')
            StatText = ['Reference Data -> SNR(' SNR '): ' num2str(MRSCont.QM.SNR.(which)(kk)) '; FWHM: '...
                        num2str(MRSCont.QM.FWHM.(which)(kk)) ' / ' (num2str(MRSCont.QM.FWHM.(which)(kk)*MRSCont.processed.(which){kk}.txfrq/1e6))...
                        ' ppm / Hz'];
                else
                    StatText = ['Water Data -> SNR(' SNR '): ' num2str(MRSCont.QM.SNR.(which)(kk)) '; FWHM: '...
                                num2str(MRSCont.QM.FWHM.(which)(kk)) '/' (num2str(MRSCont.QM.FWHM.(which)(kk)*MRSCont.processed.(which){kk}.txfrq/1e6))...
                                ' ppm / Hz'];
                end
            end
            InfoText  = uicontrol('Parent',Info,'style','text','FontSize', 12, 'FontName', 'Arial',...
                                         'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                         'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);

 %%% 4. VISUALIZATION PART OF THIS TAB %%%
 %osp_plotProcess is used to visualize the processed spectra
            temp = osp_plotProcess(MRSCont, kk,which,1 ); % Create figure
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
            
            outputFile      = [filename '_OspreyProcess_' which '.eps'];
        case 'OspreyFit' %Fit
             outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyFit');
            [~,filename,~]  = fileparts(MRSCont.files{kk});

            Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',colormapfig.Background);
            set(input_figure, 'Heights', [-0.12 -0.88]);
            % Get parameter from file to fill the info panel
            if  ~strcmp (which, 'ref') && ~strcmp (which, 'w') %Metabolite data?
                StatText = ['Metabolite Data -> Sequence: ' Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' which,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.range(1)) ' to ' num2str(MRSCont.opts.fit.range(2)) ' ppm; Baseline knot spacing: ' num2str(MRSCont.opts.fit.bLineKnotSpace) ...
                        ' ppm\nNumber of metabolites: ' num2str(MRSCont.fit.resBasisSet.(which){1,MRSCont.info.A.unique_ndatapoint_indsort(kk)}.nMets) '; Number of macro moclecules: ' num2str(MRSCont.fit.resBasisSet.(which){1,MRSCont.info.A.unique_ndatapoint_indsort(kk)}.nMM) ...
                        ' scale: '  num2str(MRSCont.fit.scale{kk})];
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
                                        'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                        'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
            Results = uix.Panel('Parent', Plot,...
                                       'Title', ['Raw Amplitudes'],'FontName', 'Arial','HighlightColor', colormapfig.Foreground,...
                                       'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground, 'ShadowColor', colormapfig.Foreground);
            if  ~strcmp (MRSCont.opts.fit.style, 'Concatenated') ||  strcmp(which, 'ref') || strcmp(which, 'w') %Is not concateneted or is reference/water fit 
                Style = which;
                RawAmpl = MRSCont.fit.results.(which).fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk};
            else %Is concatenated and not water/reference
                Style = 'conc';
                RawAmpl = MRSCont.fit.results.(Style).fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk};
            end
            if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                if ~(strcmp(Style, 'ref') || strcmp(Style, 'w')) %Metabolite fit
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl) %Names and Amplitudes
                        NameText = [NameText, [MRSCont.fit.resBasisSet.(Style){1,MRSCont.info.A.unique_ndatapoint_indsort(kk)}.name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                else %Water/reference fit but this should never happen in this loop
                   NameText = ['Water: ' ];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                end
                set(Results, 'Title', ['Raw Amplitudes']);
                    FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                    FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                    FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
            else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                if ~(strcmp(Style, 'ref') || strcmp(Style, 'w')) %Metabolite fit
                    if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results.ref.fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk});
                    else
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results.water.fitParams{1,kk}.ampl .* MRSCont.fit.scale{kk});
                    end
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl) %Names and Amplitudes
                        NameText = [NameText, [MRSCont.fit.resBasisSet.(Style){1,MRSCont.info.A.unique_ndatapoint_indsort(kk)}.name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                    set(Results, 'Title', ['Raw Water Ratio']);
                    FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                    FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                    FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                else %Water/reference fit
                   NameText = ['Water: ' ];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                   set(Results, 'Title', ['Raw Amplitudes']);
                   FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                   FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                   'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                   'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                   FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                   'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                   'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                end
            end
%%%  5. VISUALIZATION PART OF THIS TAB %%%
%osp_plotFit is used to visualize the fits (off,diff1,diff2,sum,ref,water)
            temp = figure( 'Visible', 'off' );
            temp = osp_plotFit(MRSCont, kk,Style,1,which);
            ViewAxes = gca();
            set(ViewAxes, 'Parent', Plot );
            close( temp );

            set(Plot,'Widths', [-0.16 -0.84]);
            set(Plot.Children(2), 'Units', 'normalized');
            set(Plot.Children(2), 'OuterPosition', [0.17,0.02,0.75,0.98])
            outputFile      = [filename '_OspreyFit_' Style '_' which '.eps'];
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
                'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);

        %%% 2.VISUALIZATION PART OF THIS TAB %%%
        % In this case osp_plotCoreg or osp_plotSegment is used to visualize the
        % coregistration or the segmentation
            Results = uix.VBox('Parent', Plot,'BackgroundColor',colormapfig.Background);
            temp = figure( 'Visible', 'off' );
            if MRSCont.flags.didSeg %Did segment. In this case coreg has already been performed. Visualize both
                osp_plotCoreg(MRSCont, kk, 1);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', Results );
                colormap(Results.Children,'gray');
                close( temp );
                temp = figure( 'Visible', 'off' );
                osp_plotSegment(MRSCont, kk, 1);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', Results );
                colormap(Results.Children(1),'gray');
                close( temp );
            else % Only coreg has been run
                osp_plotCoreg(MRSCont, kk, 1);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', Results );
                colormap(Results.Children,'gray');
                close( temp );
            end
            outputFile      = [filename '_OspreyCoregSeg.eps'];
        case 'OspreySpecOverview' %SpecOverview
            set(Info,'Title', 'Descriptive Information');
            groupString = '';
            for g = 1 : MRSCont.overview.NoGroups
                groupString = [groupString MRSCont.overview.groupNames{g} ' with ' num2str(sum(MRSCont.overview.groups(:) == g)) ' subjects; '];
            end
            StatText = ['Sequence: ' Seq '; Number of subjects: ' num2str(MRSCont.nDatasets) '; Number of Groups: ' num2str(MRSCont.overview.NoGroups) '\n'...
                         'Distribution: ' groupString '\n'];

           InfoText  = uicontrol('Parent',Info,'style','text',...
                'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
            Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',colormapfig.Background);
            set(input_figure, 'Heights', [-0.1 -0.9]);    
            outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview','Individual');
            outputFile  = [which '.eps']; 
            for g = 1 :  MRSCont.overview.NoGroups %Loop over groups
                temp = osp_plotOverviewSpec(MRSCont, which,1, g, 0.1);
                if g == 1
                    temp = get(temp,'Parent');
                    fig_hold = get(temp,'Parent');
                else
                    ax=get(temp,'Parent');
                    copyobj(ax.Children, fig_hold.Children(1));
                    close_fig= get(ax,'Parent');
                    close(close_fig);
                end

            end
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
                    'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                set(input_figure, 'Heights', [-0.1 -0.9]);                
                    outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Mean');
                    outputFile  = [which '.eps'];                    
                    for g = 1 :  MRSCont.overview.NoGroups %Loop over groups
                        if gui.overview.Number.Groups > 1
                            temp = osp_plotMeanSpec(MRSCont, which,1,g,0.1,1);
                            if g == 1
                                fig_hold = temp;
                            else
                                ax=get(temp,'Children');
                                lines=get(ax,'Children');
                                copyobj(lines, fig_hold.Children(1));
                                close(temp);
                            end   
                        else
                            fig_hold = osp_plotMeanSpec(MRSCont, which,0,g);
                        end
                    end
                    set(fig_hold.Children, 'Parent', Plot );                    
        case 'OspreyRaincloudOverview' %Raincloud plot
            set(Info,'Title', 'Descriptive Information');
                groupString = '';
                for g = 1 : MRSCont.overview.NoGroups
                    groupString = [groupString MRSCont.overview.groupNames{g} ' with ' num2str(sum(MRSCont.overview.groups(:) == g)) ' subjects; '];
                end
                StatText = ['Sequence: ' Seq '; Number of subjects: ' num2str(MRSCont.nDatasets) '; Number of Groups: ' num2str(MRSCont.overview.NoGroups) '\n'...
                             'Distribution: ' groupString '\n'];

               InfoText  = uicontrol('Parent',Info,'style','text',...
                    'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                set(input_figure, 'Heights', [-0.1 -0.9]);   
            outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Raincloud');
            split_which = strsplit(which,'-');
            if strcmp(split_which{2},'AlphaCorrWaterScaled') || strcmp(split_which{2},'AlphaCorrWaterScaledGroupNormed')
                metab = 'GABA';
            end
            fig_hold = osp_plotRaincloud(MRSCont,split_which{1},split_which{2},metab,'Raincloud plot',1);
            delete( fig_hold.Children(1));
            set( fig_hold.Children, 'Parent', Plot );
            set(out.Children.Children(1).Children,'Children',flipud(out.Children.Children(1).Children.Children));
            close(fig_hold);
            outputFile  = [metab '_' split_which{1} '_' split_which{2} '.eps'];  
        case 'OspreyScatterOverview' %Correlation plot
            set(Info,'Title', 'Descriptive Information');
                groupString = '';
                for g = 1 : MRSCont.overview.NoGroups
                    groupString = [groupString MRSCont.overview.groupNames{g} ' with ' num2str(sum(MRSCont.overview.groups(:) == g)) ' subjects; '];
                end
                StatText = ['Sequence: ' Seq '; Number of subjects: ' num2str(MRSCont.nDatasets) '; Number of Groups: ' num2str(MRSCont.overview.NoGroups) '\n'...
                             'Distribution: ' groupString '\n'];

               InfoText  = uicontrol('Parent',Info,'style','text',...
                    'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                'BackgroundColor',colormapfig.Background,'ForegroundColor', colormapfig.Foreground);
                Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',colormapfig.Background);
                set(input_figure, 'Heights', [-0.1 -0.9]);   
            outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Correlation');
            split_which = strsplit(which,'-');
            if strcmp(split_which{2},'AlphaCorrWaterScaled') || strcmp(split_which{2},'AlphaCorrWaterScaledGroupNormed')
                metab = 'GABA';
            end
            split_corr = strsplit(corr,'-');
            if strcmp(split_corr{1},'corr')
                fig_hold = osp_plotScatter(MRSCont,split_which{1},split_which{2},metab,MRSCont.overview.corr.Meas{str2num(split_corr{2})},MRSCont.overview.corr.Names{str2num(split_corr{2})},0);
                outputFile  = [metab '_' split_which{1} '_' split_which{2} '_'  MRSCont.overview.corr.Names{str2num(split_corr{2})} '.eps'];
            else if strcmp(split_corr{1},'metab')
                fig_hold = osp_plotScatter(MRSCont,split_which{1},split_which{2},metab,MRSCont.overview.corr.Names{str2num(split_corr{1})},metab,0);
                outputFile  = [metab '_' split_which{1} '_' split_which{2} '_'  metab '.eps'];
                else if strcmp(split_corr{1},'SNR')
                        fig_hold = osp_plotScatter(MRSCont,split_which{1},split_which{2},metab,MRSCont.QM.SNR.A','SNR Subspectrum A',0);
                        outputFile  = [metab '_' split_which{1} '_' split_which{2} '_SNR.eps'];
                    else if  strcmp(split_corr{1},'FWHM')
                        fig_hold = osp_plotScatter(MRSCont,split_which{1},split_which{2},metab,MRSCont.QM.FWHM.A','FWHM (Hz)',0);
                        outputFile  = [metab '_' split_which{1} '_' split_which{2} '_FWHM.eps'];
                        end
                    end
                end
            end
            delete( fig_hold.Children(1));
            delete( fig_hold.Children(1));
            set(fig_hold.Children, 'Parent', Plot );
            set(out.Children.Children(1).Children,'Children',flipud(out.Children.Children(1).Children.Children));
            close(fig_hold);
            
    end    
    set(out,'Renderer','painters','Menu','none','Toolbar','none');
%% Clean up and save

% Save the figure to the output folder
% Determine output folder
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
saveas(out,fullfile(outputFolder,outputFile),'epsc');
end