function osp_onPrint( ~, ~ ,gui)
%% osp_onPrint
%   Callback function on print figure button click.
%
%
%   USAGE:
%       osp_onPrint( ~, ~ ,gui);
%
%   INPUT:      gui      = gui class containing all handles and the MRSCont 
%
%   OUTPUT:     Changes in gui parameters and MRSCont are written into the
%               gui class
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-01-16)
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
    MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
    selectedTab = get(gui.layout.tabs, 'Selection');
    screenSize      = [1, 1, 1000, 900];
    canvasSize      = screenSize;
    canvasSize(4)   = screenSize(4) * 0.7;
    canvasSize(3)   = canvasSize(4) * (11/8.5);
    canvasSize(2)   = (screenSize(4) - canvasSize(4))/2;
    canvasSize(1)   = (screenSize(3) - canvasSize(3))/2;
    out = figure('NumberTitle', 'off', 'Visible', 'on', 'Menu', 'none','Position', canvasSize,...
                    'ToolBar', 'none', 'HandleVisibility', 'off', 'Renderer', 'painters', 'Color', gui.colormap.Background);
    switch selectedTab
        case 1
            Title = [MRSCont.ver.Osp ' ' MRSCont.ver.Load];
        case 2
            Title = [MRSCont.ver.Osp ' ' MRSCont.ver.Pro];
        case 3
            Title = [MRSCont.ver.Osp ' ' MRSCont.ver.Fit];
        case 4
            Title = [MRSCont.ver.Osp ' ' MRSCont.ver.Coreg];
        case 5
            Title = [MRSCont.ver.Osp ' ' MRSCont.ver.Over];
        otherwise
            Title = '';
    end
            
    Frame = uix.Panel('Parent',out, 'Padding', 1, 'Title', Title,...
                                 'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                 'HighlightColor', gui.colormap.Background, 'ShadowColor', gui.colormap.Background);
    input_figure = uix.VBox('Parent', Frame,  'BackgroundColor',gui.colormap.Background, 'Spacing', 5);                            
    box = uix.HBox('Parent', input_figure,'BackgroundColor',gui.colormap.Background, 'Spacing',6);
    Info = uix.Panel('Parent',box, 'Padding', 5, 'Title', MRSCont.files{gui.controls.Selected},...
                                 'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
    LogoFig = figure('Visible','off'); 
    [I, map] = imread('osprey.gif','gif');
    axes(LogoFig, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
    ViewAxes = gca();
    set( ViewAxes,'Box','off','XColor', 'none','YColor','none', 'Parent', box );

    set(box, 'Width', [-0.9 -0.1]);
    switch selectedTab
        case 1 %Load tab
            outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyLoad');
            [~,filename,~]  = fileparts(MRSCont.files{gui.controls.Selected});

                % Grid for Plot and Data control sliders
                if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) % HBox for HERMES/HERCULES
                    Plot = uix.HBox('Parent', input_figure, 'BackgroundColor',gui.colormap.Background, 'Units', 'normalized');
                else
                    Plot= uix.VBox('Parent', input_figure, 'BackgroundColor',gui.colormap.Background, 'Units', 'normalized');
                end
                InfoText  = uicontrol('Parent',Info,'style','text',...
                                              'FontSize', 12, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                              'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
            % Get parameter from file to fill the info panel
            if gui.load.Selected == 1 %Is metabolite data?
                StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                             '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                             num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
            else if gui.load.Selected == 2 %Is water or ref data?
            StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                             num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
                else
                    StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                             num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
                end
            end
            set(InfoText, 'String', sprintf(StatText));
     %%% 4. VISUALIZATION PART OF THIS TAB %%%
     %osp_plotLoad is used to visualize the raw data. Number of subplots
     %depends on the number of subspectra of the seuqence
            if gui.load.Selected == 1 %Metabolite data/tab
                outputFile      = [filename '_OspreyLoad_mets.pdf'];
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mets');
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
                    gui.layout.multiACload = uix.VBox('Parent', Plot, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                        gui.layout.multiAload = uix.VBox('Parent', gui.layout.multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                        gui.layout.multiCload = uix.VBox('Parent', gui.layout.multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    gui.layout.multiBDload = uix.VBox('Parent', Plot,'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                        gui.layout.multiBload = uix.VBox('Parent', gui.layout.multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                        gui.layout.multiDload = uix.VBox('Parent', gui.layout.multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    set( temp.Children(1), 'Parent', gui.layout.multiDload );
                    set( temp.Children(1), 'Parent', gui.layout.multiCload );
                    set( temp.Children(1), 'Parent', gui.layout.multiBload );
                    set( temp.Children(1), 'Parent', gui.layout.multiAload );
                    set(Plot,'Width', [-0.49 -0.49]);
                    set(gui.layout.multiDload.Children(1), 'Units', 'normalized')
                    set(gui.layout.multiDload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(gui.layout.multiCload.Children(1), 'Units', 'normalized')
                    set(gui.layout.multiCload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(gui.layout.multiBload.Children(1), 'Units', 'normalized')
                    set(gui.layout.multiBload.Children(1), 'OuterPosition', [0,0,1,1])
                    set(gui.layout.multiAload.Children(1), 'Units', 'normalized')
                    set(gui.layout.multiAload.Children(1), 'OuterPosition', [0,0,1,1])

                end
            else if gui.load.Selected == 2 %ref data/tab
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref');
                    ViewAxes = gca();
                    set( ViewAxes, 'Parent', Plot );
                    outputFile      = [filename '_OspreyLoad_ref.pdf'];
                else %water data/tab has only one window all the time
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w');
                    ViewAxes = gca();
                    set(ViewAxes, 'Parent', Plot );
                    outputFile      = [filename '_OspreyLoad_w.pdf'];
                end
            end
            set(input_figure, 'Heights', [-0.1 -0.9]);
            % Get rid of the Load figure
            close( temp );

        case 2 %Process Tab
            outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyProcess');
            [~,filename,~]  = fileparts(MRSCont.files{gui.controls.Selected});
            Selection = gui.process.Names{gui.process.Selected};            

            Plot = uix.HBox('Parent', input_figure, ...
                'Padding', 5,'BackgroundColor', gui.colormap.Background);
            set(input_figure, 'Heights', [-0.11 -0.89]);
            % Get parameter from file to fill the info panel
            if (strcmp(gui.process.Names{gui.process.Selected},'A') || strcmp(gui.process.Names{gui.process.Selected},'B') || strcmp(gui.process.Names{gui.process.Selected},'C') || strcmp(gui.process.Names{gui.process.Selected},'D') || strcmp(gui.process.Names{gui.process.Selected},'diff1') || strcmp(gui.process.Names{gui.process.Selected},'diff2') || strcmp(gui.process.Names{gui.process.Selected},'sum') )
                StatText = ['Metabolite Data -> SNR(' gui.process.SNR{gui.process.Selected} '): '  num2str(MRSCont.QM.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                            num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM.freqShift.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq/1e6)...
                            ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq/1e6) ' Hz'];
            else if (strcmp(gui.process.Names{gui.process.Selected},'ref') || strcmp(gui.process.Names{gui.process.Selected},'mm'))
            StatText = ['Reference Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                        num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                        ' Hz / ppm'];
                else
                    StatText = ['Water Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                                num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '/' (num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm'];
                end
            end
            InfoText  = uicontrol('Parent',Info,'style','text','FontSize', 12, 'FontName', 'Arial',...
                                         'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                         'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

 %%% 4. VISUALIZATION PART OF THIS TAB %%%
 %osp_plotProcess is used to visualize the processed spectra
            temp = osp_plotProcess(MRSCont, gui.controls.Selected,gui.process.Names{gui.process.Selected}); % Create figure
            %Subplots are distributed here
                proSpecs = uix.VBox('Parent', Plot, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    proPre = uix.VBox('Parent', proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    proPost = uix.VBox('Parent', proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                proOut = uix.VBox('Parent', Plot,'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    proDrift = uix.VBox('Parent', proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    proAlgn = uix.VBox('Parent', proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);

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
            
            outputFile      = [filename '_OspreyProcess_' Selection '.pdf'];
        case 3 %Fit
             outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyFit');
            [~,filename,~]  = fileparts(MRSCont.files{gui.controls.Selected});
            Selection = gui.fit.Names{gui.fit.Selected};
            Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
            set(input_figure, 'Heights', [-0.12 -0.88]);
            if  ~strcmp (MRSCont.opts.fit.style, 'Concatenated') ||  strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') || strcmp(gui.fit.Names{gui.fit.Selected}, 'w') %Is not concateneted or is reference/water fit 
            gui.fit.Style = Selection;
            else %Is concatenated and not water/reference
                gui.fit.Style = 'conc';
            end
            if ~strcmp(gui.fit.Names{gui.fit.Selected}, 'ref') && ~strcmp(gui.fit.Names{gui.fit.Selected}, 'w')
                refShift = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refShift;
                refFWHM = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.refFWHM;
                ph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph0;
                ph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ph1;
                iniph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph0;
                iniph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph1;  
            end
            RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};

            % Get parameter from file to fill the info panel
            if  ~strcmp (Selection, 'ref') && ~strcmp (Selection, 'w') %Metabolite data
            iniph0 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph0;
            iniph1 = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.prelimParams.ph1;
            StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' Selection,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.range(1)) ' to ' num2str(MRSCont.opts.fit.range(2)) ' ppm; Baseline knot spacing: ' num2str(MRSCont.opts.fit.bLineKnotSpace) ' ppm; ph0: ' num2str(ph0,'%1.2f') '°; ph1: ' num2str(ph1,'%1.2f') '°; refShift: ' num2str(refShift,'%1.2f') ' Hz; refFWHM: ' num2str(refFWHM,'%1.2f')...
                        ' ppm\nNumber of metabolites: ' num2str(MRSCont.fit.resBasisSet.(gui.fit.Style){1,MRSCont.info.A.unique_ndatapoint_indsort(gui.controls.Selected)}.nMets) '; Number of MM/lipids: ' num2str(MRSCont.fit.resBasisSet.(gui.fit.Style){1,MRSCont.info.A.unique_ndatapoint_indsort(gui.controls.Selected)}.nMM) ...
                        ' scale: '  num2str(MRSCont.fit.scale{gui.controls.Selected}) '; initial ph0: ' num2str(iniph0,'%1.2f') '°; initial ph1: ' num2str(iniph1,'%1.2f') '°'];
            else if strcmp (Selection, 'ref') %Reference data?
            StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' Selection,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.rangeWater(1)) ' to ' num2str(MRSCont.opts.fit.rangeWater(2)) ' ppm'];
                else %Is water data
                    StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; Fitting algorithm: ' MRSCont.opts.fit.method  '; Fitting Style: ' MRSCont.opts.fit.style '; Selected subspecs: ' Selection,...
                        '\nFitting range: ' num2str(MRSCont.opts.fit.rangeWater(1)) ' to ' num2str(MRSCont.opts.fit.rangeWater(2)) ' ppm'];
                end
            end
 %%% 4. FILLING FITTED AMPLITUDE PANEL %%%
 % Creates the panel on the right side with the fitted ammplitudes
            InfoText  = uicontrol('Parent',Info,'style','text',...
                                        'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            Results = uix.Panel('Parent', Plot,...
                                       'Title', ['Raw Amplitudes'],'FontName', 'Arial','HighlightColor', gui.colormap.Foreground,...
                                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl) %Names and Amplitudes
                        NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style){1,MRSCont.info.A.unique_ndatapoint_indsort(gui.controls.Selected)}.name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                else %Water/reference fit but this should never happen in this loop
                   NameText = ['Water: ' ];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                end
                set(Results, 'Title', ['Raw Amplitudes']);
                    FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                    if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                    else
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results.w.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                    end
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl) %Names and Amplitudes
                        NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style){1,MRSCont.info.A.unique_ndatapoint_indsort(gui.controls.Selected)}.name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                    set(Results, 'Title', ['Raw Water Ratio']);
                    FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                else %Water/reference fit
                   NameText = ['Water: ' ];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                   set(Results, 'Title', ['Raw Amplitudes']);
                   FitText = uix.HBox('Parent', Results, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                   FitTextNames  = uicontrol('Parent',FitText,'style','text',...
                   'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                   'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                   FitTextAmpl  = uicontrol('Parent',FitText,'style','text',...
                   'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                   'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                end
            end
%%%  5. VISUALIZATION PART OF THIS TAB %%%
%osp_plotFit is used to visualize the fits (off,diff1,diff2,sum,ref,water)
            temp = figure( 'Visible', 'off' );
            temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,Selection);
            ViewAxes = gca();
            set(ViewAxes, 'Parent', Plot );
            close( temp );

            set(Plot,'Widths', [-0.16 -0.84]);
            set(Plot.Children(2), 'Units', 'normalized');
            set(Plot.Children(2), 'OuterPosition', [0.17,0.02,0.75,0.98])
            outputFile      = [filename '_OspreyFit_' gui.fit.Style '_' Selection '.pdf'];
        case 4 %Coreg/Seg
            outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyCoregSeg');
            addpath(genpath([gui.folder.spmversion filesep])); % Add SPM  path
            [~,filename,~]  = fileparts(MRSCont.files{gui.controls.Selected});
            
            % Creates layout for plotting and data control
            Plot = uix.HBox('Parent', input_figure,'BackgroundColor',gui.colormap.Background);
            set(input_figure, 'Heights', [-0.1 -0.9]);
            % Get parameter from file to fill the info panel

            StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '\nraw subspecs: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];

           InfoText  = uicontrol('Parent',Info,'style','text',...
                'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

        %%% 2.VISUALIZATION PART OF THIS TAB %%%
        % In this case osp_plotCoreg or osp_plotSegment is used to visualize the
        % coregistration or the segmentation
            Results = uix.VBox('Parent', Plot,'BackgroundColor',gui.colormap.Background);
            temp = figure( 'Visible', 'off' );
            if MRSCont.flags.didSeg %Did segment. In this case coreg has already been performed. Visualize both
                osp_plotCoreg(MRSCont, gui.controls.Selected);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', Results );
                colormap(Results.Children,'gray')
                close( temp );
                temp = figure( 'Visible', 'off' );
                osp_plotSegment(MRSCont, gui.controls.Selected);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', Results );
                colormap(Results.Children(1),'gray');
                close( temp );
            else % Only coreg has been run
                osp_plotCoreg(MRSCont, gui.controls.Selected);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', Results );
                colormap(Results.Children,'gray');
                close( temp );
            end            
            outputFile      = [filename '_OspreyCoregSeg.pdf'];
            rmpath(genpath([gui.folder.spmversion filesep])); %Remove SPM path
        case 6 %Overview
            ovSelection = get(gui.layout.overviewTab, 'Selection');
            set(Info,'Title', 'Descriptive Information');
            groupString = '';
            for g = 1 : MRSCont.overview.NoGroups
                groupString = [groupString MRSCont.overview.groupNames{g} ' with ' num2str(sum(MRSCont.overview.groups(:) == g)) ' subjects; '];
            end
            StatText = ['Sequence: ' gui.load.Names.Seq '; Number of subjects: ' num2str(MRSCont.nDatasets) '; Number of Groups: ' num2str(MRSCont.overview.NoGroups) '\n'...
                         'Distribution: ' groupString '\n'];

           InfoText  = uicontrol('Parent',Info,'style','text',...
                'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            Plot = uix.HBox('Parent', input_figure, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
            set(input_figure, 'Heights', [-0.1 -0.9]);
            switch ovSelection
                case 1 %SpecOverview
                    Selection = gui.controls.pop_specsOvPlot.String(gui.process.Selected);
                    outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview','Individual');
                    if gui.controls.GM == 0
                        outputFile  = [Selection{1} '.pdf']; 
                        for g = 1 :  gui.overview.Number.Groups %Loop over groups
                            temp = osp_plotOverviewSpec(MRSCont, Selection{1},g, gui.layout.shiftind);
                            if g == 1
                                temp = get(temp,'Parent');
                                fig_hold = get(temp,'Parent');
                            else
                                ax=get(temp,'Parent');
                                copyobj(ax.Children, fig_hold.Children(1));
                                set(fig_hold.Children, 'Parent', Plot );
                                close_fig= get(ax,'Parent');
                            end

                        end
                    else
                        outputFile  = [Selection{1} 'Grand_mean.pdf']; 
                       temp = osp_plotOverviewSpec(MRSCont, Selection{1},'GMean', gui.layout.shiftind);
                       temp = get(temp,'Parent');
                       fig_hold = get(temp,'Parent'); 
                       set(fig_hold.Children, 'Parent', Plot );
                    end                   
                    close(fig_hold);
                case 2 %MeanOverview
                    Selection = gui.controls.pop_meanOvPlot.String(gui.process.Selected);
                    outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Mean');
                    if gui.controls.GM == 0
                        outputFile  = [Selection{1} '.pdf'];
                        for g = 1 :  gui.overview.Number.Groups
                            if gui.overview.Number.Groups > 1
                                temp = osp_plotMeanSpec(MRSCont, Selection{1},g,1,1/gui.overview.Number.Groups);
                                if g == 1
                                    fig_hold = temp;
                                else
                                    ax=get(temp,'Children');
                                    lines=get(ax,'Children');
                                    copyobj(lines, fig_hold.Children(1));
                                    close(temp);
                                end   
                            else
                                fig_hold = osp_plotMeanSpec(MRSCont, Selection{1},g);
                            end
                        end
                    else
                        outputFile  = [Selection{1} 'Grand_mean.pdf'];
                        fig_hold = osp_plotMeanSpec(MRSCont, Selection{1},'GMean', 1);
                    end
                    set(fig_hold.Children, 'Parent', Plot );
                    close(fig_hold);
                    
                case 4 %Raincloud plot
                    outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Raincloud');
                    Selection = gui.quant.popMenuNames{gui.quant.Selected.Quant};
                    if ~strcmp(Selection,'Quality')  
                        split_Selection = strsplit(Selection,'-');
                        if strcmp(split_Selection{2},'AlphaCorrWaterScaled') || strcmp(split_Selection{2},'AlphaCorrWaterScaledGroupNormed')
                            metab = 'GABA';
                        else
                            metab = MRSCont.quantify.metabs{gui.overview.Selected.Metab};
                        end
                        if ~gui.controls.GM
                            fig_hold = osp_plotRaincloud(MRSCont,split_Selection{1},split_Selection{2},metab,'Raincloud plot');
                        else
                            fig_hold = osp_plotRaincloud(MRSCont,split_Selection{1},split_Selection{2},metab,'Raincloud plot',1);
                        end
                    else
                       quality = {'SNR','FWHM','freqShift'};
                       if ~gui.controls.GM
                            fig_hold = osp_plotRaincloud(MRSCont,'Quality','Quality',quality{gui.overview.Selected.Metab},'Raincloud plot');  
                       else
                            fig_hold = osp_plotRaincloud(MRSCont,'Quality','Quality',quality{gui.overview.Selected.Metab},'Raincloud plot',1); 
                       end
                        split_Selection{2}=quality{gui.overview.Selected.Metab};
                        split_Selection{1}='Quality';
                        metab = 'Spectral';                        
                    end                    
                    delete( fig_hold.Children(1));
                    set( fig_hold.Children, 'Parent', Plot );
                    set(out.Children.Children(1).Children(1).Children,'Children',flipud(out.Children.Children(1).Children(1).Children.Children));
                    close(fig_hold);
                    if ~gui.controls.GM
                        outputFile  = [metab '_' split_Selection{1} '_' split_Selection{2} '.pdf'];  
                    else
                        outputFile  = [metab '_' split_Selection{1} '_' split_Selection{2} 'Grand_mean.pdf'];  
                    end
                case 5 %Correlation plot
                    outputFolder    = fullfile(MRSCont.outputFolder,'Figures','OspreyOverview', 'Correlation');
                    Selection = gui.quant.popMenuNames{gui.quant.Selected.Quant};
                    split_Selection = strsplit(Selection,'-');
                    MRSCont.flags.isGUI =0;
                    if strcmp(split_Selection{2},'AlphaCorrWaterScaled') || strcmp(split_Selection{2},'AlphaCorrWaterScaledGroupNormed')
                        metab = 'GABA';
                    else
                        metab = MRSCont.quantify.metabs{gui.overview.Selected.Metab};
                    end 
                    if gui.overview.Selected.CorrChoice == 1
                        fig_hold = osp_plotScatter(MRSCont,split_Selection{1},split_Selection{2},metab,gui.overview.CorrMeas{gui.overview.Selected.Corr},gui.overview.Names.Corr{gui.overview.Selected.Corr});
                        outputFile  = [metab '_' split_Selection{1} '_' split_Selection{2} '_'  gui.overview.Names.Corr{gui.overview.Selected.Corr} '.pdf'];
                    else if gui.overview.Selected.CorrChoice == 2
                        fig_hold = osp_plotScatter(MRSCont,split_Selection{1},gui.quant.Names.Quants{gui.quant.Selected.Quant},MRSCont.quantify.metabs{gui.overview.Selected.Metab},metab,metab);
                        outputFile  = [metab '_' split_Selection{1} '_' split_Selection{2} '_'  MRSCont.quantify.metabs{gui.overview.Selected.Metab} '.pdf'];
                        else
                            switch gui.overview.Selected.Corr
                                case 1
                                fig_hold = osp_plotScatter(MRSCont,split_Selection{1},split_Selection{2},metab,MRSCont.QM.SNR.A',gui.overview.Names.QM{gui.overview.Selected.Corr});
                                outputFile  = [metab '_' split_Selection{1} '_' split_Selection{2} '_'  gui.overview.Names.QM{gui.overview.Selected.Corr} '.pdf'];
                                case 2
                                fig_hold = osp_plotScatter(MRSCont,split_Selection{1},split_Selection{2},metab,MRSCont.QM.FWHM.A',gui.overview.Names.QM{gui.overview.Selected.Corr});
                                outputFile  = [metab '_' split_Selection{1} '_' split_Selection{2} '_'  gui.overview.Names.QM{gui.overview.Selected.Corr} '.pdf'];
                            end
                        end
                    end
                    delete( fig_hold.Children(1));
                    delete( fig_hold.Children(1));
                    set(fig_hold.Children, 'Parent', Plot );
                    set(out.Children.Children.Children(1).Children,'Children',flipud(out.Children.Children.Children(1).Children.Children));
                    close(fig_hold);
                    MRSCont.flags.isGUI =1;
                    end
            
    end    
    set(out,'Renderer','painters','Menu','none','Toolbar','none');
    setappdata(gui.figure,'MRSCont',MRSCont);   % Write MRSCont into hidden container in gui class 
%% Clean up and save

% Save the figure to the output folder
% Determine output folder
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
fig_pos = out.PaperPosition;
out.PaperSize = [fig_pos(3) fig_pos(4)];

saveas(out,fullfile(outputFolder,outputFile),'pdf');
close(out);
end