function osp_iniFitWindow(gui)
%% osp_iniFitWindow
%   This function creates the inital fitting window in the gui.
%
%
%   USAGE:
%       osp_iniFitWindow(gui);
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
%%% 1. GET HANDLES %%%
% This function initializes the fit tab 
        MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
        gui.layout.tabs.TabEnables{3} = 'on';
        gui.layout.tabs.Selection  = 3;
        gui.layout.EmptyFitPlot = 0;
%%% 2. CREATING SUB TABS FOR THIS TAB %%%%
% In this case one tab fo each fit (off,sum,diff1,diff2,ref,water)
         gui.layout.fitTab.TabWidth   = 115;
         for t = 1 : gui.fit.Number %Create tabs depending on the number of fits
                gui.layout.(['fitTab' gui.fit.Names{t}]) = uix.VBox('Parent', gui.layout.fitTab,...
                                                            'BackgroundColor',gui.colormap.Background);
                gui.layout.fitTabhandles{t} = ['fitTab' gui.fit.Names{t}];
         end
        gui.layout.fitTab.TabTitles  = gui.fit.Names;
%%% 3. FILLING INFO PANEL FOR THIS TAB %%%%
% All the information from the Raw data is read out here
        for t = 1 : gui.fit.Number %Loop over fits
            % Parameter shown in the info panel on top
            gui.Info.fit = uix.Panel('Parent',  gui.layout.(gui.layout.fitTabhandles{t}), ...
                                    'Padding', 5, 'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}],...
                                    'FontName', 'Arial','HighlightColor', gui.colormap.Foreground,'BackgroundColor',...
                                    gui.colormap.Background,'ForegroundColor',gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            % Creates layout for plotting and data control
            gui.Plot.fit = uix.HBox('Parent', gui.layout.(gui.layout.fitTabhandles{t}), ...
                                   'Padding', 5,'BackgroundColor',gui.colormap.Background);
            set(gui.layout.(gui.layout.fitTabhandles{t}), 'Heights', [-0.1 -0.9]);
            % Get parameter from file to fill the info panel
            if gui.load.Selected == 1 %Is metabolite data?
                StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                             '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                             num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
            else if gui.load.Selected == 2 %Is reference data?
            StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                             num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
                else %Is water data
                    StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                             '\nraw subspecs: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                             num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
                end
            end
 %%% 4. FILLING FITTED AMPLITUDE PANEL %%%
 % Creates the panel on the right side with the fitted ammplitudes
            gui.InfoText.fit  = uicontrol('Parent',gui.Info.fit,'style','text',...
                                        'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            gui.Results.fit = uix.Panel('Parent', gui.Plot.fit,...
                                       'Title', ['Raw Amplitudes'],'FontName', 'Arial','HighlightColor', gui.colormap.Foreground,...
                                       'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            if  ~strcmp (MRSCont.opts.fit.style, 'Concatenated') ||  strcmp(gui.fit.Names{t}, 'ref') || strcmp(gui.fit.Names{t}, 'w') %Is not concateneted or is reference/water fit 
                gui.fit.Style = gui.fit.Names{t};
                RawAmpl = MRSCont.fit.results.(gui.fit.Names{t}).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
            else %Is concatenated and not water/reference
                gui.fit.Style = 'conc';
                RawAmpl = MRSCont.fit.results.(gui.fit.Style).fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected};
            end
            if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater) %Raw amplitudes are reported as no water/reference fitting was performed
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl) %Names and Amplitudes
                        NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fitStyle){1,gui.controls.Selected}.name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                else %Water/reference fit but this should never happen in this loop
                   NameText = ['Water: ' ];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                end
                set(gui.fitResults, 'Title', ['Raw Amplitudes']);
                    gui.Results.FitText = uix.HBox('Parent', gui.Results.fit, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            else %If water/reference data is fitted Raw amplitudes are calculated with regard to water
                if ~(strcmp(gui.fit.Style, 'ref') || strcmp(gui.fit.Style, 'w')) %Metabolite fit
                    if MRSCont.flags.hasRef %Calculate Raw Water Scaled amplitudes
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results.ref.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                    else
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results.water.fitParams{1,gui.controls.Selected}.ampl .* MRSCont.fit.scale{gui.controls.Selected});
                    end
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl) %Names and Amplitudes
                        NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.fit.Style){1,gui.controls.Selected}.name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                    set(gui.Results.fit, 'Title', ['Raw Water Ratio']);
                    gui.Results.FitText = uix.HBox('Parent', gui.Results.fit, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                else %Water/reference fit
                   NameText = ['Water: ' ];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                   set(gui.Results.fit, 'Title', ['Raw Amplitudes']);
                   gui.Results.FitText = uix.HBox('Parent', gui.Results.fit, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                   gui.Results.FitTextNames  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                   'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                   'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                   gui.Results.FitTextAmpl  = uicontrol('Parent',gui.Results.FitText,'style','text',...
                   'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                   'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                end
            end
%%%  5. VISUALIZATION PART OF THIS TAB %%%
%osp_plotFit is used to visualize the fits (off,diff1,diff2,sum,ref,water)
            temp = figure( 'Visible', 'off' );
            temp = osp_plotFit(MRSCont, gui.controls.Selected,gui.fit.Style,1,gui.fit.Names{t});
            ViewAxes = gca();
            set(ViewAxes, 'Parent', gui.Plot.fit );
            close( temp );
%%% 6. DATA CONTROLS FOR THIS TAB %%%
%Creates the scroll bar for the datasets of the MRSCont on the right side
            set(gui.Plot.fit,'Widths', [-0.16 -0.84]);
        end
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end