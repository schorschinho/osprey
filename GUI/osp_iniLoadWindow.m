function osp_iniLoadWindow(gui) 
%% osp_iniLoadWindow
%   This function creates the inital load window in the gui.
%
%
%   USAGE:
%       osp_iniLoadWindow(gui);
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
% This function creates the initial load window
        MRSCont = getappdata(gui.figure,'MRSCont');  % Get MRSCont from hidden container in gui class
        if MRSCont.flags.hasRef %Get variables regarding Subspectra
            gui.controls.Number = gui.controls.Number + 1;
            gui.load.Names.Spec{2} = 'reference';
        if MRSCont.flags.hasWater
            gui.controls.Number = gui.controls.Number + 1;
            gui.load.Names.Spec{3} = 'water';
        end
        else if MRSCont.flags.hasWater
            gui.controls.Number = gui.controls.Number + 1;
            gui.load.Names.Spec{2} = 'water';
            end
        end
        gui.layout.tabs.TabEnables{1} = 'on';
        gui.layout.tabs.Selection  = 1;
        gui.layout.EmptydataPlot = 0;
 %%% 2. CREATING SUB TABS FOR THIS TAB %%%
 % In this case one tab for each subspec (A,B,C,D,ref,water)
            gui.layout.metabLoTab = uix.VBox('Parent', gui.layout.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.layout.rawTab.TabWidth   = 115;
            gui.layout.rawTab.Selection  = 1;
            gui.layout.rawTabhandles = {'metabLoTab'};
            if gui.controls.Number == 2
                if MRSCont.flags.hasRef
                    gui.layout.refLoTab = uix.VBox('Parent', gui.layout.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.layout.rawTab.TabTitles  = gui.load.Names.Spec;
                    gui.layout.rawTab.TabEnables = {'on', 'on'};
                    gui.layout.rawTabhandles = {'metabLoTab', 'refLoTab'};
                end
                if MRSCont.flags.hasWater
                    gui.layout.wLoTab = uix.VBox('Parent', gui.layout.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.layout.rawTab.TabTitles  = gui.load.Names.Spec;
                    gui.layout.rawTab.TabEnables = {'on', 'on'};
                    gui.layout.rawTabhandles = {'metabLoTab', 'wLoTab'};
                end
            end
            if gui.controls.Number == 3
                gui.layout.refLoTab = uix.VBox('Parent', gui.layout.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                gui.layout.wLoTab = uix.VBox('Parent', gui.layout.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                gui.layout.rawTab.TabTitles  = gui.load.Names.Spec;
                gui.layout.rawTab.TabEnables = {'on', 'on','on'};
                gui.layout.rawTabhandles = {'metabLoTab', 'refLoTab', 'wLoTab'};
            end
 %%% 3. FILLING INFO PANEL FOR THIS TAB %%%
 % All the information from the Raw data is read out here
        for t = gui.controls.Number : -1 : 1 % Loop over subspectra & tabs
            % Parameter shown in the info panel on top
            gui.Info.data = uix.Panel('Parent', gui.layout.(gui.layout.rawTabhandles{t}), ...
                                     'Padding', 5, 'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}],...
                                     'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                     'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            % Grid for Plot and Data control sliders
            gui.Plot.data = uix.VBox('Parent', gui.layout.(gui.layout.rawTabhandles{t}), 'Padding', 5, 'BackgroundColor',gui.colormap.Background, 'Units', 'normalized');
            gui.InfoText.data  = uicontrol('Parent',gui.Info.data,'style','text',...
                                          'FontSize', 12, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
            set(gui.layout.(gui.layout.rawTabhandles{t}), 'Heights', [-0.1 -0.9]);

        % Get parameter from file to fill the info panel
        if gui.load.Selected == 1 %Is metabolite data?
            StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
        else if gui.load.Selected == 2 %Is water or ref data?
        StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
            else
                StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
            end
        end
        set(gui.InfoText.data, 'String', sprintf(StatText));
 %%% 4. VISUALIZATION PART OF THIS TAB %%%
 %osp_plotLoad is used to visualize the raw data. Number of subplots
 %depends on the number of subspectra of the seuqence
        temp = figure( 'Visible', 'off' );
        if t == 1 %Metabolite data/tab
            temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mets',1 );
            if MRSCont.flags.isUnEdited % One window for UnEdited
                ViewAxes = gca();
                set( ViewAxes, 'Parent', gui.Plot.data );
            end
            if MRSCont.flags.isMEGA %Two windows for MEGA
                set( temp.Children(2), 'Parent', gui.Plot.data );
                set( temp.Children(1), 'Parent', gui.Plot.data );
                set(gui.Plot.data,'Heights', [-0.49 -0.49]);
                set(gui.Plot.data.Children(2), 'Units', 'normalized')
                set(gui.Plot.data.Children(2), 'OuterPosition', [0,0.5,1,0.5])
                set(gui.Plot.data.Children(1), 'Units', 'normalized')
                set(gui.Plot.data.Children(1), 'OuterPosition', [0,0,1,0.5])
            end
        else if t == 2 %ref data/tab
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref',1 );
                ViewAxes = gca();
                set( ViewAxes, 'Parent', gui.Plot.data );
            else %water data/tab has only one window all the time
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',1 );
                ViewAxes = gca();
                set(ViewAxes, 'Parent', gui.Plot.data );
            end
        end

        % Get rid of the Load figure
        close( temp );
        end
        setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end