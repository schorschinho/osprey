function osp_iniLoadWindow(gui) 
%% osp_iniLoadWindow
%   This function creates the initial load window in the gui.
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
        if MRSCont.flags.hasMM %Get variables regarding Subspectra %re_mm
            gui.controls.Number = gui.controls.Number + 1;%re_mm
            gui.load.Names.Spec{2} = 'MM';%re_mm
        end  %re_mm
        if MRSCont.flags.hasRef %Get variables regarding Subspectra
            gui.controls.Number = gui.controls.Number + 1;
            if MRSCont.flags.hasMM%re_mm
                gui.load.Names.Spec{3} = 'reference';%re_mm
            else%re_mm
            gui.load.Names.Spec{2} = 'reference';
            end%re_mm
        if MRSCont.flags.hasWater
            gui.controls.Number = gui.controls.Number + 1;
            if MRSCont.flags.hasMM%re_mm
                gui.load.Names.Spec{4} = 'water';%re_mm
            else%re_mm
            gui.load.Names.Spec{3} = 'water';
            end%re_mm
        end
        else if MRSCont.flags.hasWater
            gui.controls.Number = gui.controls.Number + 1;
            if MRSCont.flags.hasMM%re_mm
                gui.load.Names.Spec{3} = 'water';%re_mm
            else%re_mm
                gui.load.Names.Spec{2} = 'water';
            end%re_mm
            end
        end
        gui.layout.tabs.TabEnables{1} = 'on';
        gui.layout.tabs.Selection  = 1;
        gui.layout.EmptydataPlot = 0;
 %%% 2. CREATING SUB TABS FOR THIS TAB %%%
 % In this case one tab for each subspec (A,B,C,D,ref,water)
            gui.layout.metabLoTab = uix.VBox('Parent', gui.layout.rawTab, 'BackgroundColor',gui.colormap.Background,'Spacing',5);
            gui.layout.rawTab.TabWidth   = 115;
            gui.layout.rawTab.Selection  = 1;
            gui.layout.rawTabhandles = {'metabLoTab'};
            if gui.controls.Number == 1
                    gui.layout.rawTab.TabTitles  = gui.load.Names.Spec;
                    gui.layout.rawTab.TabEnables = {'on'};
            end
            if gui.controls.Number == 2
                if MRSCont.flags.hasRef
                    gui.layout.refLoTab = uix.VBox('Parent', gui.layout.rawTab, 'BackgroundColor',gui.colormap.Background);
                    gui.layout.rawTab.TabTitles  = gui.load.Names.Spec;
                    gui.layout.rawTab.TabEnables = {'on', 'on'};
                    gui.layout.rawTabhandles = {'metabLoTab', 'refLoTab'};
                end
                if MRSCont.flags.hasWater
                    gui.layout.wLoTab = uix.VBox('Parent', gui.layout.rawTab, 'BackgroundColor',gui.colormap.Background);
                    gui.layout.rawTab.TabTitles  = gui.load.Names.Spec;
                    gui.layout.rawTab.TabEnables = {'on', 'on'};
                    gui.layout.rawTabhandles = {'metabLoTab', 'wLoTab'};
                end
            end
            if gui.controls.Number == 3
                gui.layout.refLoTab = uix.VBox('Parent', gui.layout.rawTab,  'BackgroundColor',gui.colormap.Background);
                gui.layout.wLoTab = uix.VBox('Parent', gui.layout.rawTab, 'BackgroundColor',gui.colormap.Background);
                gui.layout.rawTab.TabTitles  = gui.load.Names.Spec;
                gui.layout.rawTab.TabEnables = {'on', 'on','on'};
                gui.layout.rawTabhandles = {'metabLoTab', 'refLoTab', 'wLoTab'};
            end
            if gui.controls.Number == 4 %re_mm
                gui.layout.mmLoTab = uix.VBox('Parent', gui.layout.rawTab,  'BackgroundColor',gui.colormap.Background);%re_mm
                gui.layout.refLoTab = uix.VBox('Parent', gui.layout.rawTab,  'BackgroundColor',gui.colormap.Background);%re_mm
                gui.layout.wLoTab = uix.VBox('Parent', gui.layout.rawTab, 'BackgroundColor',gui.colormap.Background);%re_mm
                gui.layout.rawTab.TabTitles  = gui.load.Names.Spec;%re_mm
                gui.layout.rawTab.TabEnables = {'on', 'on','on','on'};%re_mm
                gui.layout.rawTabhandles = {'metabLoTab', 'mmLoTab', 'refLoTab', 'wLoTab'};%re_mm
            end%re_mm
 %%% 3. FILLING INFO PANEL FOR THIS TAB %%%
 % All the information from the Raw data is read out here
        for t = gui.controls.Number : -1 : 1 % Loop over subspectra & tabs
            % Parameter shown in the info panel on top
            gui.upperBox.data.box = uix.HBox('Parent', gui.layout.(gui.layout.rawTabhandles{t}),'BackgroundColor',gui.colormap.Background, 'Spacing',5);
            if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                gui.upperBox.data.upperLeftButtons = uix.Panel('Parent', gui.upperBox.data.box, ...
                                         'Padding', 5, 'Title', ['Navigate voxel'],...
                                         'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                         'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
                gui.controls.Buttonbox = uix.HBox('Parent',gui.upperBox.data.upperLeftButtons, 'BackgroundColor',gui.colormap.Background);
                gui.controls.navigate_RawTab = uix.Grid('Parent',gui.controls.Buttonbox,'BackgroundColor',gui.colormap.Background);
                gui.controls.text_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','X:',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Y:',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Z:',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.b_left_x = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
                gui.controls.b_left_y = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
                gui.controls.b_left_z = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
                set(gui.controls.b_left_x,'Callback',{@osp_onLeftX,gui});
                set(gui.controls.b_left_y,'Callback',{@osp_onLeftY,gui});
                set(gui.controls.b_left_z,'Callback',{@osp_onLeftZ,gui});
                if gui.info.nXvoxels <= 1
                    gui.controls.b_left_x.Enable = 'off';
                end
                if gui.info.nYvoxels <= 1
                    gui.controls.b_left_y.Enable = 'off';
                end
                if gui.info.nZvoxels <= 1
                    gui.controls.b_left_z.Enable = 'off';
                end
                gui.controls.text_act_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_act_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_act_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.b_right_x = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
                gui.controls.b_right_y = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
                gui.controls.b_right_z = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
                set(gui.controls.b_right_x,'Callback',{@osp_onRightX,gui});
                set(gui.controls.b_right_y,'Callback',{@osp_onRightY,gui});
                set(gui.controls.b_right_z,'Callback',{@osp_onRightZ,gui});   
                if gui.info.nXvoxels <= 1
                    gui.controls.b_right_x.Enable = 'off';
                end
                if gui.info.nYvoxels <= 1
                    gui.controls.b_right_y.Enable = 'off';
                end
                if gui.info.nZvoxels <= 1
                    gui.controls.b_right_z.Enable = 'off';
                end                
                set( gui.controls.navigate_RawTab, 'Widths', [-20 -30 -20 -30], 'Heights', [-33 -33 -33] );
            end
            gui.upperBox.data.Info = uix.Panel('Parent', gui.upperBox.data.box, ...
                                     'Padding', 5, 'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}],...
                                     'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                     'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            gui.upperBox.data.upperButtons = uix.Panel('Parent', gui.upperBox.data.box, ...
                                     'Padding', 5, 'Title', ['Save'],...
                                     'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                     'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            gui.controls.b_save_RawTab = uicontrol('Parent',gui.upperBox.data.upperButtons,'Style','PushButton');
            [img, ~, ~] = imread('Printer.png', 'BackgroundColor', gui.colormap.Background);
            [img2] = imresize(img, 0.10);
            set(gui.controls.b_save_RawTab,'CData', img2, 'TooltipString', 'Create EPS figure from current file');
            set(gui.controls.b_save_RawTab,'Callback',{@osp_onPrint,gui});
            if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                set(gui.upperBox.data.box, 'Width', [-0.12 -0.78 -0.1]);
            else
                set(gui.upperBox.data.box, 'Width', [-0.9 -0.1]);
            end
            % Grid for Plot and Data control sliders
            if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) % HBox for HERMES/HERCULES
                gui.Plot.data = uix.HBox('Parent', gui.layout.(gui.layout.rawTabhandles{t}), 'BackgroundColor',gui.colormap.Background, 'Units', 'normalized');
            else
                gui.Plot.data = uix.VBox('Parent', gui.layout.(gui.layout.rawTabhandles{t}), 'BackgroundColor',gui.colormap.Background, 'Units', 'normalized');
            end
            gui.InfoText.data  = uicontrol('Parent',gui.upperBox.data.Info,'style','text',...
                                          'FontSize', 12, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
            set(gui.layout.(gui.layout.rawTabhandles{t}), 'Heights', [-0.1 -0.9]);


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
        if ~isfield(MRSCont.flags,'isPRIAM') && ~isfield(MRSCont.flags,'isMRSI') && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
           set(gui.InfoText.data, 'String', sprintf(StatText));
        else
            StatText = ['Voxel ' num2str(gui.controls.act_x) ': '  StatText];
            set(gui.InfoText.data, 'String', sprintf(StatText));
        end
        
 %%% 4. VISUALIZATION PART OF THIS TAB %%%
 %osp_plotLoad is used to visualize the raw data. Number of subplots
 %depends on the number of subspectra of the seuqence
        temp = figure( 'Visible', 'off' );
        if t == 1 %Metabolite data/tab
            temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mets');
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
            if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) %Four windows for HERMES/HERCULES
                gui.layout.multiACload = uix.VBox('Parent', gui.Plot.data, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.layout.multiAload = uix.VBox('Parent', gui.layout.multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    gui.layout.multiCload = uix.VBox('Parent', gui.layout.multiACload,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                gui.layout.multiBDload = uix.VBox('Parent', gui.Plot.data,'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.layout.multiBload = uix.VBox('Parent', gui.layout.multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    gui.layout.multiDload = uix.VBox('Parent', gui.layout.multiBDload, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                set( temp.Children(1), 'Parent', gui.layout.multiDload );
                set( temp.Children(1), 'Parent', gui.layout.multiCload );
                set( temp.Children(1), 'Parent', gui.layout.multiBload );
                set( temp.Children(1), 'Parent', gui.layout.multiAload );
                set(gui.Plot.data,'Width', [-0.49 -0.49]);
                set(gui.layout.multiDload.Children(1), 'Units', 'normalized')
                set(gui.layout.multiDload.Children(1), 'OuterPosition', [0,0,1,1])
                set(gui.layout.multiCload.Children(1), 'Units', 'normalized')
                set(gui.layout.multiCload.Children(1), 'OuterPosition', [0,0,1,1])
                set(gui.layout.multiBload.Children(1), 'Units', 'normalized')
                set(gui.layout.multiBload.Children(1), 'OuterPosition', [0,0,1,1])
                set(gui.layout.multiAload.Children(1), 'Units', 'normalized')
                set(gui.layout.multiAload.Children(1), 'OuterPosition', [0,0,1,1])
                
            end
        else
            if MRSCont.flags.hasMM %re_mm
                if t == 2 %ref data/tab %re_mm
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm'); %re_mm
                    ViewAxes = gca(); %re_mm
                    set( ViewAxes, 'Parent', gui.Plot.data ); %re_mm
                end %re_mm
                if t == 3 %ref data/tab %re_mm
                    if MRSCont.flags.hasRef%re_mm
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref'); %re_mm
                    else%re_mm
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w'); %re_mm
                    end%re_mm
                    ViewAxes = gca(); %re_mm
                    set( ViewAxes, 'Parent', gui.Plot.data ); %re_mm
                end %re_mm
                if t == 4 %ref data/tab %re_mm
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w'); %re_mm
                    ViewAxes = gca(); %re_mm
                    set( ViewAxes, 'Parent', gui.Plot.data ); %re_mm
                end %re_mm
            else %re_mm
                if t == 2 %ref data/tab
                    if MRSCont.flags.hasRef
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref');
                        ViewAxes = gca();
                        set( ViewAxes, 'Parent', gui.Plot.data );
                    elseif MRSCont.flags.hasWater
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w');
                        ViewAxes = gca();
                        set( ViewAxes, 'Parent', gui.Plot.data );
                    end
                else %water data/tab has only one window all the time
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w');
                    ViewAxes = gca();
                    set(ViewAxes, 'Parent', gui.Plot.data );
                end
            end %re_mm
        end

        % Get rid of the Load figure
        close( temp );
        end
        setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end