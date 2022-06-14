function osp_iniCoregWindow(gui)
%% osp_iniCoregWindow
%   This function creates the initial coregistration window in the gui.
%
%
%   USAGE:
%       osp_iniCoregWindow(gui);
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
% This function creates the initial coreg/seg window
    warning('off','MATLAB:handle_graphics:exceptions:SceneNode');
    MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class

    % Get variables regarding secondary T1 and PET images
    if MRSCont.flags.hasSecondT1
        gui.controls.NumberImages = gui.controls.NumberImages + 1;
        gui.load.Names.Images{2} = '2nd structural';
        if MRSCont.flags.hasPET
            gui.controls.NumberImages = gui.controls.NumberImages + 1;
            gui.load.Names.Images{3} = 'PET';
        end
    else if MRSCont.flags.hasPET
            gui.controls.NumberImages = gui.controls.NumberImages + 1;
            gui.load.Names.Images{2} = 'PET';
        end
    end

    gui.layout.tabs.TabEnables{4} = 'on';
    gui.layout.tabs.Selection  = 4;
    gui.layout.EmptyQuantPlot = 0;

    %%% 2. CREATING SUB TABS FOR THIS TAB %%%
    % In this case one tab for each image (first structural, second
    % structural, PET)
    gui.layout.structuralLoTab = uix.VBox('Parent', gui.layout.coregTab, 'BackgroundColor',gui.colormap.Background,'Spacing',5);
    gui.layout.coregTab.TabWidth   = 115;
    gui.layout.coregTab.Selection  = 1;
    gui.layout.coregTabhandles = {'structuralLoTab'};
    if gui.controls.NumberImages == 1
        gui.layout.coregTab.TabTitles  = gui.load.Names.Images;
        gui.layout.coregTab.TabEnables = {'on'};
    end
    if gui.controls.NumberImages == 2
        if MRSCont.flags.hasRef
            gui.layout.secondstructuralLoTab = uix.VBox('Parent', gui.layout.coregTab, 'BackgroundColor',gui.colormap.Background);
            gui.layout.coregTab.TabTitles  = gui.load.Names.Images;
            gui.layout.coregTab.TabEnables = {'on', 'on'};
            gui.layout.coregTabhandles = {'structuralLoTab', 'secondstructuralLoTab'};
        end
        if MRSCont.flags.hasWater
            gui.layout.petLoTab = uix.VBox('Parent', gui.layout.coregTab, 'BackgroundColor',gui.colormap.Background);
            gui.layout.coregTab.TabTitles  = gui.load.Names.Images;
            gui.layout.coregTab.TabEnables = {'on', 'on'};
            gui.layout.coregTabhandles = {'structuralLoTab', 'petLoTab'};
        end
    end
    if gui.controls.NumberImages == 3
        gui.layout.secondstructuralLoTab = uix.VBox('Parent', gui.layout.coregTab,  'BackgroundColor',gui.colormap.Background);
        gui.layout.petLoTab = uix.VBox('Parent', gui.layout.coregTab, 'BackgroundColor',gui.colormap.Background);
        gui.layout.coregTab.TabTitles  = gui.load.Names.Images;
        gui.layout.coregTab.TabEnables = {'on', 'on','on'};
        gui.layout.coregTabhandles = {'structuralLoTab', 'secondstructuralLoTab', 'petLoTab'};
    end

%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
for t = gui.controls.NumberImages : -1 : 1 % Loop over tabs
    gui.upperBox.coreg.box = uix.HBox('Parent', gui.layout.(gui.layout.coregTabhandles{t}),'BackgroundColor',gui.colormap.Background,'Spacing',5);
    if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                gui.upperBox.coreg.upperLeftButtons = uix.Panel('Parent', gui.upperBox.coreg.box, ...
                                         'Padding', 5, 'Title', ['Navigate voxel'],...
                                         'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                         'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
                gui.controls.Buttonbox = uix.HBox('Parent',gui.upperBox.coreg.upperLeftButtons, 'BackgroundColor',gui.colormap.Background);
                gui.controls.navigate_RawTab = uix.Grid('Parent',gui.controls.Buttonbox,'BackgroundColor',gui.colormap.Background);
                gui.controls.text_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','X:',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Y:',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Z:',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
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
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_act_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_act_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
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
    gui.upperBox.coreg.upperButtons1 = uix.Panel('Parent', gui.upperBox.coreg.box, ...
                                     'Padding', 5, 'Title', ['Nii Viewer'],...
                                     'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                     'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
    gui.controls.b_nii_coregTab = uicontrol('Parent',gui.upperBox.coreg.upperButtons1,'Style','PushButton');
    [img, ~, ~] = imread('Nii.png', 'BackgroundColor', gui.colormap.Background);
    [img2] = imresize(img, 0.1);
    set(gui.controls.b_nii_coregTab,'CData', img2, 'TooltipString', 'Start external Nifti viewer');
    set(gui.controls.b_nii_coregTab,'Callback',{@osp_onNii,gui});
    gui.upperBox.coreg.Info = uix.Panel('Parent', gui.upperBox.coreg.box, 'Padding', 5, ...
                              'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}],...
                              'FontName', gui.font,'HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,...
                              'ForegroundColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
    gui.upperBox.coreg.upperButtons = uix.Panel('Parent', gui.upperBox.coreg.box, ...
                                     'Padding', 5, 'Title', ['Save'],...
                                     'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                     'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
    gui.controls.b_save_coregTab = uicontrol('Parent',gui.upperBox.coreg.upperButtons,'Style','PushButton');
    [img, ~, ~] = imread('Printer.png', 'BackgroundColor', gui.colormap.Background);
    [img2] = imresize(img, 0.1);
    set(gui.controls.b_save_coregTab,'CData', img2, 'TooltipString', 'Create EPS figure from current file');
    set(gui.controls.b_save_coregTab,'Callback',{@osp_onPrint,gui});
    if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
        set(gui.upperBox.coreg.box, 'Width', [-0.12 -0.1 -0.68 -0.1]);
    else
        set(gui.upperBox.coreg.box, 'Width', [-0.1 -0.8 -0.1]);
    end
    % Creates layout for plotting and data control
    gui.Plot.coreg = uix.HBox('Parent', gui.layout.(gui.layout.coregTabhandles{t}),'BackgroundColor',gui.colormap.Background);
    set(gui.layout.(gui.layout.coregTabhandles{t}), 'Heights', [-0.1 -0.9]);
    % Get parameter from file to fill the info panel

        StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.controls.Selected}.spectralwidth) ' Hz'...
            '\nraw subspecs: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
            '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
            num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];

   gui.InfoText.coreg  = uicontrol('Parent',gui.upperBox.coreg.Info,'style','text',...
        'FontSize', 12, 'FontName', gui.font,'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

%%% 2.VISUALIZATION PART OF THIS TAB %%%
% In this case osp_plotCoreg or osp_plotSegment is used to visualize the
% coregistration or the segmentation
    gui.Results.coreg = uix.VBox('Parent', gui.Plot.coreg,'BackgroundColor',gui.colormap.Background);
    temp = figure( 'Visible', 'off' );
    if (isfield(MRSCont.flags, 'isPRIAM')) &&  MRSCont.flags.isPRIAM
        if MRSCont.flags.didSeg %Did segment. In this case coreg has already been performed. Visualize both
            osp_plotCoreg(MRSCont, gui.controls.Selected,gui.controls.act_x);
            ViewAxes = gca();
            set(ViewAxes, 'Parent', gui.Results.coreg );
            colormap(gui.Results.coreg.Children,'gray')
            close( temp );
            temp = figure( 'Visible', 'off' );
            osp_plotSegment(MRSCont, gui.controls.Selected,gui.controls.act_x);
            ViewAxes = gca();
            set(ViewAxes, 'Parent', gui.Results.coreg );
            colormap(gui.Results.coreg.Children(1),'gray');
            close( temp );
        else % Only coreg has been run
            osp_plotCoreg(MRSCont, gui.controls.Selected,gui.controls.act_x);
            ViewAxes = gca();
            set(ViewAxes, 'Parent', gui.Results.coreg );
            colormap(gui.Results.coreg.Children,'gray');
            close( temp );
        end

    else
        if t == 1 %First structural tab
            if MRSCont.flags.didSeg %Did segment. In this case coreg has already been performed. Visualize both
                osp_plotCoreg(MRSCont, gui.controls.Selected);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', gui.Results.coreg );
                colormap(gui.Results.coreg.Children,'gray')
                close( temp );
                temp = figure( 'Visible', 'off' );
                osp_plotSegment(MRSCont, gui.controls.Selected);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', gui.Results.coreg );
                colormap(gui.Results.coreg.Children(1),'gray');
                close( temp );
            else % Only coreg has been run
                osp_plotCoreg(MRSCont, gui.controls.Selected);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', gui.Results.coreg );
                colormap(gui.Results.coreg.Children,'gray');
                close( temp );
            end
        elseif t == 2 % second structural tab
            % If voxel has been registered to a second T1, add it here
            if MRSCont.flags.didSeg && MRSCont.flags.hasSecondT1
                temp = figure( 'Visible', 'off' );
                osp_plotCoregSecond(MRSCont, gui.controls.Selected);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', gui.Results.coreg );
                colormap(gui.Results.coreg.Children,'gray')
                close( temp );
                temp = figure( 'Visible', 'off' );
                osp_plotSegmentSecond(MRSCont, gui.controls.Selected);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', gui.Results.coreg );
                colormap(gui.Results.coreg.Children(1),'gray');
                close( temp );
            else
                temp = figure( 'Visible', 'off' );
                osp_plotCoregSecond(MRSCont, gui.controls.Selected);
                ViewAxes = gca();
                set(ViewAxes, 'Parent', gui.Results.coreg );
                colormap(gui.Results.coreg.Children(1),'gray');
                close( temp );
            end
        elseif t == 3 % PET tab
            % If voxel has been registered to a PET image, add it here
            if MRSCont.flags.hasPET
                temp = figure( 'Visible', 'off' );
                temp = osp_plotCoregPET(MRSCont, gui.controls.Selected);
                set( temp.Children(2), 'Parent', gui.Results.coreg );
                set( temp.Children(1), 'Parent', gui.Results.coreg );
                set(gui.Results.coreg,'Heights', [-0.49 -0.49]);
                set(gui.Results.coreg.Children(2), 'Units', 'normalized')
                set(gui.Results.coreg.Children(2), 'OuterPosition', [0,0.5,1,0.5])
                set(gui.Results.coreg.Children(1), 'Units', 'normalized')
                set(gui.Results.coreg.Children(1), 'OuterPosition', [0,0,1,0.5])
                colormap(gui.Results.coreg.Children(2),'gray');
                close( temp );
            end
        end
    end
    h = findall(groot,'Type','figure');
    for ff = 1 : length(h)
        if ~(strcmp(h(ff).Tag, 'Osprey') ||  strcmp(h(ff).Tag, 'TMWWaitbar'))
            close(h(ff))
        end
    end
    warning('on','MATLAB:handle_graphics:exceptions:SceneNode');
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end
