function OspreyGUI(MRSCont)
%% OspreyGUI(MRSCont)
%   This function creates a one-in-all figure with visualizations of the
%   processed data (spectra in the frequency domain), voxel coregistration
%   and segmentation, and quantification tables.
%
%   The figure contains several tabs, not all of which may be available at
%   all times:
%       - Raw Data
%       - Processed Data    
%       - Coregistration and segmentation
%       - Fit
%       - Quantification
%
%   As an example, if coregistration and segmentation have not been
%   performed, the respective tab will be grayed out.
%
%   USAGE:
%       OspreyGUI(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-06-30)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-06-30: First version of the code.

% Close any remaining open figures
close all;


% Create the user interface for the application and return a
% structure of handles for global use.
gui = struct();
gui.SelectedFile = 1;
% Create the overall figure
gui.window = figure(...
    'Name', 'Osprey', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'HandleVisibility', 'off', ...
    'Renderer', 'painters');
% Resize such that width is 1.2941 * height (1.2941 is the ratio
% between width and height of standard US letter size (11x8.5 in).
screenSize      = get(0,'ScreenSize');
canvasSize      = screenSize;
canvasSize(4)   = screenSize(4) * 0.9;
canvasSize(3)   = canvasSize(4) * (11/8.5);
canvasSize(2)   = (screenSize(4) - canvasSize(4))/2;
canvasSize(1)   = (screenSize(3) - canvasSize(3))/2;
set(gui.window, 'Position', canvasSize);

% Create the main horizontal box division between menu (left) and display
% panel tabs (right)
mainLayout = uix.HBox(...
    'Parent', gui.window);
    % Create the left-side menu
    leftMenu = uix.VBox(...
        'Parent',mainLayout, ...
        'Padding',5);
        % Divide into the upper panel containing the Gannet logo button and the
        % lower panel containing the menu buttons
        p1 = uix.Panel(...
            'Parent', leftMenu, ...
            'Padding', 5);
            % Place Gannet logo on a push button
            b_about = uicontrol(p1,'Style','PushButton');
            set(b_about,'Units','Normalized','Position',[0 0 1 1],'BackgroundColor',[0.93 0.93 0.93]);
            [img, ~, ~] = imread('osprey.png', 'BackgroundColor', [0.93 0.93 0.93]);
            [img2] = imresize(img, 0.09);
            set(b_about,'CData', img2);
            logoFcn = @()imread('osprey.png', 'BackgroundColor', [0.93 0.93 0.93]);
logoBanner = uiw.utility.loadIcon(logoFcn);

d = uiw.dialog.About(...
    'Name', 'Osprey',...
    'Version','0.0.1',...
    'Date', 'October 6, 2019',...
    'Timeout', 3,...
    'CustomText', 'Osprey is provided by Johns Hopkins University.',...
    'ContactInfo', 'gabamrs@gmail.com',...
    'LogoCData', logoBanner);

        p2 = uix.VButtonBox(...
            'Parent', leftMenu, ...
            'Padding', 5, ...
            'Spacing', 10, ...
            'ButtonSize', [180 50]);
        set(leftMenu, 'Heights', [-0.2 -0.8]);
        % Data input button
        b_input = uicontrol('Parent', p2,'Style','PushButton','String','Data input');
        set(b_input,'Units','Normalized','Position',[0.1 0.9 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        set(b_input,'Callback',@gui_input);
        % Load button
        b_load = uicontrol('Parent', p2,'Style','PushButton','String','Load data','Enable','off');
        set(b_load,'Units','Normalized','Position',[0.1 0.75 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % Fit button
        b_fit = uicontrol('Parent', p2,'Style','PushButton','String','Fit data','Enable','off');
        set(b_fit,'Units','Normalized','Position',[0.1 0.67 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % Coregister button
        b_coreg = uicontrol('Parent', p2,'Style','PushButton','String','CoRegister','Enable','off');
        set(b_coreg,'Units','Normalized','Position',[0.1 0.59 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % Segment button
        b_segm = uicontrol('Parent', p2,'Style','PushButton','String','Segment','Enable','off');
        set(b_segm,'Units','Normalized','Position',[0.1 0.51 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % Quantify button
        b_quant = uicontrol('Parent', p2,'Style','PushButton','String','Quantify','Enable','off');
        set(b_quant,'Units','Normalized','Position',[0.1 0.43 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % CoRegStandAlone button
        b_coregSA = uicontrol('Parent', p2,'Style','PushButton','String','CoRegStandAlone','Enable','off');
        set(b_coregSA,'Units','Normalized','Position',[0.1 0.28 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        % DeIdentify button
        b_deid = uicontrol('Parent', p2,'Style','PushButton','String','DeIdentify','Enable','off');
        set(b_deid,'Units','Normalized','Position',[0.1 0.2 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        set(b_deid,'Callback',@gui_deid);
        % Settings button
        b_settings = uicontrol('Parent', p2,'Style','PushButton','String','Settings');
        set(b_settings,'Units','Normalized','Position',[0.1 0.05 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        set(b_settings,'Callback',@gui_settings);
        % Exit button
        b_exit = uicontrol('Parent', p2,'Style','PushButton','String','Exit');
        set(b_exit,'Units','Normalized','Position',[0.1 0 0.8 0.08], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
        set(b_exit,'Callback',@onExit);
        
    % Create list of files
    controlPanel = uix.Panel('Parent', leftMenu, 'Title', 'Select a File:');
    set(controlPanel,'Units','Normalized','Position',[0.5 0 0.66 0.1], 'FontSize', 16, 'FontName', 'Century Gothic', 'FontWeight', 'Bold');
    fileList = MRSCont.files;
    SepFileList = cell(1,MRSCont.nDatasets);
    RedFileList = cell(1,MRSCont.nDatasets);
    for i = 1 : MRSCont.nDatasets
        SepFileList{i} =  split(fileList(i), filesep);
        RedFileList{i} = [filesep SepFileList{i}{end-2} filesep SepFileList{i}{end-1} filesep SepFileList{i}{end}];        
    end
    clear SepFileList 
    ListBox = uicontrol( 'Style', 'list', ...
    'BackgroundColor', 'w', ...
    'Parent', controlPanel, ...
    'String',RedFileList(:) , ...
    'Value', gui.SelectedFile, ...
    'Callback', @onListSelection);
    
    % Create the display panel tab row
    tabs = uix.TabPanel(...
        'Parent', mainLayout, ...
        'Padding', 5, ...
        'FontName', 'Century Gothic', ...
        'FontSize', 16);
    rawTab     = uix.VBox('Parent', tabs, 'Padding', 5);
    proTab     = uix.VBox('Parent', tabs, 'Padding', 5);
    coregTab    = uix.Grid('Parent', tabs, 'Spacing', 5);
    fitTab      = uix.VBox('Parent', tabs, 'Padding', 5);
    quantifyTab = uix.Panel('Parent', tabs, 'Padding', 5);
    tabs.TabTitles  = {'Raw', 'Processed', 'Coreg / Segment', 'Fit', 'Quantification'};
    tabs.TabWidth   = 200;
    tabs.Selection  = 1;
    tabs.TabEnables = {'off', 'off', 'off', 'off', 'off'};
    
    set( mainLayout, 'Widths', [-0.12  -0.88], 'Spacing', 5 );
    % Now enable the display tabs depending on which processing steps have
    % been completed:
    if MRSCont.flags.didLoadData % Was data loaded at all that can be looked at?
        tabs.TabEnables{1} = 'on';
        initialLoadWindow()
    
    end
    
    if MRSCont.flags.didCoreg % Have coreg/segment masks been created?
        tabs.TabEnables{2} = 'on';
        
        %%% COREG+SEGMENT DISPLAY PANEL %%%
        % Fill the coreg display panel
        
        S.sl_plotCoreg = uicontrol('Parent',coregTab,'style','slide',...
            'min',1,'max',MRSCont.nDatasets,'value',1,...
            'sliderstep',[1 10],...
            'callback',{@sl_plotCoreg_Call,S});
        if MRSCont.nDatasets == 1
            S.sl_plotCoreg.Visible = 'off';
        end
        set(coregTab, 'Widths',[-1 -1 10], 'Heights', [-1 -1]);
        %%% /COREG+SEGMENT DISPLAY PANEL %%%
    
    end
    if MRSCont.flags.didFit % Has data fitting been run?
        
        %%% FIT DISPLAY PANEL %%%
        tabs.TabEnables{3} = 'on';
        % Fill the fit display panel
        fitPlot = uix.HBox(...
            'Parent', fitTab, ...
            'Padding', 5);
        
        gui.dataAxes = axes('Parent', fitPlot);  
        % Create the controls
        S.fitControls = uix.HBox(...
            'Parent', fitPlot, ...
            'Padding', 5);
        % Create buttons to switch individual plots on/off
        S.fitButtonList = uix.VButtonBox(...
            'Parent', S.fitControls, ...
            'Padding', 5);
        for kk = 1:length(MRSCont.fit.basisSet.name)
            S.fitMetabButton{kk} = uicontrol( 'Parent', S.fitButtonList, 'String', MRSCont.fit.basisSet.name{kk}, 'Style', 'checkbox', 'FontSize', 12, 'FontName', 'Century Gothic', 'Value', 1);
        end
        set(S.fitButtonList, 'ButtonSize', [100 15], 'Spacing', 10);
        
        % Slider to select datasets on the far right
        S.sl_plotFit = uicontrol('Parent',S.fitControls,'style','slide',...
            'min',1,'max',MRSCont.nDatasets,'value',1,...
            'Units', 'Normalized', 'Position', [0 0 1 1], ...
            'sliderstep',[1 1],...
            'callback',{@sl_plotFit_Call,S});
        set(S.fitControls, 'Widths', [-1 10]);
        
        if MRSCont.nDatasets == 1
            S.sl_plotFit.Visible = 'off';
        end
        
        set(fitPlot, 'Widths', [-1 -0.2]);
        %%% /FIT DISPLAY PANEL %%%
        
        
        tabs.TabEnables{4} = 'on';
        
        
    end
    
set(mainLayout, 'Widths', [-0.2 -0.8]);
    

function initialLoadWindow()
  
        %%% DATA DISPLAY PANEL %%%
        % Fill the data display panel
        gui.dataStats = uix.Panel(...
            'Parent', rawTab, ...
            'Padding', 5);
        gui.dataPlot = uix.HBox(...
            'Parent', rawTab, ...
            'Padding', 5);
        set(rawTab, 'Heights', [-0.1 -0.9]);


        temp = figure( 'Visible', 'off' );
        temp = osp_plotLoad(MRSCont, gui.SelectedFile,'mets',1 );
        gui.ViewAxes = gca();
        set( gui.ViewAxes, 'Parent', gui.dataPlot );
        % Get rid of the Load figure
        close( temp );
        gui.dataControls = uix.HBox(...
            'Parent', gui.dataPlot, ...
            'Padding', 5);
        gui.sl_plotData = uicontrol('Parent',gui.dataControls,'style','slide',...
            'min',1,'max',MRSCont.nDatasets,'value',1,...
            'Units', 'Normalized', 'Position', [0 0 1 1], ...
            'sliderstep',[1 1],...
            'callback',{@sl_plotData_Call,gui});
        addlistener(gui.sl_plotData, 'Value', 'PreSet',@sl_plotData_realtime_Call);

        if MRSCont.nDatasets == 1
            gui.sl_plotData.Visible = 'off';
        end
        set(gui.dataPlot, 'Widths', [-0.97 -0.03]);
        %%% /DATA DISPLAY PANEL %%%
end

function updateLoadWindow()
        if ishandle( gui.dataPlot )
            delete( gui.dataPlot );
        end
        if ishandle( gui.dataStats )
            delete( gui.dataStats );
        end    
        %%% DATA DISPLAY PANEL %%%
        % Fill the data display panel
        gui.dataStats = uix.Panel(...
            'Parent', rawTab, ...
            'Padding', 5);
        gui.dataPlot = uix.HBox(...
            'Parent', rawTab, ...
            'Padding', 5);
        set(rawTab, 'Heights', [-0.1 -0.9]);


        temp = figure( 'Visible', 'off' );
        temp = osp_plotLoad(MRSCont, gui.SelectedFile,'mets',1 );
        gui.ViewAxes = gca();
        set( gui.ViewAxes, 'Parent', gui.dataPlot );
        % Get rid of the Load figure
        close( temp );
        gui.dataControls = uix.HBox(...
            'Parent', gui.dataPlot, ...
            'Padding', 5);
        gui.sl_plotData = uicontrol('Parent',gui.dataControls,'style','slide',...
            'min',1,'max',MRSCont.nDatasets,'value',gui.SelectedFile,...
            'Units', 'Normalized', 'Position', [0 0 1 1], ...
            'sliderstep',[1 1],...
            'callback',{@sl_plotData_Call,gui});
        addlistener(gui.sl_plotData,'ContinuousValueChange',@sl_plotData_realtime_Call);
%         addlistener(gui.sl_plotData, 'Value', 'PreSet',@sl_plotData_realtime_Call);
        if MRSCont.nDatasets == 1
            gui.sl_plotData.Visible = 'off';
        end
        set(gui.dataPlot, 'Widths', [-0.97 -0.03]);
        %%% /DATA DISPLAY PANEL %%%
end

%-------------------------------------------------------------------------%
    %------- CALLBACK FUNCTIONS FOR THE MAIN MENU -------% 
    function onExit( ~, ~ )
        % User wants to quit out of the application
        delete( gui.window );
    end % onExit
    %------- /CALLBACK FUNCTIONS FOR THE MAIN MENU ------% 
    %------- CALLBACK FUNCTIONS FOR THE LIST MENU -------%   

    function onListSelection( src, ~ )
        % User selected a demo from the list - update "data" and refresh
        gui.SelectedFile = get( src, 'Value' );
        set(gui.sl_plotData,'value', gui.SelectedFile);
        set(ListBox, 'value', gui.SelectedFile);
        updateLoadWindow()
    end % onListSelection
     %------- CALLBACK FUNCTIONS FOR THE LIST MENU -------%   
     
    function [] = sl_plotFit_Call(varargin)
        % Callback for the slider.
        [h,S] = varargin{[1,3]};  % calling handle and data structure.
        idx = get(h,'value');

    end

 function sl_plotData_realtime_Call(~,~)
        % Callback for the slider.
        slider_value=get(gui.sl_plotData,'value');
        idx=round(slider_value);
        set(ListBox, 'value', idx);
    end


 function [h,S] = sl_plotData_Call(varargin)
        % Callback for the slider.
        [h,S] = varargin{[1,3]};  % calling handle and data structure.
        idx=round(h.Value);
        h.Value=idx;
        gui.SelectedFile = idx;
        set(ListBox, 'value', idx);
        updateLoadWindow()
        
    end
end

