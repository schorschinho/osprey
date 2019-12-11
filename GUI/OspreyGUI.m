function OspreyGUI(MRSCont)
%% OspreyGUI(MRSCont)
%   This function creates a one-in-all figure with visualizations of the
%   processed data (spectra in the frequency domain), voxel coregistration
%   and segmentation, and quantification tables.
%
%   The figure contains several tabs, not all of which may be available at
%   all times:
%       - Data
%       - Processed
%       - Coregistration and segmentation
%       - Fit
%       - Quantification
%       - Overview
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
%   AUTHORS:
%       Helge Zöllner (Johns Hopkins University, 2019-11-07)
%       hzoelln2@jhmi.edu
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
%       2019-07-11: First version of the code.
%



 
%% Initialize global variables

    % Close any remaining open figures & specifiy spmfolder
    close all;
    [settingsFolder,~,~] = fileparts(which('OspreySettings.m'));
    allFolders      = strsplit(settingsFolder, filesep);
    ospFolder       = strjoin(allFolders(1:end-1), filesep); % parent folder (= Osprey folder)
    matlabFolder    = strjoin(allFolders(1:end-2), filesep); % parent-parent folder (usually MATLAB folder)
    % SPM
    addpath(genpath([ospFolder filesep 'spm12' filesep]));    % SPM path
    % Check if SPM12 is installed
    spmversion = fileparts(which('spm'));
    if isempty(spmversion)
        error('SPM not found! Please install SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12) and set the path in OspreySettings.');
    elseif strcmpi(spmversion(end-3:end),'spm8')
        error(['SPM8 detected, but only SPM12 is supported. ' ...
               'Please install SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12) and set the path in OspreySettings.']);
    end

    % Create the user interface for the application and return a
    % structure of handles and variables for global use.
    gui = struct();
    %Here we set up the color layout
    %[0 0 0] Black
    %[51/255 51/255 51/255] DarkGray
    %[71/255 71/255 71/255] LightGray
    %[110/255 136/255 164/255] LightBlue
    %[11/255 71/255 111/255] DarkBlue
    %[244/255 244/255 242/255] OffWhite
    %[217/255 224/255 33/255] Yellow
    %[1 1 1] White

    %Darkmode colormmap
%     gui.colormap.Background = [71/255 71/255 71/255];
%     gui.colormap.LightAccent = [110/255 136/255 164/255];
%     gui.colormap.Foreground = [244/255 244/255 242/255];
%     gui.colormap.Accent = [217/255 224/255 33/255];
    %Bluemode colormap
%     gui.colormap.Background = [110/255 136/255 164/255];
%     gui.colormap.LightAccent = [51/255 51/255 51/255];
%     gui.colormap.Foreground = [244/255 244/255 242/255];
%     gui.colormap.Accent = [217/255 224/255 33/255];
    %Bluemode colormap
    gui.colormap.Background = [255/255 254/255 254/255];
    gui.colormap.LightAccent = [110/255 136/255 164/255];
    gui.colormap.Foreground = [11/255 71/255 111/255];
    gui.colormap.Accent = [11/255 71/255 111/255];
    
    
    MRSCont.colormap = gui.colormap;

    gui.SelectedDataset = 1;
    gui.SelectedSubFile = 1;
    gui.SelectedSubSpec = 1;
    gui.SelectedFit = 1;
    gui.SelectedModel = 1;
    gui.SelectedQuant = 1;
    gui.NoSpecs = 1;
    gui.SpecNames = {'metabolites'};
    gui.SelectedMetab = 1;
    gui.KeyPress = 0;
    if MRSCont.flags.didLoadData
        if ~isempty(MRSCont.raw{1,gui.SelectedDataset}.seq)
            if strcmp(sprintf('\n'),MRSCont.raw{1,gui.SelectedDataset}.seq(end)) %Clean up Sequence Name if needed
                SeqName = MRSCont.raw{1,gui.SelectedDataset}.seq(1:end-1);
            else
                SeqName = MRSCont.raw{1,gui.SelectedDataset}.seq;
            end
        else
            SeqName ='';
        end
        gui.GeometryNames = fieldnames(MRSCont.raw{1,1}.geometry.size); %Get variables regarding voxel geometry
    end
    if MRSCont.flags.didProcess %Get variables regarding the processing
        gui.NoPro = length(fieldnames(MRSCont.processed));
        gui.ProNames = fieldnames(MRSCont.processed);
        gui.ProNames = sort(gui.ProNames);
    end
    if MRSCont.flags.didFit %Get variables regarding the fitting
        gui.NoFits = length(fieldnames(MRSCont.fit.results));
        gui.FitNames = fieldnames(MRSCont.fit.results);
    end

    if MRSCont.flags.didQuantify %Get variables regarding the quantification
        gui.NoQuantModels = length(fieldnames(MRSCont.quantify.tables));
        gui.QuantModelNames = fieldnames(MRSCont.quantify.tables);
        gui.NoQuants = length(fieldnames(MRSCont.quantify.tables.(gui.QuantModelNames{1})));
        gui.QuantNames = fieldnames(MRSCont.quantify.tables.(gui.QuantModelNames{1}));
        gui.NoQuantMetabs = length(MRSCont.quantify.metabs);
    end
    if MRSCont.flags.didOverview %Get variables for the overview tab
        gui.NAAnormed = 1;
        gui.NoGroups = MRSCont.overview.NoGroups;
        [cb] = cbrewer('qual', 'Dark2', 12, 'pchip');
        temp = cb(3,:);
        cb(3,:) = cb(4,:);
        cb(4,:) = temp;
        if isfield(MRSCont.overview, 'corr')
            gui.CorrNames = MRSCont.overview.corr.Names{1};
            gui.CorrMeas = MRSCont.overview.corr.Meas;
        end
        gui.SelectedCorr = 1;
        gui.SelectedCorrChoice = 1;
        gui.QMNames = {'SNR','FWHM (ppm)'};
    end
%% Create the overall figure
gui.window = figure('Name', 'Osprey', 'NumberTitle', 'off', 'MenuBar', 'none', 'Visible', 'on',...
                    'ToolBar', 'none', 'HandleVisibility', 'off', 'Renderer', 'painters', 'Color', gui.colormap.Background);
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
    gui.mainLayout = uix.HBox('Parent', gui.window,'BackgroundColor',gui.colormap.Background);
    % Create the left-side menu
    gui.leftMenu = uix.VBox('Parent',gui.mainLayout, 'Padding',5,'BackgroundColor',gui.colormap.Background);
    % Divide into the upper panel containing the Gannet logo button and the
    % lower panel containing the menu buttons
   gui.b_about = uicontrol(gui.leftMenu,'Style','PushButton');
        set(gui.b_about,'Units','Normalized','Position',[0 0 1 1],'BackgroundColor',gui.colormap.Background);
        [img, ~, ~] = imread('osprey.png', 'BackgroundColor', gui.colormap.Background);
        [img2] = imresize(img, 0.09);
        set(gui.b_about,'CData', img2);

%% Create left menu
    gui.p2 = uix.VButtonBox(...
            'Parent', gui.leftMenu, 'Padding', 5, 'Spacing', 5, ...
            'BackgroundColor',gui.colormap.Background);
    set(gui.leftMenu, 'Heights', [-0.2 -0.8]);
    set(gui.p2, 'ButtonSize', [300 60]);
    % Load button
    gui.b_load = uicontrol('Parent', gui.p2,'Style','PushButton','String','Load data','Enable','on','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.b_load,'Units','Normalized','Position',[0.1 0.75 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.b_load,'Callback',{@onLoad});
    if MRSCont.flags.didLoadData
        gui.b_load.Enable = 'off';
    end
    % Process button
    gui.b_proc = uicontrol('Parent', gui.p2,'Style','PushButton','String','Process data','Enable','on','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.b_proc,'Units','Normalized','Position',[0.1 0.75 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.b_proc,'Callback',{@onProc});
    if MRSCont.flags.didProcess
        gui.b_proc.Enable = 'off';
    else if ~MRSCont.flags.didLoadData
            gui.b_proc.Enable = 'off';
        end 
    end
    % Fit button
    gui.b_fit = uicontrol('Parent', gui.p2,'Style','PushButton','String','Fit data','Enable','on','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.b_fit,'Units','Normalized','Position',[0.1 0.67 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.b_fit,'Callback',{@onFit});
    if MRSCont.flags.didFit
        gui.b_fit.Enable = 'off';
    else if ~MRSCont.flags.didProcess
            gui.b_fit.Enable = 'off';
        end
    end
    % Coregister button
    gui.b_coreg = uicontrol('Parent', gui.p2,'Style','PushButton','String','CoRegister','Enable','off','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.b_coreg,'Units','Normalized','Position',[0.1 0.59 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.b_coreg,'Callback',{@onCoreg});
    if ~isempty(MRSCont.files_nii) && ~MRSCont.flags.didCoreg && MRSCont.flags.didLoadData
        gui.b_coreg.Enable = 'on';
    end
    % Segment button
    gui.b_segm = uicontrol('Parent', gui.p2,'Style','PushButton','String','Segment','Enable','off','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.b_segm,'Units','Normalized','Position',[0.1 0.51 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
        set(gui.b_segm,'Callback',{@onSeg});
    if ~isempty(MRSCont.files_nii) && ~MRSCont.flags.didSeg && MRSCont.flags.didCoreg
        gui.b_segm.Enable = 'on';
    end
    % Quantify button
    gui.b_quant = uicontrol('Parent', gui.p2,'Style','PushButton','String','Quantify','Enable','on','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.b_quant,'Units','Normalized','Position',[0.1 0.43 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.b_quant,'Callback',{@onQuant});
    if MRSCont.flags.didQuantify
        gui.b_quant.Enable = 'off';
    else if ~MRSCont.flags.didFit
            gui.b_quant.Enable = 'off';
        end
    end
    % DeIdentify button
    gui.b_deid = uicontrol('Parent', gui.p2,'Style','PushButton','String','DeIdentify','Enable','off','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.b_deid,'Units','Normalized','Position',[0.1 0.2 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.b_deid,'Callback',{@gui_deid});
    % Exit button
    gui.b_exit = uicontrol('Parent', gui.p2,'Style','PushButton','String','Exit','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.b_exit,'Units','Normalized','Position',[0.1 0 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.b_exit,'Callback',{@onExit});

    % Create list of files for the Listbox
    gui.controlPanel = uix.Panel('Parent', gui.leftMenu, 'Title', 'Select a File:','BackgroundColor',gui.colormap.Background);
    set(gui.controlPanel,'Units','Normalized','Position',[0.5 0 0.66 0.1], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold', 'ForegroundColor',gui.colormap.Foreground, 'HighlightColor',gui.colormap.Foreground, 'ShadowColor',gui.colormap.Foreground);
    gui.fileList = MRSCont.files;
    SepFileList = cell(1,length(MRSCont.files));
    gui.RedFileList = cell(1,length(MRSCont.files));
    gui.OnlyFileList = cell(1,length(MRSCont.files));
    for i = 1 : length(MRSCont.files)
        SepFileList{i} =  split(gui.fileList(i), filesep);
        if length(SepFileList{i}) == 1
            SepFileList{i} =  split(gui.fileList(i), '\');
        end
        gui.RedFileList{i} = [filesep SepFileList{i}{end-2} filesep SepFileList{i}{end-1} filesep SepFileList{i}{end}];
        gui.OnlyFileList{i} = [SepFileList{i}{end}];
    end
    clear SepFileList
    gui.ListBox = uicontrol('Style', 'list','BackgroundColor', 'w','FontName', 'Arial','BackgroundColor',gui.colormap.Background, ...
                            'Parent', gui.controlPanel, 'String',gui.RedFileList(:) , ...
                            'Value', gui.SelectedDataset, 'Interruptible', 'on', 'BusyAction', 'cancel', ...
                            'ForegroundColor',gui.colormap.Foreground,...
                            'Callback', {@onListSelection},'KeyPressFcn',@WindowKeyDown, 'KeyReleaseFcn', @WindowKeyUp);

%% Create the display panel tab row
    gui.tabs = uix.TabPanel('Parent', gui.mainLayout, 'Padding', 5, 'FontName', 'Arial',...
                            'FontSize', 16,'SelectionChangedFcn',@(src,event)SelectionChangedFcn(src,event),'BackgroundColor',gui.colormap.Background,...
                            'ForegroundColor', gui.colormap.Foreground, 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
        gui.rawTab      = uix.TabPanel('Parent', gui.tabs, 'Padding', 5,'Padding', 5, 'BackgroundColor',gui.colormap.Background,...
                                        'ForegroundColor', gui.colormap.Foreground, 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground,...
                                        'FontName', 'Arial', 'TabLocation',...
                                       'bottom','FontSize', 10,'SelectionChangedFcn',@(src,event)RawTabChangeFcn(src,event));
        gui.proTab      = uix.TabPanel('Parent', gui.tabs, 'Padding', 5,'Padding', 5, 'BackgroundColor',gui.colormap.Background,...
                                        'ForegroundColor', gui.colormap.Foreground, 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground,...
                                       'FontName', 'Arial', 'TabLocation','bottom',...
                                       'FontSize', 10);
        gui.coregTab    = uix.VBox('Parent', gui.tabs, 'Spacing', 5, 'BackgroundColor',gui.colormap.Background);
        gui.fitTab      = uix.TabPanel('Parent', gui.tabs, 'Padding', 5,'Padding', 5, 'BackgroundColor',gui.colormap.Background,...
                                        'ForegroundColor', gui.colormap.Foreground, 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground,...
                                       'FontName', 'Arial', 'TabLocation','bottom',...
                                       'FontSize', 10,'SelectionChangedFcn',@(src,event)FitTabChangeFcn(src,event));
        gui.quantifyTab = uix.VBox('Parent', gui.tabs, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
        gui.overviewTab = uix.TabPanel('Parent', gui.tabs, 'Padding', 5,'Padding', 5, 'BackgroundColor',gui.colormap.Background,...
                                        'ForegroundColor', gui.colormap.Foreground, 'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground,...
                                       'FontName', 'Arial', 'TabLocation','bottom',...
                                       'FontSize', 10, 'SelectionChangedFcn',@(src,event)OverviewTabChangedFcn(src,event));
    gui.tabs.TabTitles  = {'Raw', 'Processed', 'Cor/Reg', 'Fit', 'Quantified','Overview'};
    gui.tabs.TabWidth   = 115;
    gui.tabs.Selection  = 1;
    gui.tabs.TabEnables = {'off', 'off', 'off', 'off', 'off', 'off'};
    set( gui.mainLayout, 'Widths', [-0.2  -0.8], 'Spacing', 5 );
%% Here we create the inital setup of the tabs
    % Now enable the display tabs depending on which processing steps have
    % been completed:
if MRSCont.flags.didLoadData % Was data loaded at all that can be looked at?
    osp_iniLoadWindow();
end
if MRSCont.flags.didProcess % Has data fitting been run?
    osp_iniProcessWindow();
end
if MRSCont.flags.didCoreg % Have coreg/segment masks been created?
    osp_iniCoregWindow();
end
if MRSCont.flags.didFit % Has data fitting been run?
    osp_iniFitWindow();
end
if MRSCont.flags.didQuantify % Has data fitting been run?
    osp_iniQuantifyWindow();
end
if MRSCont.flags.didOverview % Has data fitting been run?
    osp_iniOverviewWindow();
end
gui.tabs.Selection  = 1;
if ~MRSCont.flags.didLoadData    
    gui.ListBox.Enable = 'off';
end
set(gui.window, 'Visible','on');
%% FUNCTIONS FOR THE LEFT MENU
function updateListBox()
    gui.SelectedDataset = get(gui.ListBox, 'value');
    switch gui.tabs.Selection
        case 1
            osp_updateLoadWindow();
        case 2
            gui.proInfoText = gui.(gui.proTabhandles{gui.SelectedSubFile}).Children(2).Children;
            % Grid for Plot and Data control sliders
            gui.proPlot = gui.(gui.proTabhandles{gui.SelectedSubFile});
            osp_updateProWindow();
        case 3
            osp_updateCoregWindow();
        case 4
            osp_updateFitWindow();
        case 5
            osp_updateQuantifyWindow();
    end
end
%% FUNCTIONS FOR THE ALL TAB
function osp_processingWindow()
    switch gui.tabs.Selection
        case 1
            gui.tabs.TabEnables{1} = 'on';
            gui.dummy = uix.VBox('Parent', gui.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.dummyTxt = uicontrol('Parent',gui.dummy,'style','text', ...
                                          'FontSize', 16, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
            gui.rawTab.TabTitles  = {'...'};
        case 2 
            gui.tabs.TabEnables{2} = 'on';
            gui.dummy = uix.VBox('Parent', gui.proTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.dummyTxt = uicontrol('Parent',gui.dummy,'style','text', ...
                                          'FontSize', 16, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
            gui.proTab.TabTitles  = {'...'};
        case 3 
            gui.tabs.TabEnables{3} = 'on';
            gui.dummy = uix.VBox('Parent', gui.coregTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.dummyTxt = uicontrol('Parent',gui.dummy,'style','text', ...
                                          'FontSize', 16, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
        case 4 
            gui.tabs.TabEnables{4} = 'on';
            gui.dummy = uix.VBox('Parent', gui.fitTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.dummyTxt = uicontrol('Parent',gui.dummy,'style','text', ...
                                          'FontSize', 16, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
            gui.fitTab.TabTitles  = {'...'};
        case 5 
            gui.tabs.TabEnables{5} = 'on';
            gui.dummy = uix.VBox('Parent', gui.quantifyTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.dummyTxt = uicontrol('Parent',gui.dummy,'style','text', ...
                                          'FontSize', 16, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
    end
    set(gui.dummyTxt, 'String', sprintf('In progress...'));        
end
%% FUNCTIONS FOR THE LOAD TAB
function osp_iniLoadWindow()
        if MRSCont.flags.hasRef %Get variables regarding Subspectra
            gui.NoSpecs = gui.NoSpecs + 1;
            gui.SpecNames{2} = 'reference';
        if MRSCont.flags.hasWater
            gui.NoSpecs = gui.NoSpecs + 1;
            gui.SpecNames{3} = 'water';
        end
        else if MRSCont.flags.hasWater
            gui.NoSpecs = gui.NoSpecs + 1;
            gui.SpecNames{2} = 'water';
            end
        end
        gui.tabs.TabEnables{1} = 'on';
        gui.tabs.Selection  = 1;
        gui.EmptydataPlot = 0;
 %%%%%%%%%%%%%%%%%%CREATING SUB TABS FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
 % In this case one tab fo each subspec (A,B,C,D,ref,water)
            gui.metabLoTab = uix.VBox('Parent', gui.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.rawTab.TabWidth   = 115;
            gui.rawTab.Selection  = 1;
            gui.rawTabhandles = {'metabLoTab'};
            if gui.NoSpecs == 2
                if MRSCont.flags.hasRef
                    gui.refLoTab = uix.VBox('Parent', gui.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.rawTab.TabTitles  = gui.SpecNames;
                    gui.rawTab.TabEnables = {'on', 'on'};
                    gui.rawTabhandles = {'metabLoTab', 'refLoTab'};
                end
                if MRSCont.flags.hasWater
                    gui.wLoTab = uix.VBox('Parent', gui.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.rawTab.TabTitles  = gui.SpecNames;
                    gui.rawTab.TabEnables = {'on', 'on'};
                    gui.rawTabhandles = {'metabLoTab', 'wLoTab'};
                end
            end
            if gui.NoSpecs == 3
                gui.refLoTab = uix.VBox('Parent', gui.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                gui.wLoTab = uix.VBox('Parent', gui.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                gui.rawTab.TabTitles  = gui.SpecNames;
                gui.rawTab.TabEnables = {'on', 'on','on'};
                gui.rawTabhandles = {'metabLoTab', 'refLoTab', 'wLoTab'};
            end
 %%%%%%%%%%%%%%%%%%FILLING INFO PANEL FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
 % All the information from the Raw data is read out here
            for t = gui.NoSpecs : -1 : 1
            % Parameter shown in the info panel on top
            gui.dataInfo = uix.Panel('Parent', gui.(gui.rawTabhandles{t}), ...
                                     'Padding', 5, 'Title', ['Actual file: ' MRSCont.files{gui.SelectedDataset}],...
                                     'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                     'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            % Grid for Plot and Data control sliders
            gui.dataPlot = uix.VBox('Parent', gui.(gui.rawTabhandles{t}), 'Padding', 5, 'BackgroundColor',gui.colormap.Background, 'Units', 'normalized');
            gui.dataInfoText  = uicontrol('Parent',gui.dataInfo,'style','text',...
                                          'FontSize', 12, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
            set(gui.(gui.rawTabhandles{t}), 'Heights', [-0.1 -0.9]);

        % Get parameter from file to fill the info panel
        if gui.SelectedSubFile == 1
            StatText = ['Metabolite Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
        else if gui.SelectedSubFile == 2
        StatText = ['Reference Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
            else
                StatText = ['Water Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
            end
        end
        set(gui.dataInfoText, 'String', sprintf(StatText));
 %%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
 %osp_plotLoad is used to visualize the raw data
        temp = figure( 'Visible', 'off' );
        if t == 1
            temp = osp_plotLoad(MRSCont, gui.SelectedDataset,'mets',1 );
            if MRSCont.flags.isUnEdited
                gui.ViewAxes = gca();
                set( gui.ViewAxes, 'Parent', gui.dataPlot );
            end
            if MRSCont.flags.isMEGA
                set( temp.Children(2), 'Parent', gui.dataPlot );
                set( temp.Children(1), 'Parent', gui.dataPlot );
                set(gui.dataPlot,'Heights', [-0.49 -0.49]);
                set(gui.dataPlot.Children(2), 'Units', 'normalized')
                set(gui.dataPlot.Children(2), 'OuterPosition', [0,0.5,1,0.5])
                set(gui.dataPlot.Children(1), 'Units', 'normalized')
                set(gui.dataPlot.Children(1), 'OuterPosition', [0,0,1,0.5])
            end
        else if t == 2
                temp = osp_plotLoad(MRSCont, gui.SelectedDataset,'ref',1 );
                if MRSCont.flags.isUnEdited
                    gui.ViewAxes = gca();
                    set( gui.ViewAxes, 'Parent', gui.dataPlot );
                end
                if MRSCont.flags.isMEGA
                    set( temp.Children(2), 'Parent', gui.dataPlot );
                    set( temp.Children(1), 'Parent', gui.dataPlot );
                    set(gui.dataPlot,'Heights', [-0.49 -0.49]);
                    set(gui.dataPlot.Children(2), 'Units', 'normalized')
                    set(gui.dataPlot.Children(2), 'OuterPosition', [0,0.5,1,0.5])
                    set(gui.dataPlot.Children(1), 'Units', 'normalized')
                    set(gui.dataPlot.Children(1), 'OuterPosition', [0,0,1,0.5])
                end
            else
                temp = osp_plotLoad(MRSCont, gui.SelectedDataset,'w',1 );
                gui.ViewAxes = gca();
                set( gui.ViewAxes, 'Parent', gui.dataPlot );
            end
        end

        % Get rid of the Load figure
        close( temp );
        end
end
function osp_updateLoadWindow()
        gui.dataInfo = gui.(gui.rawTabhandles{gui.SelectedSubFile}).Children(2);
        gui.dataInfoText = gui.(gui.rawTabhandles{gui.SelectedSubFile}).Children(2).Children;
        % Grid for Plot and Data control sliders
        gui.dataPlot = gui.(gui.rawTabhandles{gui.SelectedSubFile});
        gui.EmptydataPlot = 0;
%%%%%%%%%%%%%%%%%%FILLING INFO PANEL FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% All the information from the Raw data is read out here
        if gui.SelectedSubFile == 1
            StatText = ['Metabolite Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
        else if gui.SelectedSubFile == 2
            StatText = ['Reference Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
            else
                StatText = ['Water Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
            end
        end
        set(gui.dataInfoText, 'String',sprintf(StatText))
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
        temp = figure( 'Visible', 'off' );
        if gui.SelectedSubFile == 1
            temp = osp_plotLoad(MRSCont, gui.SelectedDataset,'mets',1 );
            if MRSCont.flags.isUnEdited
                gui.ViewAxes = gca();
                delete(gui.dataPlot.Children(1).Children(1).Children)
                set( gui.ViewAxes.Children, 'Parent', gui.dataPlot.Children(1).Children(1));
                set(  gui.dataPlot.Children(1).Children(1).Title, 'String', gui.ViewAxes.Title.String)
                set(  gui.dataPlot.Children(1).Children(1), 'XLim', gui.ViewAxes.XLim)
            end
            if MRSCont.flags.isMEGA               
                delete(gui.dataPlot.Children(1).Children(1).Children)
                delete(gui.dataPlot.Children(1).Children(2).Children)
                set( temp.Children(2).Children, 'Parent', gui.dataPlot.Children(1).Children(2));
                set( temp.Children(1).Children, 'Parent', gui.dataPlot.Children(1).Children(1));
                set(  gui.dataPlot.Children(1).Children(2).Title, 'String', temp.Children(2).Title.String)
                set(  gui.dataPlot.Children(1).Children(2), 'XLim', temp.Children(2).XLim)
                set(  gui.dataPlot.Children(1).Children(1), 'XLim', temp.Children(1).XLim)
            end
        else if gui.SelectedSubFile == 2
                temp = osp_plotLoad(MRSCont, gui.SelectedDataset,'ref',1 );
                if MRSCont.flags.isUnEdited
                    gui.ViewAxes = gca();
                    delete(gui.dataPlot.Children(1).Children(1).Children)
                    set( gui.ViewAxes.Children, 'Parent', gui.dataPlot.Children(1).Children(1));
                    set(  gui.dataPlot.Children(1).Children(1).Title, 'String', gui.ViewAxes.Title.String)
                    set(  gui.dataPlot.Children(1).Children(1), 'XLim', gui.ViewAxes.XLim)
                end
                if MRSCont.flags.isMEGA               
                    delete(gui.dataPlot.Children(1).Children(1).Children)
                    delete(gui.dataPlot.Children(1).Children(2).Children)
                    set( temp.Children(2).Children, 'Parent', gui.dataPlot.Children(1).Children(2));
                    set( temp.Children(1).Children, 'Parent', gui.dataPlot.Children(1).Children(1));
                    set(  gui.dataPlot.Children(1).Children(2).Title, 'String', temp.Children(2).Title.String)
                    set(  gui.dataPlot.Children(1).Children(2), 'XLim', temp.Children(2).XLim)
                    set(  gui.dataPlot.Children(1).Children(1), 'XLim', temp.Children(1).XLim)
                end
            else
                temp = osp_plotLoad(MRSCont, gui.SelectedDataset,'w',1 );
                gui.ViewAxes = gca();
                delete(gui.dataPlot.Children(1).Children(1).Children)
                set( gui.ViewAxes.Children, 'Parent', gui.dataPlot.Children(1).Children(1));
                set(  gui.dataPlot.Children(1).Children(1).Title, 'String', gui.ViewAxes.Title.String)
                set(  gui.dataPlot.Children(1).Children(1), 'XLim', gui.ViewAxes.XLim)
            end
        end
        % Get rid of the Load figure
        close( temp );    
        set(gui.dataInfo,'Title', ['Actual file: ' MRSCont.files{gui.SelectedDataset}] )
end
%% FUNCTIONS FOR THE PROCESS TAB
function osp_iniProcessWindow()
    gui.tabs.TabEnables{2} = 'on';
        gui.tabs.Selection  = 2;
        gui.EmptyProPlot = 0;
%%%%%%%%%%%%%%%%%%CREATING SUB TABS FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% In this case one tab fo each subspec (A,B,C,D,ref,water)
        gui.AProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
        gui.proTab.TabWidth   = 115;
        gui.proTabhandles = {'AProTab'};
% Set up tabs with regard to the sequence type
        if MRSCont.flags.isUnEdited
            if (MRSCont.flags.hasRef && MRSCont.flags.hasWater)
                gui.refProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.wProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.proTab.TabTitles  = {'metabolites','reference','water'};
                gui.proTab.TabEnables = {'on', 'on','on'};
                gui.proTabhandles = {'AProTab', 'refProTab', 'wProTab'};
            else
                if MRSCont.flags.hasRef
                    gui.refProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.proTab.TabTitles  = {'metabolites','reference'};
                    gui.proTab.TabEnables = {'on', 'on'};
                    gui.proTabhandles = {'AProTab', 'refProTab'};
                end
                if MRSCont.flags.hasWater
                    gui.wProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.proTab.TabTitles  = {'metabolites','water'};
                    gui.proTab.TabEnables = {'on', 'on'};
                    gui.proTabhandles = {'AProTab', 'wProTab'};
                end
            end
        end
        if MRSCont.flags.isMEGA
            if (MRSCont.flags.hasRef && MRSCont.flags.hasWater)
                gui.BProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.diff1ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.sumProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.refProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.wProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.proTab.TabTitles  = {'on','off','diff','sum','reference','water'};
                gui.proTab.TabEnables = {'on', 'on','on', 'on', 'on', 'on'};
                gui.proTabhandles = {'AProTab','BProTab', 'diff1ProTab','sumProTab', 'refProTab', 'wProTab'};
            else
                if MRSCont.flags.hasRef
                    gui.BProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.diff1ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.sumProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.refProTab = uix.VBox('Parent', gui.rawTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.proTab.TabTitles  = {'on','off','diff','sum','reference'};
                    gui.proTab.TabEnables = {'on', 'on','on', 'on', 'on'};
                    gui.proTabhandles = {'AProTab','BProTab', 'diff1ProTab','sumProTab', 'refProTab'};
                end
                if MRSCont.flags.hasWater
                    gui.BProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.diff1ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.sumProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.wProTab = uix.VBox('Parent', gui.rawTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.proTab.TabTitles  = {'on','off','diff','sum','water'};
                    gui.proTab.TabEnables = {'on', 'on','on', 'on', 'on'};
                    gui.proTabhandles = {'AProTab','BProTab', 'diff1ProTab','sumProTab', 'wProTab'};
                end
            end
        end
        if MRSCont.flags.isHERMES
            if (MRSCont.flags.hasRef && MRSCont.flags.hasWater)
                gui.BProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.CProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.DProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.diff1ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.diff2ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.sumProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.refProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.wProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.proTab.TabTitles  = {'A','B','C','D','diff1', 'diff2','sum','reference','water'};
                gui.proTab.TabEnables = {'on', 'on','on', 'on', 'on', 'on', 'on', 'on', 'on'};
                gui.proTabhandles = {'AProTab','BProTab','CProTab','DProTab', 'diff1ProTab', 'diff2ProTab', 'sumProTab', 'refProTab', 'wProTab'};
            else
                if MRSCont.flags.hasRef
                    gui.BProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.CProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.DProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.diff1ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.diff2ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.sumProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);                    
                    gui.refProTab = uix.VBox('Parent', gui.rawTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.proTab.TabTitles  = {'A','B','C','D','diff1', 'diff2','sum','reference'};
                    gui.proTab.TabEnables = {'on', 'on','on','on', 'on', 'on', 'on'};
                    gui.proTabhandles = {'AProTab','BProTab','CProTab','DProTab','diff1ProTab', 'diff2ProTab', 'sumProTab','refProTab'};
                end
                if MRSCont.flags.hasWater
                    gui.BProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.CProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.DProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.diff1ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.diff2ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.sumProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);                    
                    gui.wProTab = uix.VBox('Parent', gui.rawTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.proTab.TabTitles  = {'A','B','C','D','diff1', 'diff2','sum','water'};
                    gui.proTab.TabEnables = {'on', 'on','on','on', 'on', 'on', 'on'};
                    gui.proTabhandles = {'AProTab','BProTab','CProTab','DProTab','diff1ProTab', 'diff2ProTab', 'sumProTab','wProTab'};
                end
            end
        end
        if MRSCont.flags.isHERCULES
            if (MRSCont.flags.hasRef && MRSCont.flags.hasWater)
                gui.BProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.CProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.DProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);   
                gui.diff1ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.diff2ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.sumProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);                
                gui.refProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.wProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                gui.proTab.TabTitles  = {'A','B','C','D','reference','water'};
                gui.proTab.TabEnables = {'on', 'on','on', 'on', 'on', 'on'};
                gui.proTabhandles = {'AProTab','BProTab','CProTab','DProTab','diff1ProTab', 'diff2ProTab', 'sumProTab', 'refProTab', 'wProTab'};
            else
                if MRSCont.flags.hasRef
                    gui.BProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.CProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.DProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.diff1ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.diff2ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.sumProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);                    
                    gui.refProTab = uix.VBox('Parent', gui.rawTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.proTab.TabTitles  = {'A','B','C','D','reference'};
                    gui.proTab.TabEnables = {'on', 'on','on','on'};
                    gui.proTabhandles = {'AProTab','BProTab','CProTab','DProTab','diff1ProTab', 'diff2ProTab', 'sumProTab','refProTab'};
                end
                if MRSCont.flags.hasWater
                    gui.BProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.CProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.DProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.diff1ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.diff2ProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.sumProTab = uix.VBox('Parent', gui.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);                    
                    gui.wProTab = uix.VBox('Parent', gui.rawTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.proTab.TabTitles  = {'A','B','C','D','water'};
                    gui.proTab.TabEnables = {'on', 'on','on','on'};
                    gui.proTabhandles = {'AProTab','BProTab','CProTab','DProTab','diff1ProTab', 'diff2ProTab', 'sumProTab','wProTab'};
                end
            end
        end

%%%%%%%%%%%%%%%%%%FILLING INFO PANEL FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% All the information from the Raw data is read out here
        for t = length(gui.proTabhandles) : -1 : 1
            ind=find(ismember(gui.proTabhandles,[gui.ProNames{t} 'ProTab']));
            gui.proTab.Selection  = ind;
            gui.proInfo = uix.Panel('Parent', gui.(gui.proTabhandles{ind}), ...
                'Padding', 5, 'Title', ['Actual file: ' MRSCont.files{gui.SelectedDataset}],...
                'HighlightColor', gui.colormap.Foreground,'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            % Creates layout for plotting and data control
            gui.proPlot = uix.HBox('Parent', gui.(gui.proTabhandles{ind}), ...
                'Padding', 5,'BackgroundColor', gui.colormap.Background);
            set(gui.(gui.proTabhandles{ind}), 'Heights', [-0.1 -0.9]);
            % Get parameter from file to fill the info panel
            if (strcmp(gui.ProNames{t},'A') || strcmp(gui.ProNames{t},'B') || strcmp(gui.ProNames{t},'C') || strcmp(gui.ProNames{t},'D') || strcmp(gui.ProNames{t},'diff1') || strcmp(gui.ProNames{t},'diff2') || strcmp(gui.ProNames{t},'sum'))
                StatText = ['Metabolite Data -> SNR: ' num2str(MRSCont.QM.SNR.(gui.ProNames{t})(gui.SelectedDataset)) '; FWHM: '...
                            num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)*MRSCont.processed.(gui.ProNames{t}){gui.SelectedDataset}.txfrq/1e6))...
                            ' ppm / Hz \nReference shift: ' num2str(MRSCont.QM.freqShift.(gui.ProNames{t})(gui.SelectedDataset)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.(gui.ProNames{t})(gui.SelectedDataset)*MRSCont.processed.(gui.ProNames{t}){gui.SelectedDataset}.txfrq/1e6)...
                            ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.(gui.ProNames{t})(gui.SelectedDataset)*MRSCont.processed.(gui.ProNames{t}){gui.SelectedDataset}.txfrq/1e6) ' Hz'];
            else if strcmp(gui.ProNames{t},'ref')
            StatText = ['Reference Data -> SNR: ' num2str(MRSCont.QM.SNR.(gui.ProNames{t})(gui.SelectedDataset)) '; FWHM: '...
                        num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)*MRSCont.processed.(gui.ProNames{t}){gui.SelectedDataset}.txfrq/1e6))...
                        ' ppm / Hz'];
                else
                    StatText = ['Water Data -> SNR: ' num2str(MRSCont.QM.SNR.(gui.ProNames{t})(gui.SelectedDataset)) '; FWHM: '...
                                num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)) '/' (num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)*MRSCont.processed.(gui.ProNames{t}){gui.SelectedDataset}.txfrq/1e6))...
                                ' ppm / Hz'];
                end
            end
            gui.proInfoText  = uicontrol('Parent',gui.proInfo,'style','text',...
                                         'FontSize', 12, 'FontName', 'Arial',...
                                         'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                         'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

 %%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
 %osp_plotProcess is used to visualize the processed spectra
            temp = figure( 'Visible', 'off' );
            temp = osp_plotProcess(MRSCont, gui.SelectedDataset,gui.ProNames{t},1 );
                gui.proSpecs = uix.VBox('Parent', gui.proPlot, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.proPre = uix.VBox('Parent', gui.proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    gui.proPost = uix.VBox('Parent', gui.proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                gui.proOut = uix.VBox('Parent', gui.proPlot,'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.proDrift = uix.VBox('Parent', gui.proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    gui.proAlgn = uix.VBox('Parent', gui.proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);

            set( temp.Children(1), 'Parent', gui.proDrift );
            set( temp.Children(1), 'Parent', gui.proAlgn );
            set( temp.Children(1), 'Parent', gui.proPost );
            set( temp.Children(1), 'Parent', gui.proPre );
            close( temp );
%%%%%%%%%%%%%%%%%%DATA CONTROLS FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(gui.proPlot,'Widths', [-0.49 -0.49]);

            set(gui.proPre.Children(1), 'Units', 'normalized')
            set(gui.proPre.Children(1), 'OuterPosition', [0,0,1,1])
            set(gui.proPost.Children(1), 'Units', 'normalized')
            set(gui.proPost.Children(1), 'OuterPosition', [0,0,1,1])
            set(gui.proDrift.Children(1), 'Units', 'normalized')
            set(gui.proDrift.Children(1), 'OuterPosition', [0,0,1,1])
            set(gui.proAlgn.Children(1), 'Units', 'normalized')
            set(gui.proAlgn.Children(1), 'OuterPosition', [0,0,1,1])
        end
    set(gui.proTab,'SelectionChangedFcn',@(src,event)ProTabChangeFcn(src,event));
end
function osp_updateProWindow()
        gui.EmptyProPlot = 0;
        if (MRSCont.flags.isMEGA || MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
            if (gui.SelectedSubSpec == 1 || gui.SelectedSubSpec == 2 || gui.SelectedSubSpec == 3 || gui.SelectedSubSpec == 6)
                t = gui.SelectedSubSpec;
            else if gui.SelectedSubSpec == 4
                    t = 5;
                else
                    t = 4;
                end
            end
        else
            t = gui.SelectedSubSpec;
        end
%%%%%%%%%%%%%%%%%%FILLING INFO PANEL FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% All the information from the Raw data is read out here
        if (strcmp(gui.ProNames{t},'A') || strcmp(gui.ProNames{t},'B') || strcmp(gui.ProNames{t},'C') || strcmp(gui.ProNames{t},'D') || strcmp(gui.ProNames{t},'diff1') || strcmp(gui.ProNames{t},'diff2') || strcmp(gui.ProNames{t},'sum'))
                               StatText = ['Metabolite Data -> SNR: ' num2str(MRSCont.QM.SNR.(gui.ProNames{t})(gui.SelectedDataset)) '; FWHM: '...
                            num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)*MRSCont.processed.(gui.ProNames{t}){gui.SelectedDataset}.txfrq/1e6))...
                            ' ppm / Hz \nReference shift: ' num2str(MRSCont.QM.freqShift.(gui.ProNames{t})(gui.SelectedDataset)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.(gui.ProNames{t})(gui.SelectedDataset)*MRSCont.processed.(gui.ProNames{t}){gui.SelectedDataset}.txfrq/1e6)...
                            ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.(gui.ProNames{t})(gui.SelectedDataset)*MRSCont.processed.(gui.ProNames{t}){gui.SelectedDataset}.txfrq/1e6) ' Hz'];
        else if strcmp(gui.ProNames{t},'ref')
        StatText = ['Reference Data -> SNR: ' num2str(MRSCont.QM.SNR.(gui.ProNames{t})(gui.SelectedDataset)) '; FWHM: '...
                    num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)*MRSCont.processed.(gui.ProNames{t}){t}.txfrq/1e6))...
                    ' ppm / Hz'];
            else
                StatText = ['Water Data -> SNR: ' num2str(MRSCont.QM.SNR.(gui.ProNames{gui.SelectedSubSpec})(gui.SelectedDataset)) '; FWHM: '...
                            num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)) '/' (num2str(MRSCont.QM.FWHM.(gui.ProNames{t})(gui.SelectedDataset)*MRSCont.processed.(gui.ProNames{t}){t}.txfrq/1e6))...
                            ' ppm / Hz'];
            end
        end
        set(gui.proInfoText, 'String',sprintf(StatText))
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
        temp = figure( 'Visible', 'off' );
        temp = osp_plotProcess(MRSCont, gui.SelectedDataset,gui.ProNames{t},1 );
        delete(gui.proDrift.Children.Children)
        delete(gui.proAlgn.Children.Children)
        delete(gui.proPost.Children.Children)
        delete(gui.proPre.Children.Children)
        set( temp.Children(1).Children, 'Parent', gui.proDrift.Children ); % Update drift plot
        set(  gui.proDrift.Children, 'YLim', temp.Children(1).YLim);
        set( temp.Children(2).Children, 'Parent', gui.proAlgn.Children ); % Update aligned and averaged plot
        set(  gui.proAlgn.Children, 'XLim', temp.Children(2).XLim);
        set(  gui.proAlgn.Children, 'YLim', temp.Children(2).YLim);
        set( temp.Children(3).Children, 'Parent', gui.proPost.Children ); % Update post alignment plot
        set(  gui.proPost.Children, 'XLim', temp.Children(3).XLim);
        set(  gui.proPost.Children, 'YLim', temp.Children(3).YLim);
        set( temp.Children(4).Children, 'Parent', gui.proPre.Children ); % Update pre alignment plot
        set(  gui.proPre.Children, 'XLim', temp.Children(4).XLim);
        set(  gui.proPre.Children, 'YLim', temp.Children(4).YLim);
        close( temp );
        set(gui.proInfo,'Title', ['Actual file: ' MRSCont.files{gui.SelectedDataset}] )
end
%% FUNCTIONS FOR THE COREG/SEG TAB
function osp_iniCoregWindow()
    addpath(genpath([spmversion filesep]));
    gui.tabs.TabEnables{3} = 'on';
    gui.tabs.Selection  = 3;
    gui.EmptyQuantPlot = 0;
%%%%%%%%%%%%%%%%%%FILLING INFO PANEL FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% All the information from the Raw data is read out here
    gui.coregInfo = uix.Panel('Parent', gui.coregTab, 'Padding', 5, ...
                              'Title', ['Actual file: ' MRSCont.files{gui.SelectedDataset}],...
                              'FontName', 'Arial','HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
    % Creates layout for plotting and data control
    gui.coregPlot = uix.HBox('Parent', gui.coregTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
    set(gui.coregTab, 'Heights', [-0.1 -0.9]);
    % Get parameter from file to fill the info panel

    StatText = ['Metabolite Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                 '; raw subspecs: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.averages)...
                 '; Sz: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                 num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];

   gui.CoregInfoText  = uicontrol('Parent',gui.coregInfo,'style','text',...
        'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% In this case osp_plotCoreg or osp_plotSegment is used to visualize the
% coregistration or the segmentation
    gui.coregResults = uix.VBox('Parent', gui.coregPlot, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
    temp = figure( 'Visible', 'off' );
    if MRSCont.flags.didSeg
        osp_plotCoreg(MRSCont, gui.SelectedDataset, 1);
        gui.ViewAxes = gca();
        set( gui.ViewAxes, 'Parent', gui.coregResults );
        colormap(gui.coregResults.Children,'gray')
        close( temp );
        temp = figure( 'Visible', 'off' );
        osp_plotSegment(MRSCont, gui.SelectedDataset, 1)
        gui.ViewAxes = gca();
        set( gui.ViewAxes, 'Parent', gui.coregResults );
        colormap(gui.coregResults.Children(1),'gray')
        close( temp );
    else
        osp_plotCoreg(MRSCont, gui.SelectedDataset, 1)
        gui.ViewAxes = gca();
        set( gui.ViewAxes, 'Parent', gui.coregResults );
        colormap(gui.coregResults.Children,'gray')
        close( temp );
    end

    rmpath(genpath([spmversion filesep]));
end

function osp_updateCoregWindow()
        addpath(genpath([spmversion filesep]));
        gui.EmptyProPlot = 0;
%%%%%%%%%%%%%%%%%%FILLING INFO PANEL FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% All the information from the Raw data is read out here
        StatText = ['Metabolite Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.sz) ';  dimensions: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
       set(gui.CoregInfoText, 'String',sprintf(StatText))
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
        if MRSCont.flags.didSeg && length(gui.coregResults.Children) == 2
            delete( gui.coregResults.Children(1) );
            delete( gui.coregResults.Children(1) );
            temp = figure( 'Visible', 'off' );
            temp = osp_plotCoreg(MRSCont, gui.SelectedDataset, 1);
            set( temp.Children, 'Parent', gui.coregResults);
            colormap(gui.coregResults.Children(1),'gray')
            close( temp );
            temp = figure( 'Visible', 'off' );
            temp = osp_plotSegment(MRSCont, gui.SelectedDataset, 1);
            set( temp.Children, 'Parent', gui.coregResults );
            colormap(gui.coregResults.Children(1),'gray')
            close( temp );
        else if MRSCont.flags.didSeg && length(gui.coregResults.Children) == 1
            delete( gui.coregResults.Children(1) );
            temp = figure( 'Visible', 'off' );
            temp = osp_plotCoreg(MRSCont, gui.SelectedDataset, 1);
            set( temp.Children, 'Parent', gui.coregResults);
            colormap(gui.coregResults.Children(1),'gray')
            close( temp );
            temp = figure( 'Visible', 'off' );
            temp = osp_plotSegment(MRSCont, gui.SelectedDataset, 1);
            set( temp.Children, 'Parent', gui.coregResults );
            colormap(gui.coregResults.Children(1),'gray')
            close( temp );
            else if length(gui.coregResults.Children) == 1
                    temp = figure( 'Visible', 'off' );
                    temp = osp_plotCoreg(MRSCont, gui.SelectedDataset, 1);
                    delete( gui.coregResults.Children(1) );
                    set( temp.Children, 'Parent', gui.coregResults );
                    colormap(gui.coregResults.Children,'gray')
                    close( temp );
                else
                    temp = figure( 'Visible', 'off' );
                    temp = osp_plotCoreg(MRSCont, gui.SelectedDataset, 1);
                    set( temp.Children, 'Parent', gui.coregResults );
                    colormap(gui.coregResults.Children,'gray')
                    close( temp ); 
                end
            end
        end

        rmpath(genpath([spmversion filesep]));
end
%% FUNCTIONS FOR THE FIT TAB
function osp_iniFitWindow()
    gui.tabs.TabEnables{4} = 'on';
        gui.tabs.Selection  = 4;
        gui.EmptyFitPlot = 0;
%%%%%%%%%%%%%%%%%%CREATING SUB TABS FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% In this case one tab fo each fit (off,sum,diff1,diff2,ref,water)
         gui.fitTab.TabWidth   = 115;
         for t = 1 : gui.NoFits
                gui.(['fitTab' gui.FitNames{t}]) = uix.VBox('Parent', gui.fitTab, 'Padding', 5,...
                                                            'BackgroundColor',gui.colormap.Background);
                gui.fitTabhandles{t} = ['fitTab' gui.FitNames{t}];
         end
        gui.fitTab.TabTitles  = gui.FitNames;
%%%%%%%%%%%%%%%%%%FILLING INFO PANEL FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% All the information from the Raw data is read out here
        for t = 1 : gui.NoFits
            % Parameter shown in the info panel on top
            gui.fitInfo = uix.Panel('Parent',  gui.(gui.fitTabhandles{t}), ...
                                    'Padding', 5, 'Title', ['Actual file: ' MRSCont.files{gui.SelectedDataset}],...
                                    'FontName', 'Arial','HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,'ForegroundColor',gui.colormap.Foreground);
            % Creates layout for plotting and data control
            gui.fitPlot = uix.HBox('Parent', gui.(gui.fitTabhandles{t}), ...
                                   'Padding', 5,'BackgroundColor',gui.colormap.Background);
            set(gui.(gui.fitTabhandles{t}), 'Heights', [-0.1 -0.9]);
            % Get parameter from file to fill the info panel
            if gui.SelectedSubFile == 1
                StatText = ['Metabolite Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                             '; raw subspecs: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.averages)...
                             '; Sz: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                             num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
            else if gui.SelectedSubFile == 2
            StatText = ['Reference Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                             '; raw subspecs: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                             num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
                else
                    StatText = ['Water Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                             '; raw subspecs: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.averages)...
                             '; Sz: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                             num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
                end
            end
 %%%%%%%%%%%%%%%%%%FILLING FITTED AMPLITUDE PANEL %%%%%%%%%%%%%%%%%%%%%%%%
 % Creates the panel on the right side with the fitted ammplitudes
            gui.fitInfoText  = uicontrol('Parent',gui.fitInfo,'style','text',...
                                        'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                        'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            gui.fitResults = uix.Panel('Parent', gui.fitPlot, 'Padding', 5,...
                                       'Title', ['Raw Amplitudes'],'FontName', 'Arial','HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            RawAmpl = MRSCont.fit.results.(gui.FitNames{gui.SelectedFit}).fitParams{1,gui.SelectedDataset}.ampl .* MRSCont.fit.scale{gui.SelectedDataset};
            if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater)
                if ~(strcmp(gui.FitNames{gui.SelectedFit}, 'ref') || strcmp(gui.FitNames{gui.SelectedFit}, 'w'))
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl)
                        NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.FitNames{gui.SelectedFit}){1,gui.SelectedDataset}.name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                else
                   NameText = ['Water: ' ];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                end
                set(gui.fitResults, 'Title', ['Raw Amplitudes']);
                    gui.fitResultsText = uix.HBox('Parent', gui.fitResults, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.fitResultsTextNames  = uicontrol('Parent',gui.fitResultsText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    gui.fitResultsTextAmpl  = uicontrol('Parent',gui.fitResultsText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
            else
                if ~(strcmp(gui.FitNames{gui.SelectedFit}, 'ref') || strcmp(gui.FitNames{gui.SelectedFit}, 'w'))
                    if MRSCont.flags.hasRef
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results.ref.fitParams{1,gui.SelectedDataset}.ampl .* MRSCont.fit.scale{gui.SelectedDataset});
                    else
                        RawAmpl = RawAmpl ./ (MRSCont.fit.results.water.fitParams{1,gui.SelectedDataset}.ampl .* MRSCont.fit.scale{gui.SelectedDataset});
                    end
                    NameText = [''];
                    RawAmplText = [''];
                    for m = 1 : length(RawAmpl)
                        NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.FitNames{gui.SelectedFit}){1,gui.SelectedDataset}.name{m} ': \n']];
                        RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                    end
                    set(gui.fitResults, 'Title', ['Raw Water Ratio']);
                    gui.fitResultsText = uix.HBox('Parent', gui.fitResults, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                    gui.fitResultsTextNames  = uicontrol('Parent',gui.fitResultsText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                    gui.fitResultsTextAmpl  = uicontrol('Parent',gui.fitResultsText,'style','text',...
                    'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                else
                   NameText = ['Water: ' ];
                   RawAmplText = [num2str(RawAmpl,'%1.2e')];
                   set(gui.fitResults, 'Title', ['Raw Amplitudes']);
                   gui.fitResultsText = uix.HBox('Parent', gui.fitResults, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
                   gui.fitResultsTextNames  = uicontrol('Parent',gui.fitResultsText,'style','text',...
                   'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(NameText),...
                   'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                   gui.fitResultsTextAmpl  = uicontrol('Parent',gui.fitResultsText,'style','text',...
                   'FontSize', 11, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(RawAmplText),...
                   'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                end
            end
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
%osp_plotFit is used to visualize the fits (off,diff1,diff2,sum,ref,water)
            temp = figure( 'Visible', 'off' );
            temp = osp_plotFit(MRSCont, gui.SelectedDataset,gui.FitNames{gui.SelectedFit},1 );
            gui.ViewAxes = gca();
            set( gui.ViewAxes, 'Parent', gui.fitPlot );
            close( temp );
%%%%%%%%%%%%%%%%%%DATA CONTROLS FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creates the scroll bar for the datasets of the MRSCont on the right side
            set(gui.fitPlot,'Widths', [-0.16 -0.84]);
        end
end

function osp_updateFitWindow()
        gui.fitPlot = gui.(gui.fitTabhandles{gui.SelectedFit}).Children(1).Children(2);
        gui.fitInfo = gui.(gui.fitTabhandles{gui.SelectedFit}).Children(2);
        gui.fitInfoText = gui.(gui.fitTabhandles{gui.SelectedFit}).Children(2).Children;
        gui.fitResults  = gui.(gui.fitTabhandles{gui.SelectedFit}).Children(1).Children(1);
        gui.fitResultsTextAmpl = gui.(gui.fitTabhandles{gui.SelectedFit}).Children(1).Children(1).Children(1).Children(1);
        gui.fitResultsTextNames = gui.(gui.fitTabhandles{gui.SelectedFit}).Children(1).Children(1).Children(1).Children(2);
        
        gui.EmptyFitPlot = 0;
%%%%%%%%%%%%%%%%%%FILLING INFO PANEL FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% All the information from the Raw data is read out here
        if gui.SelectedFit == 1
            StatText = ['Metabolite Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.sz) ';  dimensions: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
        else if gui.SelectedFit == 2
        StatText = ['Reference Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
            else
                StatText = ['Water Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
            end
        end
        set(gui.fitInfoText, 'String',sprintf(StatText))
        % Update amplitudes for the fit results panel based on the files in the MRSCont (Raw Amplitudes or Water-scaled if ref or water supplied)
        RawAmpl = MRSCont.fit.results.(gui.FitNames{gui.SelectedFit}).fitParams{1,gui.SelectedDataset}.ampl .* MRSCont.fit.scale{gui.SelectedDataset};
        if ~(MRSCont.flags.hasRef || MRSCont.flags.hasWater)
            if ~(strcmp(gui.FitNames{gui.SelectedFit}, 'ref') || strcmp(gui.FitNames{gui.SelectedFit}, 'w'))
                NameText = [''];
                RawAmplText = [''];
                for m = 1 : length(RawAmpl)
                    NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.FitNames{gui.SelectedFit}){1,gui.SelectedDataset}.name{m} ': \n']];
                    RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                end
            else
               NameText = ['Water: ' ];
               RawAmplText = [num2str(RawAmpl,'%1.2e')];
            end
            set(gui.fitResults, 'Title', ['Raw Amplitudes']);
            set(gui.fitResultsTextNames, 'String',sprintf(NameText));
            set(gui.fitResultsTextAmpl, 'String',sprintf(RawAmplText));
        else
            if ~(strcmp(gui.FitNames{gui.SelectedFit}, 'ref') || strcmp(gui.FitNames{gui.SelectedFit}, 'w'))
                if MRSCont.flags.hasRef
                    RawAmpl = RawAmpl ./ (MRSCont.fit.results.ref.fitParams{1,gui.SelectedDataset}.ampl .* MRSCont.fit.scale{gui.SelectedDataset});
                else
                    RawAmpl = RawAmpl ./ (MRSCont.fit.results.water.fitParams{1,gui.SelectedDataset}.ampl .* MRSCont.fit.scale{gui.SelectedDataset});
                end
                NameText = [''];
                RawAmplText = [''];
                for m = 1 : length(RawAmpl)
                    NameText = [NameText, [MRSCont.fit.resBasisSet.(gui.FitNames{gui.SelectedFit}){1,gui.SelectedDataset}.name{m} ': \n']];
                    RawAmplText = [RawAmplText, [num2str(RawAmpl(m),'%1.2e') '\n']];
                end
                set(gui.fitResults, 'Title', ['Raw Water Ratio']);
                set(gui.fitResultsTextNames, 'String',sprintf(NameText));
                set(gui.fitResultsTextAmpl, 'String',sprintf(RawAmplText));
            else
               NameText = ['Water: \t'];
               RawAmplText = [num2str(RawAmpl,'%1.2e')];
               set(gui.fitResults, 'Title', ['Raw Amplitudes']);
                set(gui.fitResultsTextNames, 'String',sprintf(NameText));
                set(gui.fitResultsTextAmpl, 'String',sprintf(RawAmplText));
            end
        end
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
        temp = figure( 'Visible', 'off' );
        temp = osp_plotFit(MRSCont, gui.SelectedDataset,gui.FitNames{gui.SelectedFit},1 );
        gui.ViewAxes = gca();
        delete(gui.fitPlot.Children)
        set( gui.ViewAxes.Children, 'Parent', gui.fitPlot); %Update plot
        set(  gui.fitPlot.Title, 'String', gui.ViewAxes.Title.String) %Update title
        set(  gui.fitPlot, 'XLim', gui.ViewAxes.XLim) % Update Xlim
        set(  gui.fitPlot, 'YLim', gui.ViewAxes.YLim) % Update Ylim
        % Get rid of the Load figure
        close( temp );
        set(gui.fitInfo,'Title', ['Actual file: ' MRSCont.files{gui.SelectedDataset}] )
end
%% FUNCTIONS FOR THE QUANTIFY TAB
function osp_iniQuantifyWindow()
gui.tabs.TabEnables{5} = 'on';
        gui.tabs.Selection  = 5;
        gui.EmptyQuantPlot = 0;
%%%%%%%%%%%%%%%%%%FILLING INFO PANEL FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% All the information from the Raw data is read out here
        gui.quantInfo = uix.Panel('Parent', gui.quantifyTab, 'Padding', 5, ...
                                  'Title', ['Actual file: ' MRSCont.files{gui.SelectedDataset}],...
                                  'FontName', 'Arial','HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        % Creates layout for plotting and data control
        gui.quantPlot = uix.HBox('Parent', gui.quantifyTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
        set(gui.quantifyTab, 'Heights', [-0.1 -0.9]);
        % Get parameter from file to fill the info panel
        if gui.SelectedSubFile == 1
            StatText = ['Metabolite Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
        else if gui.SelectedSubFile == 2
        StatText = ['Reference Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
            else
                StatText = ['Water Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
            end
        end
       gui.QuantInfoText  = uicontrol('Parent',gui.quantInfo,'style','text',...
            'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
            'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% In this case a table is created based on a uicontol slider
        QuantText = cell(length(MRSCont.quantify.metabs)+1,gui.NoQuants);
        QuantText{1,1} = 'Metabolite';
        QuantText(2:end,1) = MRSCont.quantify.metabs';
            for q = 1 : gui.NoQuants
                QuantText(1,q+1) = gui.QuantNames(q);
                QuantText(2:end,q+1) = table2cell(MRSCont.quantify.tables.(gui.QuantModelNames{gui.SelectedModel}).(gui.QuantNames{q})(gui.SelectedDataset,:))';
            end
        temp=uimulticollist ( 'units', 'normalized', 'position', [0 0 1 1], 'string', QuantText,...
            'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        set ( temp, 'BackgroundColor',gui.colormap.Background);
        set( temp, 'Parent', gui.quantPlot );
end

function osp_updateQuantifyWindow()
        gui.EmptyQuantPlot = 0;
%%%%%%%%%%%%%%%%%%FILLING INFO PANEL FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
% All the information from the Raw data is read out here
        if gui.SelectedSubFile == 1
            StatText = ['Metabolite Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.sz) ';  dimensions: ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
        else if gui.SelectedSubFile == 2
        StatText = ['Reference Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_ref{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
            else
                StatText = ['Water Data -> Sequence: ' SeqName '; B0: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.te) ' / ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.spectralwidth) ' Hz'...
                         '; raw subspecs: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2})) ' x ' num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})) ' mm = '...
                         num2str(MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{1}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{2}) * MRSCont.raw_w{1,gui.SelectedDataset}.geometry.size.(gui.GeometryNames{3})/1000) ' ml'];
            end
        end
        set(gui.QuantInfoText, 'String',sprintf(StatText))
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
        QuantText = cell(length(MRSCont.quantify.metabs)+1,gui.NoQuants);
        QuantText{1,1} = 'Metabolite';
        QuantText(2:end,1) = MRSCont.quantify.metabs';
        for q = 1 : gui.NoQuants
            QuantText(1,q+1) = gui.QuantNames(q);
            QuantText(2:end,q+1) = table2cell(MRSCont.quantify.tables.(gui.QuantModelNames{gui.SelectedModel}).(gui.QuantNames{q})(gui.SelectedDataset,:))';
        end
        temp=uimulticollist ( 'units', 'normalized', 'position', [0 0 1 1], 'string', QuantText,...
            'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
         set( temp, 'BackgroundColor',gui.colormap.Background);
        delete(gui.quantPlot.Children)
        set( temp, 'Parent', gui.quantPlot ); %Update table
        set(gui.quantInfo,'Title', ['Actual file: ' MRSCont.files{gui.SelectedDataset}]) %Update info Title
end
%% FUNCTIONS FOR THE OVERVIEW TAB
function osp_iniOverviewWindow()
    gui.tabs.TabEnables{6} = 'on';
        gui.tabs.Selection  = 6;
 %%%%%%%%%%%%%%%%%%CREATING SUB TABS FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
        gui.specsOvTab = uix.HBox('Parent', gui.overviewTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
        gui.meanOvTab = uix.HBox('Parent', gui.overviewTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
        gui.quantOvTab = uix.HBox('Parent', gui.overviewTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
        gui.distrOvTab = uix.HBox('Parent', gui.overviewTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
        gui.corrOvTab = uix.HBox('Parent', gui.overviewTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
        gui.diceOvTab = uix.HBox('Parent', gui.overviewTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
        gui.overviewTab.TabTitles  = {'spectra', 'mean spectra', 'quantify table', 'distribution', 'correlation','dice ovelap'};
        gui.overviewTab.TabWidth   = 115;
        gui.overviewTab.Selection  = 1;
        gui.overviewTab.TabEnables = {'on', 'on', 'on', 'on', 'on', 'off'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Overview Panel for all specs sorted by groups
        gui.specsOvPlot = uix.VBox(...
            'Parent', gui.specsOvTab, ...
            'Padding', 5,'BackgroundColor',gui.colormap.Background);
%%%%%%%%%%%%%%%%%%DATA CONTROLS FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creates popup menu for the processed Subspectra (A,B,C,D,ref,water)
        gui.specsOvPlotControls = uix.Panel('Parent', gui.specsOvPlot,'Title', 'Actual Subspectra', ...
                                            'Padding', 5,'HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.pop_specsOvPlotControls = uicontrol('Parent',gui.specsOvPlotControls,'style','popupmenu',...
                                                'Units', 'Normalized', 'Position', [0 0 1 1],'FontName', 'Arial', ...
                                                'String',gui.proTab.TabTitles, 'Value', 1,...
                                                'callback',{@pop_specsOvPlot_Call});
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
%op_plotspec is used to visualize the processed data
        gui.shiftind = 0.2;
        for g = 1 :  gui.NoGroups
            temp = figure( 'Visible', 'off' );
            if (strcmp(gui.ProNames{gui.SelectedSubSpec},'A') || strcmp(gui.ProNames{gui.SelectedSubSpec},'D') || strcmp(gui.ProNames{gui.SelectedSubSpec},'C') || strcmp(gui.ProNames{gui.SelectedSubSpec},'D'))
                if gui.NAAnormed ==1
                    shift = gui.shiftind * (g-1);
                    temp = op_plotspec(MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(gui.ProNames{gui.SelectedSubSpec}),2,1,cb(g,:),shift,['Overview ' gui.proTab.TabTitles{gui.SelectedSubSpec}]);
                else
                    ylimmax = max(real(MRSCont.overview.all_data.(gui.ProNames{gui.SelectedSubFile}){1,1}.specs));
                    shift = ylimmax * gui.shiftind * (g-1);
                    temp = op_plotspec(MRSCont.overview.(['sort_data_g' num2str(g)]).(gui.ProNames{gui.SelectedSubSpec}),2,1,cb(g,:),shift,['Overview ' gui.proTab.TabTitles{gui.SelectedSubSpec}]);
                end
            else
                ylimmax = max(real(MRSCont.overview.all_data.(gui.ProNames{1}){1,1}.specs));
                shift = ylimmax * gui.shiftind * (g-1);
                temp = op_plotspec(MRSCont.overview.(['sort_data_g' num2str(g)]).(gui.ProNames{gui.SelectedSubSpec}),2,1,cb(g,:),shift,['Overview ' gui.proTab.TabTitles{gui.SelectedSubSpec}]);
            end
            set(gca, 'YColor', MRSCont.colormap.Background);
            set(gca,'YTickLabel',{})
            set(gca,'YTick',{});
            set(gca,'XColor',MRSCont.colormap.Foreground);
            set(gca,'Color','w');
            set(gcf,'Color','w');
            title(['Overview ' gui.proTab.TabTitles{gui.SelectedSubFile}],'Color', MRSCont.colormap.Foreground);
            if g == 1
                ax=get(temp,'Parent');
                figpl = get(ax,'Parent');
                gui.ViewAxes = gca();
                set( gui.ViewAxes, 'Parent', gui.specsOvPlot);
                % Get rid of the Load figure
                close( figpl );
            else
                ax=get(temp,'Parent');
                figpl = get(ax,'Parent');
                copyobj(ax.Children, gui.ViewAxes);
                % Get rid of the Load figure
                close( figpl );
            end
        end
        if gui.SelectedSubFile ==1
            set(gui.specsOvPlot.Children(2), 'XLim', [0.2 4.5])
        else
            set(gui.specsOvPlot.Children(2), 'XLim', [0 2*4.68])
        end
        set(gui.specsOvPlot,'Heights', [-0.07 -0.93]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Overview Panel for mean specs sorted by groups
       gui.overviewTab.Selection  = 2;
       gui.meanOvPlot = uix.VBox('Parent', gui.meanOvTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
%%%%%%%%%%%%%%%%%%DATA CONTROLS FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creates popup menu for the processed Subspectra (A,B,C,D,ref,water)
       gui.meanOvPlotControls = uix.Panel('Parent', gui.meanOvPlot,'Title', 'Actual Subspectra', ...
                                          'Padding', 5,'HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
       gui.pop_meanOvPlotControls = uicontrol('Parent',gui.meanOvPlotControls,'style','popupmenu',...
                                              'Units', 'Normalized', 'Position', [0 0 1 1],'FontName', 'Arial', ...
                                              'String',gui.proTab.TabTitles, 'Value', 1,'callback',{@pop_meanOvPlot_Call});
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
%op_plotspec is used for a dummy plot which is update later
        gui.shift = 0.5;
        temp = figure( 'Visible', 'off' );
        if gui.SelectedSubFile ==1
            temp = op_plotspec(MRSCont.overview.(['sort_data_g' num2str(1)]).(gui.ProNames{2}),2,1,cb(1,:),gui.shift*(1-1),['Overview ' gui.proTab.TabTitles{gui.SelectedSubFile}]);
        else
            temp = op_plotspec(MRSCont.overview.(['sort_data_g' num2str(1)]).(gui.ProNames{1}),2,1,cb(1,:),gui.shift*(1-1),['Overview ' gui.proTab.TabTitles{gui.SelectedSubFile}]);
        end
        set(gca, 'YColor', MRSCont.colormap.Background);
        set(gca,'YTickLabel',{})
        set(gca,'YTick',{});
        set(gca,'XColor',MRSCont.colormap.Foreground);
        set(gca,'Color','w');
        set(gcf,'Color','w');
        title(['Overview ' gui.proTab.TabTitles{gui.SelectedSubFile}],'Color', MRSCont.colormap.Foreground);
        ax=get(temp,'Parent');
        figpl = get(ax,'Parent');
        gui.ViewAxes = gca();
        set( gui.ViewAxes, 'Parent', gui.meanOvPlot);
        close( figpl );
        if gui.SelectedSubFile ==1
            set(gui.meanOvPlot.Children(2), 'XLim', [0.2 4.5])
        else
            set(gui.meanOvPlot.Children(2), 'XLim', [0 2*4.68])
        end
        osp_updatemeanOvWindow() %Update the plot with the mean and SD
        set(gui.meanOvPlot,'Heights', [-0.07 -0.93]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Overview Panel for quantifications
        gui.overviewTab.Selection  = 3;
        gui.quantOvPlot = uix.VBox('Parent', gui.quantOvTab,'Padding', 5,'BackgroundColor',gui.colormap.Background);
%%%%%%%%%%%%%%%%%%DATA CONTROLS FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creates Popup menu to change between quantifications (tCr, waterScaled etc.)
        gui.quantOvPlotControls = uix.Panel('Parent', gui.quantOvPlot,'Title', 'Actual Quantification', ...
                                            'Padding', 5,'HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.pop_quantOvPlotControls = uicontrol('Parent',gui.quantOvPlotControls,'style','popupmenu',...
                                                'Units', 'Normalized', 'Position', [0 0 1 1],'FontName', 'Arial', ...
                                                'String',gui.QuantNames, 'Value', 1,...
                                                'callback',{@pop_quantOvPlot_Call});
 %%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
 % Quantification table is created based on uicontrol
        gui.quantOvResults = uix.Panel('Parent', gui.quantOvPlot, 'Padding', 5, ...
                                        'Title', ['Results: ' (gui.QuantNames{gui.SelectedQuant})],...
                                        'FontName', 'Arial','HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        QuantTextOv = cell(MRSCont.nDatasets+1,gui.NoQuantMetabs);
        QuantTextOv(1,:) = MRSCont.quantify.metabs;
        QuantTextOv(2:end,:) = table2cell(MRSCont.quantify.tables.(gui.QuantModelNames{gui.SelectedModel}).(gui.QuantNames{gui.SelectedQuant})(:,:));
        temp=uimulticollist ( 'units', 'normalized', 'position', [0 0 1 1], 'string', QuantTextOv,...
            'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        set(temp,'BackgroundColor',gui.colormap.Background)
        set( temp, 'Parent', gui.quantOvResults );
        set(gui.quantOvPlot,'Heights', [-0.07 -0.93]);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Overview Panel for distributions/raincloud plots
        gui.overviewTab.Selection  = 4;
        gui.distrOvPlot = uix.VBox('Parent', gui.distrOvTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
%%%%%%%%%%%%%%%%%%DATA CONTROLS FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creates popup menus for differnt quantifications and metabolites
        gui.distrOvPanelControls = uix.Panel('Parent', gui.distrOvPlot,'Title', 'Actual Quantification and Metabolite', ...
                                             'Padding', 5,'HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.distrOvControls = uix.HBox('Parent', gui.distrOvPanelControls,...
                                       'Padding', 5, 'Spacing', 10,'BackgroundColor',gui.colormap.Background);
        gui.pop_distrOvQuantControls = uicontrol('Parent',gui.distrOvControls,'style','popupmenu',...
                                                 'Units', 'Normalized', 'Position', [0 0 1 1],'FontName', 'Arial',...
                                                 'String',gui.QuantNames, 'Value', 1,...
                                                 'callback',{@pop_distrOvQuant_Call});
        gui.pop_distrOvMetabControls = uicontrol('Parent',gui.distrOvControls,'style','popupmenu',...
                                                 'Units', 'Normalized', 'Position', [0 0 1 1],'FontName', 'Arial', ...
                                                 'String',MRSCont.quantify.metabs, 'Value', 1,...
                                                 'callback',{@pop_distrOvMetab_Call});
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
%osp_plotQuantifyTable to create distribution overview as raincloud plot
        temp = figure( 'Visible', 'off' );
        [temp] = osp_plotRaincloud(MRSCont,gui.QuantModelNames{gui.SelectedModel}, gui.QuantNames{gui.SelectedQuant},MRSCont.quantify.metabs{gui.SelectedMetab},'Raincloud plot',1);
        gui.ViewAxes = gca();
        set( gui.ViewAxes, 'Parent', gui.distrOvPlot);
        close( temp );
        set(gui.distrOvPlot,'Heights', [-0.07 -0.90 -0.03]);
        gui.distrOvPlot.Children(3).Legend.Location = 'North';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Overview Panel for correlation plots
        gui.overviewTab.Selection  = 5;
        gui.corrOvPlot = uix.VBox('Parent', gui.corrOvTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
 %%%%%%%%%%%%%%%%%%DATA CONTROLS FOR THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates popup menu for differnt quantification, metabolite and
% correaltion measure
        gui.corrOvPanelControls = uix.Panel('Parent', gui.corrOvPlot,'Title', 'Actual Quantification, Metabolite, and Correlation Measure', ...
                                            'Padding', 5,'HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.corrOvControls = uix.HBox('Parent', gui.corrOvPanelControls,'Padding', 5, 'Spacing', 10,'BackgroundColor',gui.colormap.Background);
        gui.pop_corrOvQuantControls = uicontrol('Parent',gui.corrOvControls,'style','popupmenu',...
                                                'Units', 'Normalized', 'Position', [0 0 1 1],'FontName', 'Arial',...
                                                'String',gui.QuantNames, 'Value', 1,'callback',{@pop_corrOvQuant_Call});
        gui.pop_corrOvMetabControls = uicontrol('Parent',gui.corrOvControls,'style','popupmenu',...
                                                'Units', 'Normalized', 'Position', [0 0 1 1],'FontName', 'Arial', ...
                                                'String',MRSCont.quantify.metabs, 'Value', 1,'callback',{@pop_corrOvMetab_Call});
        gui.pop_corrOvCorrControls = uicontrol('Parent',gui.corrOvControls,'style','popupmenu',...
                                               'Units', 'Normalized', 'Position', [0 0 1 1],'FontName', 'Arial', ...
                                               'String',gui.CorrNames, 'Value', 1,'callback',{@pop_corrOvCorr_Call});
        gui.pop_whichcorrOvCorrControls = uicontrol('Parent',gui.corrOvControls,'style','popupmenu',...
                                       'Units', 'Normalized', 'Position', [0 0 1 1],'FontName', 'Arial', ...
                                       'String',{'MRSCont.overview.corr','metabolites','QM'}, 'Value', 1,'callback',{@pop_whichcorrOvCorr_Call});
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
%osp_plotQuantifyTable is used to create a correlation plot
        temp = figure( 'Visible', 'off' );
        [temp] = osp_plotScatter(MRSCont, gui.QuantModelNames{gui.SelectedModel}, gui.QuantNames{gui.SelectedQuant},MRSCont.quantify.metabs{gui.SelectedMetab},gui.CorrMeas{gui.SelectedCorr},gui.CorrNames{gui.SelectedCorr},1);
        gui.ViewAxes = gca();
        set( gui.ViewAxes, 'Parent', gui.corrOvPlot);
        set(gui.corrOvPlot,'Heights', [-0.07 -0.90 -0.03]);
        gui.corrOvPlot.Children(3).Legend.Location = 'North';
        close( temp );
        gui.overviewTab.Selection  = 2;
end
function osp_updateSpecsOvWindow() %updates the overview tab for the spectra
        delete(gui.specsOvPlot.Children(2).Children)
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
        for g = 1 :  gui.NoGroups
            temp = figure( 'Visible', 'off' );
            if (strcmp(gui.ProNames{gui.SelectedSubSpec},'A') || strcmp(gui.ProNames{gui.SelectedSubSpec},'D') || strcmp(gui.ProNames{gui.SelectedSubSpec},'C') || strcmp(gui.ProNames{gui.SelectedSubSpec},'D'))
                if gui.NAAnormed ==1
                    shift = gui.shiftind * (g-1);
                    temp = op_plotspec(MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(gui.ProNames{gui.SelectedSubSpec}),2,1,cb(g,:),shift,['Overview ' gui.proTab.TabTitles{gui.SelectedSubSpec}]);
                else
                    ylimmax = max(real(MRSCont.overview.all_data.(gui.ProNames{gui.SelectedSubFile}){1,1}.specs));
                    shift = ylimmax * gui.shiftind * (g-1);
                    temp = op_plotspec(MRSCont.overview.(['sort_data_g' num2str(g)]).(gui.ProNames{gui.SelectedSubSpec}),2,1,cb(g,:),shift,['Overview ' gui.proTab.TabTitles{gui.SelectedSubSpec}]);
                end
            else
                ylimmax = max(real(MRSCont.overview.all_data.(gui.ProNames{1}){1,1}.specs));
                shift = ylimmax * gui.shiftind * (g-1);
                temp = op_plotspec(MRSCont.overview.(['sort_data_g' num2str(g)]).(gui.ProNames{gui.SelectedSubSpec}),2,1,cb(g,:),shift,['Overview ' gui.proTab.TabTitles{gui.SelectedSubSpec}]);
            end
            set(gca, 'YColor', MRSCont.colormap.Background);
            set(gca,'YTickLabel',{})
            set(gca,'YTick',{});
            set(gca,'XColor',MRSCont.colormap.Foreground);
            set(gca,'Color','w');
            set(gcf,'Color','w');
            title(['Overview ' gui.proTab.TabTitles{gui.SelectedSubFile}],'Color', MRSCont.colormap.Foreground);
                ax=get(temp,'Parent');
                figpl = get(ax,'Parent');
                copyobj(ax.Children, gui.specsOvPlot.Children(2));
                close( figpl );
        end
        if gui.SelectedSubSpec ==1 % Update title and limits
            set(gui.specsOvPlot.Children(2), 'XLim', [0.2 4.5])
            set(gui.specsOvPlot.Children(2).Title, 'String', ['Overview ' gui.proTab.TabTitles{gui.SelectedSubSpec}])
        else
            set(gui.specsOvPlot.Children(2), 'XLim', [0 2*4.68])
            set(gui.specsOvPlot.Children(2).Title, 'String', ['Overview ' gui.proTab.TabTitles{gui.SelectedSubSpec}])
        end
end

function osp_updatemeanOvWindow() %updates the mean and sd overview tab
        delete(gui.meanOvPlot.Children(2).Children)
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
        for g = 1 :  gui.NoGroups
            temp = figure;
            hold on
            if (strcmp(gui.ProNames{gui.SelectedSubSpec},'A') || strcmp(gui.ProNames{gui.SelectedSubSpec},'D') || strcmp(gui.ProNames{gui.SelectedSubSpec},'C') || strcmp(gui.ProNames{gui.SelectedSubSpec},'D'))
                if gui.NAAnormed ==1
                    shift = gui.shift * (g-1);
                    yu = MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.ProNames{gui.SelectedSubSpec}]) + ...
                        MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['sd_' gui.ProNames{gui.SelectedSubSpec}]);
                    yl = MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.ProNames{gui.SelectedSubSpec}]) - ...
                        MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['sd_' gui.ProNames{gui.SelectedSubSpec}]);
                    temp = fill([MRSCont.overview.(['ppm_' (gui.ProNames{gui.SelectedSubSpec})]) fliplr(MRSCont.overview.(['ppm_' (gui.ProNames{gui.SelectedSubSpec})]))], [yu+shift fliplr(yl)+shift], [0 0 0],'FaceAlpha',0.1, 'linestyle', 'none');
                    plot(MRSCont.overview.(['ppm_' (gui.ProNames{gui.SelectedSubSpec})]),MRSCont.overview.(['sort_data_g' num2str(g) '_NAAnormalized']).(['mean_' gui.ProNames{gui.SelectedSubSpec}])+shift ,'color',cb(g,:), 'LineWidth', 1);
                else
                    ylimmax = max(MRSCont.overview.(['sort_data_g' num2str(1)]).(['mean_' gui.ProNames{2}]));
                    shift = ylimmax * gui.shift * (g-1);
                    yu = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.ProNames{gui.SelectedSubSpec}]) + ...
                        MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.ProNames{gui.SelectedSubSpec}]);
                    yl = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.ProNames{gui.SelectedSubSpec}]) - ...
                        MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.ProNames{gui.SelectedSubSpec}]);
                    temp = fill([MRSCont.overview.(['ppm_' (gui.ProNames{gui.SelectedSubSpec})]) fliplr(MRSCont.overview.(['ppm_' (gui.ProNames{gui.SelectedSubSpec})]))], [yu+shift fliplr(yl)+shift], [0 0 0],'FaceAlpha',0.1, 'linestyle', 'none');
                    plot(MRSCont.overview.(['ppm_' (gui.ProNames{gui.SelectedSubSpec})]),MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.ProNames{gui.SelectedSubSpec}])+shift ,'color',cb(g,:), 'LineWidth', 1);
                end

            else
                ylimmax = max(MRSCont.overview.(['sort_data_g' num2str(1)]).(['mean_' gui.ProNames{1}]));
                shift = ylimmax * gui.shift * (g-1);
                yu = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.ProNames{gui.SelectedSubSpec}]) + ...
                    MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.ProNames{gui.SelectedSubSpec}]);
                yl = MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.ProNames{gui.SelectedSubSpec}]) - ...
                    MRSCont.overview.(['sort_data_g' num2str(g)]).(['sd_' gui.ProNames{gui.SelectedSubSpec}]);
                temp = fill([MRSCont.overview.(['ppm_' (gui.ProNames{gui.SelectedSubSpec})]) fliplr(MRSCont.overview.(['ppm_' (gui.ProNames{gui.SelectedSubSpec})]))], [yu+shift fliplr(yl)+shift], [0 0 0],'FaceAlpha',0.1, 'linestyle', 'none');
                plot(MRSCont.overview.(['ppm_' (gui.ProNames{gui.SelectedSubSpec})]),MRSCont.overview.(['sort_data_g' num2str(g)]).(['mean_' gui.ProNames{gui.SelectedSubSpec}])+shift ,'color',cb(g,:), 'LineWidth', 1);
            end
                figpl = get(temp.Parent,'Parent');
                set(temp.Parent.Children,'Parent',gui.meanOvPlot.Children(2))
                close(figpl);
        end
        if (strcmp(gui.ProNames{gui.SelectedSubSpec},'A') || strcmp(gui.ProNames{gui.SelectedSubSpec},'D') || strcmp(gui.ProNames{gui.SelectedSubSpec},'C') || strcmp(gui.ProNames{gui.SelectedSubSpec},'D'))
            set(gui.meanOvPlot.Children(2), 'XLim', [0.2 4.5])
            set(gui.meanOvPlot.Children(2).Title, 'String', ['Mean \pm SD ' gui.proTab.TabTitles{gui.SelectedSubSpec}])
        else
            set(gui.meanOvPlot.Children(2), 'XLim', [0 2*4.68])
            set(gui.meanOvPlot.Children(2).Title, 'String', ['Mean \pm SD ' gui.proTab.TabTitles{gui.SelectedSubSpec}])
        end
end

function osp_updatequantOvWindow() %updates the quantification table overview tab
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
        QuantTextOv = cell(MRSCont.nDatasets+1,gui.NoQuantMetabs);
        QuantTextOv(1,:) = MRSCont.quantify.metabs;
        QuantTextOv(2:end,:) = table2cell(MRSCont.quantify.tables.(gui.QuantModelNames{gui.SelectedModel}).(gui.QuantNames{gui.SelectedQuant})(:,:));
        temp=uimulticollist ( 'units', 'normalized', 'position', [0 0 1 1], 'string', QuantTextOv);
        set( temp, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        set( temp, 'Parent', gui.quantOvResults );
        set(gui.quantOvResults, 'Title', ['Results: ' (gui.QuantNames{gui.SelectedQuant})]);
end

function osp_updatedistrOvWindow() %updates the raincloud plot overview tab
            delete(gui.distrOvPlot.Children(3).Children)
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
            temp = figure( 'Visible', 'off' );
            [temp] = osp_plotRaincloud(MRSCont,gui.QuantModelNames{gui.SelectedModel},gui.QuantNames{gui.SelectedQuant},MRSCont.quantify.metabs{gui.SelectedMetab},'Raincloud plot',1);
            set( temp.Children(2).Children, 'Parent', gui.distrOvPlot.Children(3) );
            set(  gui.distrOvPlot.Children(3), 'XLabel', temp.Children(2).XLabel);
            set(  gui.distrOvPlot.Children(3), 'YLim', temp.Children(2).YLim);
            close(temp);
            set(gui.distrOvPlot,'Heights', [-0.07 -0.90 -0.03]);
            gui.distrOvPlot.Children(3).Legend.Location = 'North'; % Update legend
            set(gui.distrOvPlot.Children(3).Title, 'String', ['Raincloud plot: ' MRSCont.quantify.metabs{gui.SelectedMetab}]) %Update title
end

function osp_updatecorrOvWindow()   %updates the correlation overview tab
            rmpath(genpath([spmversion filesep]));
            delete(gui.corrOvPlot.Children(3).Children)
%%%%%%%%%%%%%%%%%%VISUALIZATION PART OF THIS TAB%%%%%%%%%%%%%%%%%%%%%%%%
            temp = figure( 'Visible', 'off' );
            if gui.SelectedCorrChoice == 1
                temp = osp_plotScatter(MRSCont,gui.QuantModelNames{gui.SelectedModel},gui.QuantNames{gui.SelectedQuant},MRSCont.quantify.metabs{gui.SelectedMetab},gui.CorrMeas{gui.SelectedCorr},gui.CorrNames{gui.SelectedCorr},1);
            else if gui.SelectedCorrChoice == 2
                temp = osp_plotScatter(MRSCont,gui.QuantModelNames{gui.SelectedModel},gui.QuantNames{gui.SelectedQuant},MRSCont.quantify.metabs{gui.SelectedMetab},MRSCont.quantify.metabs{gui.SelectedCorr},MRSCont.quantify.metabs{gui.SelectedCorr},1);
                else
                    switch gui.SelectedCorr
                        case 1
                        temp = osp_plotScatter(MRSCont,gui.QuantModelNames{gui.SelectedModel},gui.QuantNames{gui.SelectedQuant},MRSCont.quantify.metabs{gui.SelectedMetab},MRSCont.QM.SNR.A',gui.QMNames{gui.SelectedCorr},1);
                        case 2
                        temp = osp_plotScatter(MRSCont,gui.QuantModelNames{gui.SelectedModel},gui.QuantNames{gui.SelectedQuant},MRSCont.quantify.metabs{gui.SelectedMetab},MRSCont.QM.FWHM.A',gui.QMNames{gui.SelectedCorr},1);
                    end
                end
            end
            set(gui.corrOvPlot.Children(3), 'XLim', temp.Children(2).XLim);
            set(gui.corrOvPlot.Children(3), 'YLim', temp.Children(2).YLim);
            set(temp.Children(2).Children, 'Parent', gui.corrOvPlot.Children(3) );
            set(gui.corrOvPlot.Children(3), 'XLabel', temp.Children(2).XLabel);
            set(gui.corrOvPlot.Children(3), 'YLabel', temp.Children(2).YLabel);
            close(temp);
            set(gui.corrOvPlot,'Heights', [-0.07 -0.90 -0.03]);
            gui.corrOvPlot.Children(3).Legend.Location = 'North'; %Update legend
            if gui.SelectedCorrChoice == 1
                set(gui.corrOvPlot.Children(3).Title, 'String', [gui.CorrNames{gui.SelectedCorr} ' vs ' MRSCont.quantify.metabs{gui.SelectedMetab}]) %Update title
            else if gui.SelectedCorrChoice == 2
                set(gui.corrOvPlot.Children(3).Title, 'String', [MRSCont.quantify.metabs{gui.SelectedCorr} ' vs ' MRSCont.quantify.metabs{gui.SelectedMetab}]) %Update title
                else
                    switch gui.SelectedCorr
                        case 1
                        set(gui.corrOvPlot.Children(3).Title, 'String', [MRSCont.quantify.metabs{gui.SelectedCorr} ' vs SNR']) %Update title
                        case 2
                        set(gui.corrOvPlot.Children(3).Title, 'String', [MRSCont.quantify.metabs{gui.SelectedCorr} ' vs FHWM (ppm)']) %Update title
                    end
                end
            end
end
%% CALLBACK FUNCTIONS FOR THE LEFT MENU
function onLoad( ~, ~ ) %Callback Load button
    gui.tabs.Selection  = 1;
    osp_processingWindow()
    % User wants load data
    MRSCont = OspreyLoad(MRSCont);
    delete(gui.dummy);
    if ~isempty(MRSCont.raw{1,gui.SelectedDataset}.seq)
        if strcmp(sprintf('\n'),MRSCont.raw{1,gui.SelectedDataset}.seq(end)) %Clean up Sequence Name if needed
            SeqName = MRSCont.raw{1,gui.SelectedDataset}.seq(1:end-1);
        else
            SeqName = MRSCont.raw{1,gui.SelectedDataset}.seq;
        end
    else
            SeqName ='';
    end
    gui.GeometryNames = fieldnames(MRSCont.raw{1,1}.geometry.size); %Get variables regarding voxel geometry
    osp_iniLoadWindow();
    gui.b_load.Enable = 'off';
    gui.b_proc.Enable = 'on';
    gui.b_coreg.Enable = 'on';
    gui.ListBox.Enable = 'on';
end % onLoad

function onProc( ~, ~ ) %Callback Process button
    gui.tabs.Selection  = 2;
    osp_processingWindow()
    % User wants to process the data
    MRSCont = OspreyProcess(MRSCont);
    delete(gui.dummy);
    gui.NoPro = length(fieldnames(MRSCont.processed));
    gui.ProNames = fieldnames(MRSCont.processed);
    gui.ProNames = sort(gui.ProNames);
    osp_iniProcessWindow();
    gui.b_proc.Enable = 'off';
    gui.b_fit.Enable = 'on';
end % onProc

function onFit( ~, ~ ) %Callback Fit button
    gui.tabs.Selection  = 4;
    osp_processingWindow()
    % User wants to quit out of the application
    MRSCont = OspreyFit(MRSCont);
    delete(gui.dummy);    
    gui.NoFits = length(fieldnames(MRSCont.fit.results));
    gui.FitNames = fieldnames(MRSCont.fit.results);
    osp_iniFitWindow();
    gui.b_fit.Enable = 'off';
    gui.b_quant.Enable = 'on';
end % onFit

function onCoreg( ~, ~ ) %Callback Fit button
    if isfield(gui, 'coregResults')
        if length(gui.coregTab.Children) == 3
                delete( gui.coregTab.Children(1) );
                delete( gui.coregTab.Children(1) );
                delete( gui.coregTab.Children(1) );
        else if length(gui.coregTab.Children) == 2
                delete( gui.coregTab.Children(1) );
                delete( gui.coregTab.Children(1) );
            else
                delete( gui.coregTab.Children(1) );
            end
        end
    end
    gui.tabs.Selection  = 3;
    osp_processingWindow()   
    % User wants to quit out of the application
    addpath(genpath([spmversion filesep]));
    MRSCont = OspreyCoreg(MRSCont);    
    rmpath(genpath([spmversion filesep]));
    delete(gui.dummy);     
    osp_iniCoregWindow();
    gui.b_coreg.Enable = 'off';
    gui.b_segm.Enable = 'on';
end % onCoreg

function onSeg( ~, ~ ) %Callback Fit button
    if isfield(gui, 'coregResults')
        if length(gui.coregTab.Children) == 3
                delete( gui.coregTab.Children(1) );
                delete( gui.coregTab.Children(1) );
                delete( gui.coregTab.Children(1) );
        else if length(gui.coregTab.Children) == 2
                delete( gui.coregTab.Children(1) );
                delete( gui.coregTab.Children(1) );
            else
                delete( gui.coregTab.Children(1) );
            end
        end
    end
    gui.tabs.Selection  = 3;
    osp_processingWindow()
    % User wants to quit out of the application
    addpath(genpath([spmversion filesep]));
    MRSCont = OspreySeg(MRSCont);
    rmpath(genpath([spmversion filesep]));
    delete(gui.dummy);     
    osp_iniCoregWindow();       
    gui.b_segm.Enable = 'off';
    if MRSCont.flags.didQuantify
        MRSCont = OspreyQuantify(MRSCont);
        gui.NoQuantModels = length(fieldnames(MRSCont.quantify.tables));
        gui.QuantModelNames = fieldnames(MRSCont.quantify.tables);
        gui.NoQuants = length(fieldnames(MRSCont.quantify.tables.(gui.QuantModelNames{1})));
        gui.QuantNames = fieldnames(MRSCont.quantify.tables.(gui.QuantModelNames{1}));
        gui.NoQuantMetabs = length(MRSCont.quantify.metabs);
        gui.pop_quantOvPlotControls.String = gui.QuantNames;
        gui.pop_distrOvQuantControls.String = gui.QuantNames;
        gui.pop_corrOvQuantControls.String = gui.QuantNames;
    end
end % onSeg

function onQuant( ~, ~ ) %Callback Quantify button
    gui.tabs.Selection  = 5;
    osp_processingWindow()
    % User wants to quit out of the application
    MRSCont = OspreyQuantify(MRSCont);
    delete(gui.dummy);     
    gui.NoQuantModels = length(fieldnames(MRSCont.quantify.tables));
    gui.QuantModelNames = fieldnames(MRSCont.quantify.tables);
    gui.NoQuants = length(fieldnames(MRSCont.quantify.tables.(gui.QuantModelNames{1})));
    gui.QuantNames = fieldnames(MRSCont.quantify.tables.(gui.QuantModelNames{1}));
    gui.NoQuantMetabs = length(MRSCont.quantify.metabs);
    osp_iniQuantifyWindow();
    MRSCont = OspreyOverview(MRSCont);
    gui.NAAnormed = 1;
    gui.NoGroups = MRSCont.overview.NoGroups;
    [cb] = cbrewer('qual', 'Dark2', 12, 'pchip');
    temp = cb(3,:);
    cb(3,:) = cb(4,:);
    cb(4,:) = temp;
    if isfield(MRSCont.overview, 'corr')
        gui.CorrNames = MRSCont.overview.corr.Names{1};
        gui.CorrMeas = MRSCont.overview.corr.Meas;
    end
    gui.SelectedCorr = 1;
    gui.SelectedCorrChoice = 1;
    gui.QMNames = {'SNR','FWHM (ppm)'};
    osp_iniOverviewWindow();
    gui.b_quant.Enable = 'off';
end % onFit

function onExit( ~, ~ ) %Callback Exit button
    % User wants to quit out of the application
    delete( gui.window );
end % onExit
function onListSelection( src, ~) %Callback Listbox with all files
    % User selected file in the list refreshs active tab
    if gui.KeyPress == 0
        gui.SelectedDataset = get( src, 'Value' );
        set(gui.ListBox, 'value', gui.SelectedDataset);
        switch gui.tabs.Selection
            case 1
                osp_updateLoadWindow();
            case 2
                gui.proInfoText = gui.(gui.proTabhandles{gui.SelectedSubFile}).Children(2).Children;
                % Grid for Plot and Data control sliders
                gui.proPlot = gui.(gui.proTabhandles{gui.SelectedSubSpec});
                osp_updateProWindow();
            case 3
                osp_updateCoregWindow();
            case 4
                osp_updateFitWindow();
            case 5
                osp_updateQuantifyWindow();
        end
    else
        gui.SelectedDataset = get(gui.ListBox, 'value');
    end
end % onListSelection

function WindowKeyUp(src,EventData)
    OldValue = get( gui.ListBox,'value');
    gui.KeyPress = 0;
    updateListBox()
end

function WindowKeyDown(src,EventData)
    % 28 leftarrow
    % 29 rightarrow
    % 30 uparrow
    % 31 downarrow
    if strcmp(EventData.Key, 'uparrow')
        OldValue = get( gui.ListBox,'value');
        gui.KeyPress = 1;
        if OldValue == 1
            set(gui.ListBox, 'value', MRSCont.nDatasets );
        else
            set(gui.ListBox, 'value', OldValue-1 );
        end
    end
    if strcmp(EventData.Key, 'downarrow')
        OldValue = get( gui.ListBox,'value');
        gui.KeyPress = 1;
        if OldValue == MRSCont.nDatasets
            set(gui.ListBox, 'value', 1 );
        else
            set(gui.ListBox, 'value', OldValue+1 );
        end
    end
end
%% CALLBACK FUNCTIONS FOR THE TAB
function SelectionChangedFcn(varargin) %Callback for the main tab panel
   % User selected tab refreshs plot
   OldValue = varargin{1,2}.OldValue;
   switch OldValue
       case 1
            spec = gui.SelectedSubFile;
            if ~(spec == 2 || spec ==3)
                gui.SelectedFit = 1;
            else
                if spec == 2
                    gui.SelectedFit = 2;
                else
                    gui.SelectedFit = 3;
                end
            end
       case 2
            spec = gui.SelectedSubSpec;

            if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                if ~(spec == 5 || spec == 6)
                    gui.SelectedSubFile = 1;
                    gui.SelectedFit = 1;
                else
                    if spec == 5
                        gui.SelectedSubFile = 2;
                        gui.SelectedFit = 4;
                    else
                        gui.SelectedSubFile = 3;
                        gui.SelectedFit = 5;
                    end
                end
            elseif MRSCont.flags.isMEGA
                if ~(spec == 3 || spec == 4)
                    gui.SelectedSubFile = 1;
                    gui.SelectedFit = 1;
                else
                    if spec == 3
                        gui.SelectedSubFile = 2;
                        gui.SelectedFit = 3;
                    else
                        gui.SelectedSubFile = 3;
                        gui.SelectedFit = 4;
                    end
                end
            elseif MRSCont.flags.isUnEdited
                        gui.SelectedSubFile = spec;
                        gui.SelectedFit = spec;
            end

       case 4
            if (strcmp(gui.FitNames{gui.SelectedFit},'off') || strcmp(gui.FitNames{gui.SelectedFit},'diff1')||...
                strcmp(gui.FitNames{gui.SelectedFit},'diff2') || strcmp(gui.FitNames{gui.SelectedFit},'sum'))
                spec = 1;
            else
                if strcmp(gui.FitNames{gui.SelectedFit},'w')
                    spec = 3;
                else
                    spec = 2;
                end
            end
            if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                if ~(spec == 2 || spec == 3)
                    gui.SelectedSubFile = 1;
                    gui.SelectedSubSpec = 1;
                else
                    if spec == 2
                        gui.SelectedSubFile = 2;
                        gui.SelectedSubSpec = 5;
                    else
                        gui.SelectedSubFile = 3;
                        gui.SelectedSubSpec = 6;
                    end
                end
            elseif MRSCont.flags.isMEGA
                if ~(spec == 2 || spec == 3)
                    gui.SelectedSubFile = 1;
                    gui.SelectedSubSpec = 1;
                else
                    if spec == 2
                        gui.SelectedSubFile = 2;
                        gui.SelectedSubSpec = 3;
                    else
                        gui.SelectedSubFile = 3;
                        gui.SelectedSubSpec = 4;
                    end
                end
            elseif MRSCont.flags.isUnEdited
                        gui.SelectedSubFile = spec;
                        gui.SelectedSubSpec = spec;
            end
       case 5
       otherwise
           idx = 1;
           spec = 1;
   end


   switch gui.tabs.Selection
        case 1
        gui.ListBox.Enable = 'on';
        gui.dataInfoText = gui.(gui.rawTabhandles{gui.SelectedSubFile}).Children(2).Children;
        gui.dataPlot = gui.(gui.rawTabhandles{gui.SelectedSubFile});
        set(gui.rawTab, 'selection', gui.SelectedSubFile);
        case 2
        gui.ListBox.Enable = 'on';
        gui.proInfoText = gui.(gui.proTabhandles{gui.SelectedSubFile}).Children(2).Children;
        gui.proPlot = gui.(gui.proTabhandles{gui.SelectedSubSpec});
        set(gui.proTab, 'selection', gui.SelectedSubSpec);
        osp_updateProWindow()
        case 3
        gui.ListBox.Enable = 'on';
        gui.coregInfoText = gui.(gui.proTabhandles{gui.SelectedSubFile}).Children(2).Children;
        case 4
        gui.ListBox.Enable = 'on';
        gui.fitInfoText = gui.(gui.fitTabhandles{gui.SelectedFit}).Children(2).Children;
        gui.fitPlot = gui.(gui.fitTabhandles{gui.SelectedFit});
        set(gui.fitTab, 'selection', gui.SelectedFit);
        osp_updateFitWindow()
        case 5
        gui.ListBox.Enable = 'on';
        osp_updateQuantifyWindow()
       case 6
        gui.ListBox.Enable = 'off';
    end
end
%% CALLBACK FUNCTIONS FOR THE LOAD TAB

function RawTabChangeFcn(varargin) %Callback for tab changes in the raw tab
   % User selected tab refreshs plot
    gui.SelectedSubFile = varargin{1,2}.NewValue;
    % Parameter shown in the info panel on top
    gui.dataInfoText = gui.(gui.rawTabhandles{gui.SelectedSubFile}).Children(2).Children;
    % Grid for Plot and Data control sliders
    gui.dataPlot = gui.(gui.rawTabhandles{gui.SelectedSubFile});
    osp_updateLoadWindow();
end
%% CALLBACK FUNCTIONS FOR THE PROCESS TAB

function ProTabChangeFcn(varargin) %Callback for tab changes in the raw tab
    % User selected tab refreshs plot
    gui.SelectedSubSpec = varargin{1,2}.NewValue;
    % Parameter shown in the info panel on top
    gui.proInfoText = gui.(gui.proTabhandles{gui.SelectedSubSpec}).Children(2).Children;
    % Grid for Plot and Data control sliders gui.(gui.proTabhandles{ind})
    gui.proPlot = gui.(gui.proTabhandles{gui.SelectedSubSpec});
    gui.proPre = gui.proPlot.Children(1).Children(2).Children(2);
    gui.proPost = gui.proPlot.Children(1).Children(2).Children(1);
    gui.proDrift = gui.proPlot.Children(1).Children(1).Children(2);
    gui.proAlgn = gui.proPlot.Children(1).Children(1).Children(1);
    osp_updateProWindow();
end
%% CALLBACK FUNCTIONS FOR THE FIT TAB

function FitTabChangeFcn(varargin) %Callback for tab changes in the raw tab
    % User selected tab refreshs plot
    gui.SelectedFit = varargin{1,2}.NewValue;
    % Parameter shown in the info panel on top
    gui.fitInfo = gui.(gui.fitTabhandles{gui.SelectedFit}).Children(2);
    gui.fitResults =  gui.(gui.fitTabhandles{gui.SelectedFit}).Children(1).Children(1);
    gui.fitPlot = gui.(gui.fitTabhandles{gui.SelectedFit});
    gui.fitInfoText = gui.(gui.fitTabhandles{gui.SelectedFit}).Children(2).Children;
    gui.fitResultsText = gui.(gui.fitTabhandles{gui.SelectedFit}).Children(1).Children(1).Children;
    delete(gui.fitPlot.Children(1).Children(2).Children)
    % Grid for Plot and Data control sliders
    gui.fitPlot = gui.(gui.fitTabhandles{gui.SelectedFit});
    osp_updateFitWindow();
end

%% CALLBACK FUNCTIONS FOR OVERVIEW TAB
function pop_specsOvPlot_Call(h,~) %Callback for the popup menu of subspectra in the spectra overview tab
        idx=(h.Value);
        h.Value=idx;
        gui.SelectedSubSpec = idx;
        osp_updateSpecsOvWindow()
end % pop_specsOvPlot_Call

function pop_meanOvPlot_Call(h,~)  %Callback for the popup menu of subspectra in the mean sd overview tab
        idx=(h.Value);
        h.Value=idx;
        gui.SelectedSubSpec = idx;
        osp_updatemeanOvWindow()
end % pop_meanOvPlot_Call

function pop_quantOvPlot_Call(h,~)  %Callback for the popup menu of quantifications the quantification overview tab
        idx=(h.Value);
        h.Value=idx;
        gui.SelectedQuant = idx;
        osp_updatequantOvWindow()
end % pop_quantOvPlot_Call

function pop_distrOvQuant_Call(h,~)  %Callback for the popup menu of quantifications in the distribution tab
        idx=(h.Value);
        h.Value=idx;
        gui.SelectedQuant = idx;
        osp_updatedistrOvWindow()
end % pop_distrOvQuant_Call

function pop_distrOvMetab_Call(h,~)  %Callback for the popup menu of metabolites in the distribution tab
        idx=(h.Value);
        h.Value=idx;
        gui.SelectedMetab = idx;
        osp_updatedistrOvWindow()
end % pop_corrOvMetabControls

function pop_corrOvQuant_Call(h,~) %Callback for the popup menu of quantifications in the distribution tab
        idx=(h.Value);
        h.Value=idx;
        gui.SelectedQuant = idx;
        osp_updatecorrOvWindow()
end % pop_corrOvQuant_Call

function pop_corrOvMetab_Call(h,~) %Callback for the popup menu of metabolites in the distribution tab
        idx=(h.Value);
        h.Value=idx;
        gui.SelectedMetab = idx;
        osp_updatecorrOvWindow()
end % pop_corrOvMetab_Call

function pop_corrOvCorr_Call(h,~) %Callback for the popup menu of correlations in the correlation tab
        idx=(h.Value);
        h.Value=idx;
        gui.SelectedCorr = idx;
        osp_updatecorrOvWindow()
end % pop_corrOvCorr_Call

function pop_whichcorrOvCorr_Call(h,~) %Callback for the popup menu of correlations in the correlation tab
        idx=(h.Value);
        h.Value=idx;
        gui.SelectedCorrChoice = idx;
        gui.SelectedCorr = 1;
        if idx == 1
            set(gui.pop_corrOvCorrControls, 'String', gui.CorrNames);
            set(gui.pop_corrOvCorrControls, 'Value', gui.SelectedCorr);
        else if idx == 2
            set(gui.pop_corrOvCorrControls, 'String', MRSCont.quantify.metabs);
            set(gui.pop_corrOvCorrControls, 'Value', gui.SelectedCorr);
            else
                    set(gui.pop_corrOvCorrControls, 'String', gui.QMNames);
                set(gui.pop_corrOvCorrControls, 'Value', gui.SelectedCorr);
            end
        end
        osp_updatecorrOvWindow()
end % pop_corrOvCorr_Call

function OverviewTabChangedFcn(varargin) %Callback for the selection changes in the overview tabs
    OldValue= varargin{1,2}.OldValue;
    NewValue= varargin{1,2}.NewValue;
    switch NewValue
       case 1
            osp_updateSpecsOvWindow()
            set(gui.pop_specsOvPlotControls, 'value',gui.SelectedSubFile)
       case 2
            osp_updatemeanOvWindow()
            set(gui.pop_meanOvPlotControls, 'value',gui.SelectedSubFile)
       case 4
            osp_updatedistrOvWindow()
            set(gui.pop_distrOvQuantControls, 'value',gui.SelectedQuant)
            set(gui.pop_distrOvMetabControls, 'value',gui.SelectedMetab)
       case 5
            osp_updatecorrOvWindow()
            set(gui.pop_corrOvQuantControls, 'value',gui.SelectedQuant)
            set(gui.pop_corrOvMetabControls, 'value',gui.SelectedMetab)
      otherwise
            set(overviewTab, 'selection', 1);
    end
end
end
