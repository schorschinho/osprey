function  Osprey
%% Osprey
%   This creates the inital Osprey window. Here you pick exsisitng
%   jobFiles or MRS container to be loaded into the OspreyGui.
%
%   USAGE:
%       Osprey;
%
%
%   AUTHORS:
%       Helge Zoellner (Johns Hopkins University, 2019-11-07)
%       hzoelln2@jhmi.edu
%

%
%   HISTORY:
%       2019-07-11: First version of the code.
%% Check for available add-ons
  [~] = osp_Toolbox_Check ('OspreyGUI',0);
%% Set up Layout
%Here we set up the color layout
%default colormap
gui.colormap.Background = [255/255 254/255 254/255];
gui.colormap.LightAccent = [110/255 136/255 164/255];
gui.colormap.Foreground = [11/255 71/255 111/255];
gui.colormap.Accent = [254/255 186/255 47/255];
[settingsFolder,~,~] = fileparts(which('OspreySettings.m'));
allFolders      = strsplit(settingsFolder, filesep);
ospFolder       = strjoin(allFolders(1:end-1), filesep); % parent folder (= Osprey folder)
logoFcn = @()imread('osprey.png', 'BackgroundColor', gui.colormap.Background);
logoBanner = uiw.utility.loadIcon(logoFcn);
% Here the intro banner is created
gui.d = uiw.dialog.About('Name', 'Osprey','Version','1.0.0','Date', 'October 6, 2019',...
'Timeout', 3,'CustomText', 'Osprey is provided by Johns Hopkins University.',...
'ContactInfo', 'gabamrs@gmail.com','LogoCData', logoBanner);
%% Create the StartUp Menu
gui.out = 0;
gui.window = figure('Name', 'Osprey Startup Menu', 'NumberTitle', 'off', 'MenuBar', 'none', ...
                    'ToolBar', 'none', 'HandleVisibility', 'off', 'Renderer', 'painters', 'Color', gui.colormap.Background);
    % Resize such that width is 1.2941 * height (1.2941 is the ratio
    % between width and height of standard US letter size (11x8.5 in).
    screenSize      = get(0,'ScreenSize');
    canvasSize      = screenSize;
    canvasSize(4)   = 400;
    canvasSize(3)   = 250;
    canvasSize(2)   = (screenSize(4) - canvasSize(4))/2;
    canvasSize(1)   = (screenSize(3) - canvasSize(3))/2;
    set(gui.window, 'Position', canvasSize);
    gui.mainLayout = uix.VBox('Parent', gui.window,'BackgroundColor',gui.colormap.Background,'Padding',5);
    gui.logo = uix.Panel('Parent', gui.mainLayout, ...
                'Padding', 10,'HighlightColor', gui.colormap.Background,'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Background,'ShadowColor', gui.colormap.Background);
    [img, ~, ~] = imread('osprey.png', 'BackgroundColor', gui.colormap.Background);
    [img2] = imresize(img, 0.15);
    temp = figure('Visible','off');
    temp = imshow(img2);
    ViewAxes = gca();
    set(ViewAxes,'Box','off','XColor', 'none','YColor','none');
    set(ViewAxes, 'Parent', gui.logo);

    %Create color scheme popup 
     gui.panelcolors = uix.Panel('Parent', gui.mainLayout, 'Title','Colormap', ...
                'HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
    gui.colors = uicontrol('Parent',gui.panelcolors,'style','popupmenu',...
                                                'Units', 'Normalized', 'Position', [0 0 1 1],'FontName', 'Arial', ...
                                                'String',{'default','Blue','Dark'}, 'Value', 1);

    % Create buttons 
    gui.Buttons = uix.VButtonBox(...
            'Parent', gui.mainLayout, 'Padding', 5, 'Spacing', 5, ...
            'BackgroundColor',gui.colormap.Background);
    set(gui.Buttons, 'ButtonSize', [300 60]);
    % Create Job
    gui.CreateJob = uicontrol('Parent', gui.Buttons,'Style','PushButton','String','Create Job','ForegroundColor', gui.colormap.Foreground,...
                               'TooltipString', 'Create a OspreyJob .m-file');
    set(gui.CreateJob,'Units','Normalized','Position',[0.1 0.9 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.CreateJob,'Callback',{@onCreateJob});
    % JobFile input button
    gui.LoadJob = uicontrol('Parent', gui.Buttons,'Style','PushButton','String','Load Job file','ForegroundColor', gui.colormap.Foreground,...
                               'TooltipString', 'Load a Job (.m or .csv-file)');
    set(gui.LoadJob,'Units','Normalized','Position',[0.1 0.75 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.LoadJob,'Callback',{@onLoadJob});
    % MRSCont input button
    gui.LoadMRSCont = uicontrol('Parent', gui.Buttons,'Style','PushButton','String','Load MRSCont file','ForegroundColor', gui.colormap.Foreground,...
                               'TooltipString', 'Load a MRSContainer .mat-file');
    set(gui.LoadMRSCont,'Units','Normalized','Position',[0.1 0.67 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.LoadMRSCont,'Callback',{@onLoadMRSCont});
    % Exit button
    gui.Exit = uicontrol('Parent', gui.Buttons,'Style','PushButton','String','Exit','ForegroundColor', gui.colormap.Foreground,...
                               'TooltipString', 'See you soon!');
    set(gui.Exit,'Units','Normalized','Position',[0.1 0.67 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.Exit,'Callback',{@onExit});
    set( gui.mainLayout, 'Heights', [-0.2 -0.10  -0.7] );
%% FUNCTIONS
function loadJob(gui)    
    gui.data.MRSCont = OspreyJob(gui.out{1},1);
    idx =gui.colors.Value;
    switch idx
        case 1
            %Default colormap
            gui.colormap.Background = [255/255 254/255 254/255];
            gui.colormap.LightAccent = [110/255 136/255 164/255];
            gui.colormap.Foreground = [11/255 71/255 111/255];
            gui.colormap.Accent = [254/255 186/255 47/255];
        case 2    
            %Bluemode colormap
            gui.colormap.Background = [110/255 136/255 164/255];
            gui.colormap.LightAccent = [51/255 51/255 51/255];
            gui.colormap.Foreground = [244/255 244/255 242/255];
            gui.colormap.Accent = [217/255 224/255 33/255];
        case 3
            %Darkmode colormmap
            gui.colormap.Background = [71/255 71/255 71/255];
            gui.colormap.LightAccent = [110/255 136/255 164/255];
            gui.colormap.Foreground = [244/255 244/255 242/255];
            gui.colormap.Accent = [217/255 224/255 33/255];
    end
    gui.data.MRSCont.colormap = gui.colormap;
    gui.data.MRSCont.colormapidx = idx;
    gui.data.MRSCont.flags.isToolChecked = 1;
    OspreyGUI(gui.data.MRSCont);
    delete(gui.window);
end

function loadMRSCont(gui)
    gui.data = load(gui.out{1}, 'MRSCont');
    idx =gui.colors.Value;
    switch idx
        case 1
            %Default colormap
            gui.colormap.Background = [255/255 254/255 254/255];
            gui.colormap.LightAccent = [110/255 136/255 164/255];
            gui.colormap.Foreground = [11/255 71/255 111/255];
            gui.colormap.Accent = [254/255 186/255 47/255];
        case 2    
            %Bluemode colormap
            gui.colormap.Background = [110/255 136/255 164/255];
            gui.colormap.LightAccent = [51/255 51/255 51/255];
            gui.colormap.Foreground = [244/255 244/255 242/255];
            gui.colormap.Accent = [217/255 224/255 33/255];
        case 3
            %Darkmode colormmap
            gui.colormap.Background = [71/255 71/255 71/255];
            gui.colormap.LightAccent = [110/255 136/255 164/255];
            gui.colormap.Foreground = [244/255 244/255 242/255];
            gui.colormap.Accent = [217/255 224/255 33/255];
    end
    gui.data.MRSCont.colormap = gui.colormap;
    gui.data.MRSCont.flags.isToolChecked = 1;
    OspreyGUI(gui.data.MRSCont);
    delete(gui.window);
end

function onCreateJob( ~, ~) %Callback Create Job button
    % User wants to quit out of the application
    CreateOspreyJob_app
    delete(gui.window);
end % onCreateJob

function onLoadJob( ~, ~) %Callback Load Job button
    % User wants to quit out of the application
    gui.out =  uipickfiles('FilterSpec',[ospFolder], 'REFilter', '\.m$|\.csv$','NumFiles',1,'Prompt','Select Osprey job file (.m or .csv-file)');
    if iscell(gui.out)
        loadJob(gui);
    end
end % onLoadJob

function onLoadMRSCont( ~, ~) %Callback Load MRSCont button
    % User wants to quit out of the application
    gui.out =  uipickfiles('FilterSpec',[ospFolder filesep '*.mat'],'NumFiles',1,'Prompt','Select Osprey MRSCont file (.mat-file)');
    if iscell(gui.out)
        loadMRSCont(gui);
    end
end % onLoadMRSCont


function onExit( ~, ~ ) %Callback Exit button
    % User wants to quit out of the application
    delete( gui.window );
end % onExit
end
