function  OspreyStartUp
%% OspreyStartUp(MRSCont)
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
%

%
%   HISTORY:
%       2019-07-11: First version of the code.

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
    [settingsFolder,~,~] = fileparts(which('OspreySettings.m'));
    allFolders      = strsplit(settingsFolder, filesep);
    ospFolder       = strjoin(allFolders(1:end-1), filesep); % parent folder (= Osprey folder)
    logoFcn = @()imread('osprey.png', 'BackgroundColor', gui.colormap.Background);
    logoBanner = uiw.utility.loadIcon(logoFcn);
    % Here the intro banner is created
    gui.d = uiw.dialog.About('Name', 'Osprey','Version','0.0.1','Date', 'October 6, 2019',...
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
    canvasSize(4)   = screenSize(4) * 0.5;
    canvasSize(3)   = canvasSize(4) * 0.5;
    canvasSize(2)   = (screenSize(4) - canvasSize(4))/2;
    canvasSize(1)   = (screenSize(3) - canvasSize(3))/2;
    set(gui.window, 'Position', canvasSize);
    gui.mainLayout = uix.VBox('Parent', gui.window,'BackgroundColor',gui.colormap.Background);
    gui.b_about = uicontrol(gui.mainLayout,'Style','PushButton');
    set(gui.b_about,'Units','Normalized','Position',[0 0 1 1],'BackgroundColor',gui.colormap.Background);
    [img, ~, ~] = imread('osprey.png', 'BackgroundColor', gui.colormap.Background);
    [img2] = imresize(img, 0.09);
    set(gui.b_about,'CData', img2);
    gui.Buttons = uix.VButtonBox(...
            'Parent', gui.mainLayout, 'Padding', 5, 'Spacing', 5, ...
            'BackgroundColor',gui.colormap.Background);
    set(gui.Buttons, 'ButtonSize', [300 60]);
    % Create Job
    gui.CreateJob = uicontrol('Parent', gui.Buttons,'Style','PushButton','String','Create Job','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.CreateJob,'Units','Normalized','Position',[0.1 0.9 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    % JobFile input button
    gui.LoadJob = uicontrol('Parent', gui.Buttons,'Style','PushButton','String','Load Job File','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.LoadJob,'Units','Normalized','Position',[0.1 0.75 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.LoadJob,'Callback',{@onLoadJob});
    % MRSCont input button
    gui.LoadMRSCont = uicontrol('Parent', gui.Buttons,'Style','PushButton','String','Load MRSCont .mat-File','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.LoadMRSCont,'Units','Normalized','Position',[0.1 0.67 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.LoadMRSCont,'Callback',{@onLoadMRSCont});
    % Exit button
    gui.Exit = uicontrol('Parent', gui.Buttons,'Style','PushButton','String','Exit','ForegroundColor', gui.colormap.Foreground,'BackgroundColor', gui.colormap.Background);
    set(gui.Exit,'Units','Normalized','Position',[0.1 0.67 0.8 0.08], 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'Bold');
    set(gui.Exit,'Callback',{@onExit});
    set( gui.mainLayout, 'Heights', [-0.33  -0.66], 'Spacing', 5 );

%% FUNCTIONS
function loadJob()    
    gui.data.MRSCont = OspreyJob(gui.out{1});
    OspreyGUI(gui.data.MRSCont);
    delete(gui.window);
end
function loadMRSCont()
    gui.data = load(gui.out{1}, 'MRSCont');
    OspreyGUI(gui.data.MRSCont);
    delete(gui.window);
end
%% CALLBACK FUNCTIONS


function onLoadJob( ~, ~ ) %Callback Load Job button
    % User wants to quit out of the application
    gui.out =  uipickfiles('FilterSpec',[ospFolder], 'REFilter', '\.m$|\.csv$');
    loadJob();
end % onLoadJob

function onLoadMRSCont( ~, ~ ) %Callback Load MRSCont button
    % User wants to quit out of the application
    gui.out =  uipickfiles('FilterSpec',[ospFolder filesep '*.mat']);
    loadMRSCont();
end % onLoadMRSCont


function onExit( ~, ~ ) %Callback Exit button
    % User wants to quit out of the application
    delete( gui.window );
end % onExit
end
