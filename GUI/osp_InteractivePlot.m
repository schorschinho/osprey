function osp_InteractivePlot(MRSCont)
%demoBrowser: an example of using layouts to build a user interface
%
%   demoBrowser() opens a simple GUI that allows several of MATLAB's
%   built-in demos to be viewed. It aims to demonstrate how multiple
%   layouts can be used to create a good-looking user interface that
%   retains the correct proportions when resized. It also shows how to
%   hook-up callbacks to interpret user interaction.
%
%   See also: <a href="matlab:doc Layouts">Layouts</a>

%   Copyright 2010-2013 The MathWorks, Inc.

% Data is shared between all child functions by declaring the variables
% here (they become global to the function). We keep things tidy by putting
% all GUI stuff in one structure and all data stuff in another. As the app
% grows, we might consider making these objects rather than structures.
data = createData(MRSCont);
gui = createInterface( data.FileNames);

% Now update the GUI with the current data
updateInterface();
redrawDemo();

% Explicitly call the demo display so that it gets included if we deploy
displayEndOfDemoMessage('')

%-------------------------------------------------------------------------%
    function data = createData(MRSCont)
        % Create the shared data-structure for this application
        fileList = {MRSCont.files
            };
        selectedFile = 1;
        selectedProcess = 1;
        data = struct( ...
            'FileNames',fileList, ...
            'SelectedFile', selectedFile, ...
            'SelectedProcess', selectedProcess);
    end % createData

%-------------------------------------------------------------------------%
    function gui = createInterface( fileList )
        % Create the user interface for the application and return a
        % structure of handles for global use.
        gui = struct();
        % Open a window and add some menus
        gui.Window = figure( ...
            'Name', 'Gallery browser', ...
            'NumberTitle', 'off', ...
            'MenuBar', 'none', ...
            'Toolbar', 'none', ...
            'HandleVisibility', 'off' );
        
        % + File menu
        gui.FileMenu = uimenu( gui.Window, 'Label', 'File' );
        uimenu( gui.FileMenu, 'Label', 'Exit', 'Callback', @onExit );
        
        % + View menu
        gui.ViewMenu = uimenu( gui.Window, 'Label', 'View' );
        uimenu( gui.ViewMenu, 'Label', 'Load', 'Callback', @onMenuSelection );
        uimenu( gui.ViewMenu, 'Label', 'Process', 'Callback', @onMenuSelection );
        uimenu( gui.ViewMenu, 'Label', 'Fit', 'Callback', @onMenuSelection );

        
    
        
        % Arrange the main interface
        mainLayout = uix.HBoxFlex( 'Parent', gui.Window, 'Spacing', 3 );
        
        % + Create the panels
        controlPanel = uix.BoxPanel( ...
            'Parent', mainLayout, ...
            'Title', 'Select a File:' );
        gui.ViewPanel = uix.BoxPanel( ...
            'Parent', mainLayout, ...
            'Title', 'Viewing: ???', ...
            'HelpFcn', @onDemoHelp );
        gui.ViewContainer = uicontainer( ...
            'Parent', gui.ViewPanel );        

        % + Adjust the main layout
        set( mainLayout, 'Widths', [-1,-2]  );
        
        
        % + Create the controls
        controlLayout = uix.VBox( 'Parent', controlPanel, ...
            'Padding', 3, 'Spacing', 3 );
        gui.ListBox = uicontrol( 'Style', 'list', ...
            'BackgroundColor', 'w', ...
            'Parent', controlLayout, ...
            'String', fileList(:), ...
            'Value', 1, ...
            'Callback', @onListSelection);
        
        % + Create the view
        p = gui.ViewContainer;
        gui.ViewAxes = axes( 'Parent', p );
        
        
    end % createInterface

%-------------------------------------------------------------------------%
    function updateInterface()
        % Update various parts of the interface in response to the demo
        % being changed.
        
        % Update the list and menu to show the current demo
        set( gui.ListBox, 'Value', data.SelectedFile );
        fileName = data.FileNames{ data.SelectedFile };
        % Update the help button label
        % Update the view panel title
        set( gui.ViewPanel, 'Title', sprintf( 'Viewing: %s', fileName ) );
        % Untick all menus
        menus = get( gui.ViewMenu, 'Children' );
        set( menus, 'Checked', 'off' );
        % Use the name to work out which menu item should be ticked
        whichMenu = strcmpi( fileName, get( menus, 'Label' ) );
        set( menus(whichMenu), 'Checked', 'on' );
    end % updateInterface

%-------------------------------------------------------------------------%
    function redrawDemo()
        % Draw a demo into the axes provided
        
        % We first clear the existing axes ready to build a new one
        if ishandle( gui.ViewAxes )
            delete( gui.ViewAxes );
        end
        if ishandle( gui.ViewAxes )
            delete( gui.ViewAxes );
        end
        % Some demos create their own figure. Others don't.
        fileName = data.FileNames{data.SelectedFile};
        switch upper( data.SelectedProcess )
            case '1'
                % These demos open their own windows
                evalin( 'base', fcnName );
                gui.ViewAxes = gca();
                fig = gcf();
                set( fig, 'Visible', 'off' );
                
            case '2'
                % These demos open their own windows
                evalin( 'base', fcnName );
                gui.ViewAxes = gca();
                fig = gcf();
                set( fig, 'Visible', 'off' );
                
            case '3'
                % These demos open their own windows
                evalin( 'base', fcnName );
                gui.ViewAxes = gca();
                fig = gcf();
                set( fig, 'Visible', 'off' );
                
            otherwise
                % These demos need a window opening
                fig = figure( 'Visible', 'off' );
                evalin( 'base', fcnName );
                gui.ViewAxes = gca();
        end
        % Now copy the axes from the demo into our window and restore its
        % state.
        cmap = colormap( gui.ViewAxes );
        set( gui.ViewAxes, 'Parent', gui.ViewContainer );
        colormap( gui.ViewAxes, cmap );
        rotate3d( gui.ViewAxes, 'on' );
        % Get rid of the demo figure
        close( fig );
    end % redrawDemo

%-------------------------------------------------------------------------%
    function onListSelection( src, ~ )
        % User selected a demo from the list - update "data" and refresh
        data.SelectedFile = get( src, 'Value' );
        updateInterface();
        redrawDemo();
    end % onListSelection

%-------------------------------------------------------------------------%
    function onMenuSelection( src, ~ )

        data.selectedProcess = get( src, 'Label' );
        updateInterface();
        redrawDemo();
    end % onMenuSelection

%-------------------------------------------------------------------------%
    function onDemoHelp( ~, ~ )
        % User wnats documentation for the current demo
        showdemo( data.DemoFunctions{data.SelectedDemo} );
    end % onDemoHelp

%-------------------------------------------------------------------------%
    function onExit( ~, ~ )
        % User wants to quit out of the application
        delete( gui.Window );
    end % onExit

end % EOF