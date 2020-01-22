function osp_iniCoregWindow(gui)
%% osp_iniCoregWindow
%   This function creates the inital coregistration window in the gui.
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
% This function creates the inital coreg/seg window
    MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
    addpath(genpath([gui.folder.spmversion filesep])); % Add SPM path
    gui.layout.tabs.TabEnables{4} = 'on';
    gui.layout.tabs.Selection  = 4;
    gui.layout.EmptyQuantPlot = 0;
    
%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
    gui.Info.coreg = uix.Panel('Parent', gui.layout.coregTab, 'Padding', 5, ...
                              'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}],...
                              'FontName', 'Arial','HighlightColor', gui.colormap.Foreground,'BackgroundColor',gui.colormap.Background,...
                              'ForegroundColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
    % Creates layout for plotting and data control
    gui.Plot.coreg = uix.HBox('Parent', gui.layout.coregTab,'BackgroundColor',gui.colormap.Background);
    set(gui.layout.coregTab, 'Heights', [-0.1 -0.9]);
    % Get parameter from file to fill the info panel

    StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                 '\nraw subspecs: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                 '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                 num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];

   gui.InfoText.coreg  = uicontrol('Parent',gui.Info.coreg,'style','text',...
        'FontSize', 12, 'FontName', 'Arial','HorizontalAlignment', 'left', 'String', sprintf(StatText),...
    'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

%%% 2.VISUALIZATION PART OF THIS TAB %%%
% In this case osp_plotCoreg or osp_plotSegment is used to visualize the
% coregistration or the segmentation
    gui.Results.coreg = uix.VBox('Parent', gui.Plot.coreg,'BackgroundColor',gui.colormap.Background);
    temp = figure( 'Visible', 'off' );
    if MRSCont.flags.didSeg %Did segment. In this case coreg has already been performed. Visualize both
        osp_plotCoreg(MRSCont, gui.controls.Selected, 1);
        ViewAxes = gca();
        set(ViewAxes, 'Parent', gui.Results.coreg );
        colormap(gui.Results.coreg.Children,'gray')
        close( temp );
        temp = figure( 'Visible', 'off' );
        osp_plotSegment(MRSCont, gui.controls.Selected, 1);
        ViewAxes = gca();
        set(ViewAxes, 'Parent', gui.Results.coreg );
        colormap(gui.Results.coreg.Children(1),'gray');
        close( temp );
    else % Only coreg has been run
        osp_plotCoreg(MRSCont, gui.controls.Selected, 1);
        ViewAxes = gca();
        set(ViewAxes, 'Parent', gui.Results.coreg );
        colormap(gui.Results.coreg.Children,'gray');
        close( temp );
    end

    rmpath(genpath([gui.folder.spmversion filesep])); %Remove SPM path to avoid crash due to stupid naming conventions in SPM with internal gamma function
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end
