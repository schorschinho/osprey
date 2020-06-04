function [gui, MRSCont] = osp_processingWindow(gui,MRSCont)
%% osp_processingWindow
%This functions finds the handle to update the window in progress
%
%
%   USAGE:
%       osp_processingWindow(gui);
%
%   INPUT:  
%           gui      = gui class containing all handles and the MRSCont
%           MRSCont      = MRS container
%
%   OUTPUT:
%           gui      = gui class containing all handles and the MRSCont
%           MRSCont      = MRS container
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

    switch gui.layout.tabs.Selection
        case 1
            gui.layout.tabs.TabEnables{1} = 'on';
            gui.layout.dummy = uix.VBox('Parent', gui.layout.rawTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.layout.dummyTxt = uicontrol('Parent',gui.layout.dummy,'style','text', ...
                                          'FontSize', 16, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
            gui.layout.rawTab.TabTitles  = {'...'};
        case 2 
            gui.layout.tabs.TabEnables{2} = 'on';
            gui.layout.dummy = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.layout.dummyTxt = uicontrol('Parent',gui.layout.dummy,'style','text', ...
                                          'FontSize', 16, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
            gui.layout.proTab.TabTitles  = {'...'};
        case 3 
            gui.layout.tabs.TabEnables{3} = 'on';
            gui.layout.dummy = uix.VBox('Parent', gui.layout.fitTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.layout.dummyTxt = uicontrol('Parent',gui.layout.dummy,'style','text', ...
                                          'FontSize', 16, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
            gui.layout.fitTab.TabTitles  = {'...'};
        case 4 
            gui.layout.tabs.TabEnables{4} = 'on';
            gui.layout.dummy = uix.VBox('Parent', gui.layout.coregTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.layout.dummyTxt = uicontrol('Parent',gui.layout.dummy,'style','text', ...
                                          'FontSize', 16, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
        case 5 
            gui.layout.tabs.TabEnables{5} = 'on';
            gui.layout.dummy = uix.VBox('Parent', gui.layout.quantifyTab, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
            gui.layout.dummyTxt = uicontrol('Parent',gui.layout.dummy,'style','text', ...
                                          'FontSize', 16, 'FontName', 'Arial','ForegroundColor', gui.colormap.Foreground,...
                                          'HorizontalAlignment', 'left', 'String', '', 'BackgroundColor',gui.colormap.Background);
    end
    set(gui.layout.dummyTxt, 'String', sprintf('In progress...'));  
    MRSCont.flags.inProgress = gui.layout.dummyTxt;
end