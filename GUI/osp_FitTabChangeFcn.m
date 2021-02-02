function osp_FitTabChangeFcn(src,~,gui)
%% osp_FitTabChangeFcn
%   This function is triggered when the fit tab is changed. It refreshes
%   the figure.
%
%
%   USAGE:
%       osp_FitTabChangeFcn(src,~,gui);
%
%   INPUT:  src      = handle of the fit tabs
%           gui      = gui class containing all handles and the MRSCont             
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
    % User selected tab refreshs plot
    gui.fit.Selected = src.Selection;
    % Parameter shown in the info panel on top
%     gui.Info.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2);
    gui.Results.fit =  gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(1).Children(1);
    gui.Plot.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected});
    gui.InfoText.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children;
    set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(3).Children(1).Children.Children(4),'String',gui.controls.act_z)
    set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(3).Children(1).Children.Children(5),'String',gui.controls.act_y)
    set(gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children(3).Children(1).Children.Children(6),'String',gui.controls.act_x)
    gui.Results.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(1).Children(1).Children;
    delete(gui.Plot.fit.Children(1).Children(2).Children)
    % Grid for Plot
    gui.Plot.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected});

%%% 2. UPDATE GUI %%%    
    osp_updateFitWindow(gui);
end