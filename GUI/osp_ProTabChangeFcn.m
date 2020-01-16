function osp_ProTabChangeFcn(src,~,gui)
%% osp_ProTabChangeFcn
%   This function is triggered when the process tab is changed. It refreshes
%   the figure.
%
%
%   USAGE:
%       osp_ProTabChangeFcn(src,~,gui);
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
    gui.process.Selected = src.Selection;
    % Parameter shown in the info panel on top
    gui.InfoText.pro = gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children;
    % Grid for Plot
    gui.layout.proPlot = gui.layout.(gui.layout.proTabhandles{gui.process.Selected});
    gui.layout.proPre = gui.layout.proPlot.Children(1).Children(2).Children(2);
    gui.layout.proPost = gui.layout.proPlot.Children(1).Children(2).Children(1);
    gui.layout.proDrift = gui.layout.proPlot.Children(1).Children(1).Children(2);
    gui.layout.proAlgn = gui.layout.proPlot.Children(1).Children(1).Children(1);
%%% 2. UPDATE GUI %%%    
    osp_updateProWindow(gui);
end
