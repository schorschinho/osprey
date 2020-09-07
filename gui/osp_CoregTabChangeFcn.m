function osp_CoregTabChangeFcn(src,~,gui)
%% osp_CoregTabChangeFcn
%   This function is triggered when the coreg tab is changed. It refreshes
%   the figure.
%
%
%   USAGE:
%       osp_CoregTabChangeFcn(src,~,gui);
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
    gui.coreg.Selected = src.Selection;
    % Parameter shown in the info panel on top
%     gui.Info.coreg = gui.layout.(gui.layout.coregTabhandles{gui.coreg.Selected}).Children(2);
    gui.Results.coreg =  gui.layout.(gui.layout.coregTabhandles{gui.coreg.Selected}).Children(1).Children(1);
    gui.InfoText.coreg = gui.layout.(gui.layout.coregTabhandles{gui.coreg.Selected}).Children(2).Children;
%     gui.Results.coreg = gui.layout.(gui.layout.coregTabhandles{gui.coreg.Selected}).Children(1).Children(1).Children;
%     delete(gui.Plot.coreg.Children(1).Children(2).Children)
    % Grid for Plot
    gui.Plot.coreg = gui.layout.(gui.layout.coregTabhandles{gui.coreg.Selected});

%%% 2. UPDATE GUI %%%    
    osp_updateCoregWindow(gui);
end