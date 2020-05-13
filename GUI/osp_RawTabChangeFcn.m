function osp_RawTabChangeFcn(src,~,gui)
%% osp_RawTabChangeFcn
%   This function is triggered when the raw tab is changed. It refreshes
%   the figure.
%
%
%   USAGE:
%       osp_RawTabChangeFcn(src,~,gui);
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
    gui.load.Selected = src.Selection;
    % Parameter shown in the info panel on top
    gui.InfoText.data = gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children;
    % Grid for Plot and Data control sliders
    gui.Plot.data = gui.layout.(gui.layout.rawTabhandles{gui.load.Selected});
%%% 2. UPDATE GUI %%%    
    osp_updateLoadWindow(gui);
end