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

%%% 2. UPDATE GUI %%%    
    osp_updateFitWindow(gui);
end