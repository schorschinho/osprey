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
    MRSCont = getappdata(gui.figure,'MRSCont');  
    
    gui.fit.Selected = src.Selection;
    gui.info.nXvoxels = 1;
    gui.controls.act_x = 1;
    gui.controls.b_left_y.Enable = 'off';
    gui.controls.b_left_z.Enable = 'off';
    gui.controls.b_right_y.Enable = 'off';
    gui.controls.b_right_z.Enable = 'off';
    gui.info.nYvoxels = size(MRSCont.fit.results.(gui.fit.Names{gui.fit.Selected}).fitParams,3);
    gui.controls.act_y = 1;
    gui.info.nZvoxels = size(MRSCont.fit.results.(gui.fit.Names{gui.fit.Selected}).fitParams,1);
    gui.controls.act_z = 1;
    if size(MRSCont.fit.results.(gui.fit.Names{gui.fit.Selected}).fitParams,3) > 1
        gui.controls.b_left_y.Enable = 'on';
        gui.controls.b_right_y.Enable = 'on';
    end
    if size(MRSCont.fit.results.(gui.fit.Names{gui.fit.Selected}).fitParams,1) > 1
        gui.controls.b_left_z.Enable = 'on';
        gui.controls.b_right_z.Enable = 'on';
    end
%%% 2. UPDATE GUI %%%    
    osp_updateFitWindow(gui);
end