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
    MRSCont = getappdata(gui.figure,'MRSCont');  
    % User selected tab refreshs plot
    gui.process.Selected = src.Selection;
    
    if isfield(MRSCont,'processed')
        if MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.dims.extras > 0
            gui.info.nXvoxels = MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.sz(MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.dims.extras);
        else
            gui.info.nXvoxels = 1;
        end   
    else
        gui.info.nXvoxels = 1;
    end
    gui.controls.act_x = 1;
    if isfield(MRSCont,'processed')
        if MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.dims.subSpecs > 0
            gui.info.nYvoxels = MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.sz(MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.dims.subSpecs);
        else
            gui.info.nYvoxels = MRSCont.nDatasets(2);
        end
    else
        gui.info.nYvoxels = MRSCont.nDatasets(2);
    end
    gui.controls.act_y = 1;
            
    % Grid for Plot
    gui.layout.proPlot = gui.layout.(gui.layout.proTabhandles{gui.process.Selected});
    gui.layout.proPre = gui.layout.proPlot.Children(1).Children(2).Children(2);
    gui.layout.proPost = gui.layout.proPlot.Children(1).Children(2).Children(1);
    gui.layout.proDrift = gui.layout.proPlot.Children(1).Children(1).Children(2);
    gui.layout.proAlgn = gui.layout.proPlot.Children(1).Children(1).Children(1);
%%% 2. UPDATE GUI %%%    
    osp_updateProWindow(gui);
end
