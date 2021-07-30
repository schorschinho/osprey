function [gui] = osp_onLoad( ~, ~ ,gui)
%% osp_onLoad
%   Callback function on load button click.
%
%
%   USAGE:
%       osp_onLoad( ~, ~ ,gui);
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
%%% 1. INITIALIZE DATA %%%
    MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
    set(gui.figure,'HandleVisibility','off');
    gui.layout.tabs.Selection  = 1;
    [gui,MRSCont] = osp_processingWindow(gui,MRSCont);
%%% 2. CALL OSPREYLOAD %%%
    MRSCont = OspreyLoad(MRSCont);
    if MRSCont.flags.isPRIAM
        try
            gui.info.nXvoxels = MRSCont.raw{1,gui.controls.Selected}.nXvoxels;
            gui.info.nYvoxels = MRSCont.raw{1,gui.controls.Selected}.nYvoxels;
            gui.info.nZvoxels = MRSCont.raw{1,gui.controls.Selected}.nZvoxels;
        catch
        end
    end
    if MRSCont.flags.isMRSI
        try
            gui.info.nXvoxels = MRSCont.raw{1,gui.controls.Selected}.nXvoxels;
            gui.info.nYvoxels = MRSCont.raw{1,gui.controls.Selected}.nYvoxels;
            gui.info.nZvoxels = MRSCont.raw{1,gui.controls.Selected}.nZvoxels;
        catch
        end
    end
    delete(gui.layout.dummy);
    if ~isempty(MRSCont.raw{1,gui.controls.Selected}.seq)
        if strcmp(sprintf('\n'),MRSCont.raw{1,gui.controls.Selected}.seq(end)) %Clean up Sequence Name if needed
            gui.load.Names.Seq = MRSCont.raw{1,gui.controls.Selected}.seq(1:end-1);
        else
            gui.load.Names.Seq = MRSCont.raw{1,gui.controls.Selected}.seq;
        end
    else
            gui.load.Names.Seq ='';
    end
    gui.load.Names.Geom = fieldnames(MRSCont.raw{1,1}.geometry.size); %Get variables regarding voxel geometry
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
    set(gui.figure,'HandleVisibility','on');
%%% 3. INITIALIZE OUTPUT WINDOW %%%     
    osp_iniLoadWindow(gui);
    gui.layout.b_load.Enable = 'off';
    gui.layout.b_proc.Enable = 'on';
    if MRSCont.flags.hasSPM == 1
        gui.layout.b_coreg.Enable = 'on';
    end
    gui.layout.ListBox.Enable = 'on';

end % onLoad