function osp_updateSpecsOvWindow(gui)
%% osp_updateSpecsOvWindow
%   This function updates the specs overview tab.
%
%
%   USAGE:
%       osp_updateSpecsOvWindow(gui);
%
%   INPUT:  
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
%%% 1. INITIALIZE %%%
        MRSCont = getappdata(gui.figure,'MRSCont');  % Get MRSCont from hidden container in gui class
        selectedOvTab = get(gui.layout.overviewTab,'Selection');
        set(gui.layout.(gui.layout.overviewTabhandels{selectedOvTab}).Children(1).Children(1).Children(3).Children(1).Children(1).Children(3),'String',gui.controls.act_z)
        set(gui.layout.(gui.layout.overviewTabhandels{selectedOvTab}).Children(1).Children(1).Children(3).Children(1).Children(1).Children(4),'String',gui.controls.act_x)        
        delete(gui.Plot.specsOv.Children(2).Children)
        
%%% 2. VISUALIZATION PART OF THIS TAB %%%
        Selection = gui.controls.pop_specsOvPlot.String(gui.overview.Selected.Spec);
        if gui.controls.GM == 0
            for g = 1 :  gui.overview.Number.Groups %Loop over groups
                temp = osp_plotOverviewSpec(MRSCont, Selection{1},g, gui.layout.shiftind,'Frequency (ppm)','','',gui.controls.act_z);
                    ax=get(temp,'Children');
                    copyobj(ax.Children, gui.Plot.specsOv.Children(2));
            end
        else
           temp = osp_plotOverviewSpec(MRSCont, Selection{1},'GMean', gui.layout.shiftind,'Frequency (ppm)','','',gui.controls.act_z);
            ax=get(temp,'Children');
            copyobj(ax.Children, gui.Plot.specsOv.Children(2));          
        end
        
        which_spec_split = split(Selection);
        if length(which_spec_split) == 3
            spec = which_spec_split{2};
        else
            spec = which_spec_split{1};
        end

        switch spec
            case {'metab','mm'}        
                set(gui.Plot.specsOv.Children(2), 'XLim', [0.2 4.5])
                set(gui.Plot.specsOv.Children(2).Title, 'String', ['Overview ' Selection])
                set(gui.Plot.specsOv.Children(2), 'XTick', ax.XTick)
            case {'ref','w','mm_ref'}
                set(gui.Plot.specsOv.Children(2), 'XLim', [0 2*4.68])
                set(gui.Plot.specsOv.Children(2).Title, 'String', ['Overview ' Selection])
                set(gui.Plot.specsOv.Children(2), 'XTick', ax.XTick)
            otherwise
                set(gui.Plot.specsOv.Children(2), 'XLim', [0.2 4.5])
                set(gui.Plot.specsOv.Children(2).Title, 'String', ['Overview ' Selection])
                set(gui.Plot.specsOv.Children(2), 'XTick', ax.XTick)
        end
        h = findall(groot,'Type','figure');
        for ff = 1 : length(h)
            if ~(strcmp(h(ff).Tag, 'Osprey') ||  strcmp(h(ff).Tag, 'TMWWaitbar'))
                close(h(ff))
            end
        end
end
