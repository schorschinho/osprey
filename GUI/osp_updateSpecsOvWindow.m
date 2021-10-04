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
        delete(gui.Plot.specsOv.Children(2).Children)
%%% 2. VISUALIZATION PART OF THIS TAB %%%
        Selection = gui.controls.pop_specsOvPlot.String(gui.process.Selected);
        if gui.controls.GM == 0
            for g = 1 :  gui.overview.Number.Groups %Loop over groups
                temp = osp_plotOverviewSpec(MRSCont, Selection{1},g, gui.layout.shiftind);
                    ax=get(temp,'Parent');
                    figpl = get(ax,'Parent');
                    copyobj(ax.Children, gui.Plot.specsOv.Children(2));
                    % Get rid of the Load figure
                    close( figpl );
            end
        else
           temp = osp_plotOverviewSpec(MRSCont, Selection{1},'GMean', gui.layout.shiftind);
            ax=get(temp,'Parent');
            figpl = get(ax,'Parent');
            copyobj(ax.Children, gui.Plot.specsOv.Children(2));
            % Get rid of the Load figure
            close( figpl );           
        end
        switch Selection{1}
            case {'A','B','C','D','diff1','diff2','sum','MM','MM_clean'}        
                set(gui.Plot.specsOv.Children(2), 'XLim', [0.2 4.5])
                set(gui.Plot.specsOv.Children(2).Title, 'String', ['Overview ' Selection])
            case {'ref','w','Fit:ref','Fit:w'}
                set(gui.Plot.specsOv.Children(2), 'XLim', [0 2*4.68])
                set(gui.Plot.specsOv.Children(2).Title, 'String', ['Overview ' Selection])
            otherwise
                set(gui.Plot.specsOv.Children(2), 'XLim', [0.2 4.5])
                set(gui.Plot.specsOv.Children(2).Title, 'String', ['Overview ' Selection])
        end
end
