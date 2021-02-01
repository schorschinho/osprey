function osp_updatemeanOvWindow(gui)
%% osp_updatemeanOvWindow
%   This function updates the mean spectra overview tab.
%
%
%   USAGE:
%       osp_updatemeanOvWindow(gui);
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
        delete(gui.Plot.meanOv.Children(2).Children)
%%% 2. VISUALIZATION PART OF THIS TAB %%%
        Selection = gui.controls.pop_meanOvPlot.String(gui.process.Selected);
        if gui.controls.GM == 0
            for g = 1 :  gui.overview.Number.Groups
                if gui.overview.Number.Groups > 1
                    temp = osp_plotMeanSpec(MRSCont, Selection{1},g,1,1/gui.overview.Number.Groups);
                    ViewAxes = gca();
                    set(ViewAxes.Children,'Parent',gui.Plot.meanOv.Children(2))
                    if g < gui.overview.Number.Groups
                        close(temp);
                    end
                else
                    temp = osp_plotMeanSpec(MRSCont, Selection{1},1, g);
                    ViewAxes = gca();
                    set(ViewAxes.Children,'Parent',gui.Plot.meanOv.Children(2))
                end
            end
        else
            if ~(isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM)
                temp = osp_plotMeanSpec(MRSCont, Selection{1},'GMean', 1);
            else
                temp = osp_plotMeanSpec(MRSCont, Selection{1},1, 0.01,10);
            end
            ViewAxes = gca();
            set(ViewAxes.Children,'Parent',gui.Plot.meanOv.Children(2))
        end
        if (strcmp(Selection{1},'A') || strcmp(Selection{1},'B') || strcmp(Selection{1},'C') || strcmp(Selection{1},'D') || strcmp(Selection{1},'diff1') || strcmp(Selection{1},'diff2') || strcmp(Selection{1},'sum'))
            set(gui.Plot.meanOv.Children(2), 'XLim', ViewAxes.XLim)
            set(gui.Plot.meanOv.Children(2), 'XMinorTick', 'On')
            set(gui.Plot.meanOv.Children(2).Title, 'String', ViewAxes.Title.String)
        else
            set(gui.Plot.meanOv.Children(2), 'XLim', ViewAxes.XLim)
            set(gui.Plot.meanOv.Children(2).Title, 'String', ViewAxes.Title.String)
        end
        close(temp);
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class        
end