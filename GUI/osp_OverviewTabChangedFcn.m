function osp_OverviewTabChangedFcn(src,~,gui) 
%% osp_OverviewTabChangedFcn
%   This function is triggered when the overview tab is changed. It refreshes
%   the figure.
%
%
%   USAGE:
%       osp_OverviewTabChangedFcn(src,~,gui);
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
%%% 1. CALL THE RIGHT FUNCTION %%%
    NewValue= src.Selection;
    switch NewValue
       case 1
            osp_updateSpecsOvWindow(gui);
            set(gui.controls.pop_specsOvPlot, 'value',gui.process.Selected)
       case 2
            osp_updatemeanOvWindow(gui);
            set(gui.controls.pop_meanOvPlot, 'value',gui.process.Selected)
       case 3
            set(gui.layout.overviewTab, 'selection', 3);
       case 4
            osp_updatedistrOvWindow(gui);
            set(gui.controls.pop_distrOvQuant, 'value',gui.quant.Selected.Quant)
            set(gui.controls.pop_distrOvMetab, 'value',gui.overview.Selected.Metab)
       case 5
            osp_updatecorrOvWindow(gui);
            set(gui.controls.pop_corrOvQuant, 'value',gui.quant.Selected.Quant)
            set(gui.controls.pop_corrOvMetab, 'value',gui.overview.Selected.Metab)
      otherwise
            set(gui.layout.overviewTab, 'selection', 1);
    end
end