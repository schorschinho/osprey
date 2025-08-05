function osp_pop_multiverseOv_Call(src,~,gui) 
%% osp_pop_multiverseOv_Call
%   This function is triggered when either of the multiverse popup menus is changed. It refreshes
%   the figure.
%
%
%   USAGE:
%       osp_pop_multiverseOv_Call(src,~,gui);
%
%   INPUT:  src      = handle of the fit tabs
%           gui      = gui class containing all handles and the MRSCont             
%
%
%   AUTHORS:
%       Chris Davies-Jenkins (Johns Hopkins University, 2025-07-15)
%       cdavies9@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2025-07-15: First version of the code.

MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
idx=(src.Value);
src.Value=idx;



% If Plot type is changed, reset other menus
if iscell(src.String) && length(src.String{1})>4 && matches(src.String{1}(1:5),'Inter') %If PlotType updated, then update secondary menus:
    gui.controls.pop_multiverseOvPlotArg.Value = 1;
end

% Ascertain plot type from first drop-down menu
PlotType = gui.controls.pop_multiverseOvPlotType.String{gui.controls.pop_multiverseOvPlotType.Value};

% Handle specific plot types—enabling/disabling the appropriate menus:
switch PlotType(1:5)
    case 'Speci'
        SubSpecStr = PlotType(22:end-1);
        SubSpec = matches(MRSCont.overview.FitSpecNamesStruct.metab, SubSpecStr);
        gui.controls.pop_multiverseOvPlotArg2.Enable = "on";
        gui.controls.pop_multiverseOvPlotArg.String = unique([MRSCont.quantify.names.metab{1,SubSpec,:,:}])';
        gui.controls.pop_multiverseOvPlotArg.Enable = 'on';
    case 'Inter'
        gui.controls.pop_multiverseOvPlotArg.String = num2str((1:MRSCont.nDatasets(1))');
        gui.controls.pop_multiverseOvPlotArg2.Enable = "on";
        if gui.controls.check_multiverseOvPlotMed.Value
            gui.controls.pop_multiverseOvPlotArg.Enable = 'off';
        else
            gui.controls.pop_multiverseOvPlotArg.Enable = 'on';
        end
    case 'Model'
        gui.controls.pop_multiverseOvPlotArg.String = num2str((1:MRSCont.nDatasets(1))');
        gui.controls.pop_multiverseOvPlotArg2.Enable = "off";
        if gui.controls.check_multiverseOvPlotMed.Value
            gui.controls.pop_multiverseOvPlotArg.Enable = 'off';
        else
            gui.controls.pop_multiverseOvPlotArg.Enable = 'on';
        end
    otherwise
        error
end

osp_updatemultiverseOvWindow(gui); % Update Window

end