function osp_updatemultiverseOvWindow(gui) 
%% osp_updatemultiverseOvWindow
%   This function updates the multiverse overview tab.
%
%
%   USAGE:
%       osp_updatemultiverseOvWindow(gui);
%
%   INPUT:  
%           gui      = gui class containing all handles and the MRSCont             
%
%
%   AUTHORS:
%       Christopher Davies-Jenkins (Johns Hopkins University, 2025-07-08)
%       cdavies9@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2025-07-08: First version of the code.
%%% 1. INITIALIZE %%%
MRSCont = getappdata(gui.figure,'MRSCont');  % Get MRSCont from hidden container in gui class

if length(gui.Plot.multiverse.Children)>1
    delete(gui.Plot.multiverse.Children(1))
end

PlotType = gui.controls.pop_multiverseOvPlotType.String{gui.controls.pop_multiverseOvPlotType.Value};
Checkbox = gui.controls.check_multiverseOvPlotMed.Value;
Reference = gui.controls.pop_multiverseOvPlotArg2.String{gui.controls.pop_multiverseOvPlotArg2.Value};

%%%%% Hardcoded for now
Exp = 1;
%%%%%%%%%%%%%%%%%%%%%%%


switch PlotType(1:5)
    case "Speci" % Specification curve
        
        Metab    = gui.controls.pop_multiverseOvPlotArg.String{gui.controls.pop_multiverseOvPlotArg.Value};
        if Checkbox
            NDataset = 0;
        else
            warning("Not yet implemented...")
            NDataset = 0;
        end
        
        % ascertain sub spectrum from label
        SubSpecStr = PlotType(22:end-1);
        SubSpec = find(matches(MRSCont.overview.FitSpecNamesStruct.metab, SubSpecStr));
        
        % Embed the 2-panel Specification curve in a HBox:
        gui.Plot.multiverseHBx = uix.HBox('Parent', gui.Plot.multiverse, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
        gui.Plot.multiverseGap = uicontrol('Parent',gui.Plot.multiverseHBx,'Style','text','String','  ',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,'HorizontalAlignment','center');
        gui.Plot.multiversePlts = uix.VBox('Parent', gui.Plot.multiverseHBx, 'Padding', 45,'Spacing', 5,'BackgroundColor',gui.colormap.Background);
        
        % Plot the specification curve
        temp = osp_plotSpecificationcurve(MRSCont,Metab,Reference,NDataset,SubSpec,Exp);
        set(temp.Children(2), 'Parent', gui.Plot.multiversePlts);
        set(temp.Children(1), 'Parent', gui.Plot.multiversePlts);
        close(temp);

        % Formatting
        set(gui.Plot.multiverse,'Heights', [-0.1 -0.9]);
        set(gui.Plot.multiverseHBx, 'Widths', [-0.03 -0.97])
        set(gui.Plot.multiversePlts,'Heights', [-0.2 -0.7]);
    case "Inter" % Inter-model distributions
        if Checkbox
            NDataset = 0;
        else
            NDataset = str2num(gui.controls.pop_multiverseOvPlotArg.String(gui.controls.pop_multiverseOvPlotArg.Value));
        end

        % ascertain sub spectrum from label
        SubSpecStr = PlotType(21:end-1);
        SubSpec = find(matches(MRSCont.overview.FitSpecNamesStruct.metab, SubSpecStr));
        
        % Embed plot in a HBox
        gui.Plot.multiverseHBx = uix.HBox('Parent', gui.Plot.multiverse, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
        gui.Plot.multiverseGap = uicontrol('Parent',gui.Plot.multiverseHBx,'Style','text','String','  ',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,'HorizontalAlignment','center');
        gui.Plot.multiversePlts = uix.VBox('Parent', gui.Plot.multiverseHBx, 'Padding', 45,'Spacing', 5,'BackgroundColor',gui.colormap.Background);

        % plot distributions
        temp = osp_plotMultiverse(MRSCont,Reference,NDataset,SubSpec,Exp);
        set(temp.Children(1), 'Parent', gui.Plot.multiversePlts);
        
        % Formatting
        set(gui.Plot.multiverse,'Heights', [-0.1 -0.9]);
        set(gui.Plot.multiverseHBx, 'Widths', [-0.03 -0.97])
    case "Model" % Model hotspot plots
        if Checkbox
            NDataset = 0;
        else
            NDataset = str2num(gui.controls.pop_multiverseOvPlotArg.String(gui.controls.pop_multiverseOvPlotArg.Value));
        end

        % ascertain sub spectrum from label
        SubSpecStr = PlotType(17:end-1);
        SubSpec = find(matches(MRSCont.overview.FitSpecNamesStruct.metab, SubSpecStr));
        
        % Embed plot in HBox
        gui.Plot.multiverseHBx = uix.HBox('Parent', gui.Plot.multiverse, 'Padding', 5,'BackgroundColor',gui.colormap.Background);
        gui.Plot.multiverseGap = uicontrol('Parent',gui.Plot.multiverseHBx,'Style','text','String','  ',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,'HorizontalAlignment','center');
        gui.Plot.multiversePlts = uix.VBox('Parent', gui.Plot.multiverseHBx, 'Padding', 45,'Spacing', 5,'BackgroundColor',gui.colormap.Background);
        
        % Plot hotspot figure
        temp = osp_plotModelvariation(MRSCont,NDataset,SubSpec,Exp);
        set(temp.Children(1), 'Parent', gui.Plot.multiversePlts);
        
        % Formatting
        set(gui.Plot.multiverse,'Heights', [-0.1 -0.9]);
        set(gui.Plot.multiverseHBx, 'Widths', [-0.03 -0.97])
end
end
