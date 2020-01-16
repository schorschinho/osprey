function osp_updatedistrOvWindow(gui) 
%% osp_updatedistrOvWindow
%   This function updates the raincloud overview tab.
%
%
%   USAGE:
%       osp_updateCoregWindow(gui);
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
            delete(gui.Plot.distrOv.Children(3).Children)
%%% 2. VISUALIZATION PART OF THIS TAB %%%
            temp = figure( 'Visible', 'off' );
            if strcmp(gui.quant.Names.Quants(gui.quant.Selected.Quant),'AlphaCorrWaterScaled') || strcmp(gui.quant.Names.Quants(gui.quant.Selected.Quant),'AlphaCorrWaterScaledGroupNormed')
                metab = 'GABA';
            else
                metab = MRSCont.quantify.metabs{gui.overview.Selected.Metab};
            end
            [temp] = osp_plotRaincloud(MRSCont,gui.quant.Names.Model{gui.quant.Selected.Model},gui.quant.Names.Quants{gui.quant.Selected.Quant},metab,'Raincloud plot',1);
            set( temp.Children(2).Children, 'Parent', gui.Plot.distrOv.Children(3) );
            set(  gui.Plot.distrOv.Children(3), 'XLabel', temp.Children(2).XLabel);
            set(  gui.Plot.distrOv.Children(3), 'YLim', temp.Children(2).YLim);
            close(temp);
            set(gui.Plot.distrOv,'Heights', [-0.07 -0.90 -0.03]);
            gui.Plot.distrOv.Children(3).Legend.Location = 'North'; % Update legend
            set(gui.Plot.distrOv.Children(3).Title, 'String', ['Raincloud plot: ' metab]) %Update title
        setappdata(gui.figure,'MRSCont',MRSCont); %Write  MRSCont into hidden container in gui class            
end
