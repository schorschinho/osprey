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
        Selection = gui.quant.popMenuNames{gui.quant.Selected.Quant};
        if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            set(gui.controls.distrOvPanel,'Title',['Actual Quantification and Metabolite in Voxel ' num2str(gui.controls.act_x)]);
            set(gui.layout.distrOvTab.Children.Children(1).Children(3).Children.Children.Children(4),'String',gui.controls.act_z)
            set(gui.layout.distrOvTab.Children.Children(1).Children(3).Children.Children.Children(5),'String',gui.controls.act_y)
            set(gui.layout.distrOvTab.Children.Children(1).Children(3).Children.Children.Children(6),'String',gui.controls.act_x)
        end
%%% 2. VISUALIZATION PART OF THIS TAB %%%        
        if ~strcmp(Selection,'Quality')    
        split_Selection = strsplit(Selection,'-');
            if strcmp(split_Selection{2},'AlphaCorrWaterScaled') || strcmp(split_Selection{2},'AlphaCorrWaterScaledGroupNormed')
                metab = 'GABA';
            else
                metab = MRSCont.quantify.metabs.(split_Selection{1}){gui.overview.Selected.Metab};
            end
            if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                if ~gui.controls.GM
                    [temp] = osp_plotRaincloud(MRSCont,split_Selection{1},split_Selection{2},metab,'Raincloud plot',0,gui.controls.act_x); 
                else
                    [temp] = osp_plotRaincloud(MRSCont,split_Selection{1},split_Selection{2},metab,'Raincloud plot',1,gui.controls.act_x);                
                end
            else
                if ~gui.controls.GM
                    [temp] = osp_plotRaincloud(MRSCont,split_Selection{1},split_Selection{2},metab,'Raincloud plot'); 
                else
                    [temp] = osp_plotRaincloud(MRSCont,split_Selection{1},split_Selection{2},metab,'Raincloud plot',1);                
                end
            end
        else
            quality = {'SNR','FWHM','freqShift'};
            if ~gui.controls.GM
                [temp] = osp_plotRaincloud(MRSCont,'Quality','Quality',quality{gui.overview.Selected.Metab},'Raincloud plot');  
            else
                [temp] = osp_plotRaincloud(MRSCont,'Quality','Quality',quality{gui.overview.Selected.Metab},'Raincloud plot',1);  
            end
            split_Selection{1}=quality{gui.overview.Selected.Metab};
            metab = '';
        end

            set( temp.Children(2).Children, 'Parent', gui.Plot.distrOv.Children(3) );
            set(  gui.Plot.distrOv.Children(3), 'XLabel', temp.Children(2).XLabel);
            set(  gui.Plot.distrOv.Children(3), 'YLim', temp.Children(2).YLim);
            close(temp);
            set(gui.Plot.distrOv,'Heights', [-0.1 -0.87 -0.03]);
            if ~gui.controls.GM
                gui.Plot.distrOv.Children(3).Legend.Location = 'North'; % Update legend
            else
                gui.Plot.distrOv.Children(3).Legend.Location = 'North'; % Update legend
            end
            set(gui.Plot.distrOv.Children(3).Title, 'String', ['Raincloud plot voxel ' num2str(gui.controls.act_x) ' : ' split_Selection{1} ' ' metab]) %Update title
        setappdata(gui.figure,'MRSCont',MRSCont); %Write  MRSCont into hidden container in gui class            
end
