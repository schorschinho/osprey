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
        selectedOvTab = get(gui.layout.overviewTab,'Selection');
        if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            set(gui.controls.distrOvPanel,'Title',['Actual Quantification and Metabolite in Voxel ' num2str(gui.controls.act_x)]);
            set(gui.layout.(gui.layout.overviewTabhandels{selectedOvTab}).Children.Children(1).Children(3).Children.Children.Children(4),'String',gui.controls.act_z)
            set(gui.layout.(gui.layout.overviewTabhandels{selectedOvTab}).Children.Children(1).Children(3).Children.Children.Children(5),'String',gui.controls.act_y)
            set(gui.layout.(gui.layout.overviewTabhandels{selectedOvTab}).Children.Children(1).Children(3).Children.Children.Children(6),'String',gui.controls.act_x)
        else            
            set(gui.layout.(gui.layout.overviewTabhandels{selectedOvTab}).Children(1).Children(1).Children(3).Children(1).Children(1).Children(3),'String',gui.controls.act_z)
            set(gui.layout.(gui.layout.overviewTabhandels{selectedOvTab}).Children(1).Children(1).Children(3).Children(1).Children(1).Children(4),'String',gui.controls.act_x)  
        end
%%% 2. VISUALIZATION PART OF THIS TAB %%%        
        if ~strcmp(Selection,'Quality')    
        split_Selection = strsplit(Selection,'-');
            ind = find(strcmp(MRSCont.overview.FitSpecNamesStruct.(split_Selection{1})(1,:),split_Selection{2}));   
            ind = ind(1);
            if strcmp(split_Selection{3},'AlphaCorrWaterScaled') || strcmp(split_Selection{3},'AlphaCorrWaterScaledGroupNormed')
                metab = 'GABA';
            else
                metab = MRSCont.quantify.names.(split_Selection{1}){gui.controls.act_z,ind}{gui.overview.Selected.Metab};
            end
            if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                if ~gui.controls.GM
                    [temp] = osp_plotRaincloud(MRSCont,split_Selection{2},split_Selection{3},metab,'Raincloud plot',0,gui.controls.act_x,gui.controls.act_z); 
                else
                    [temp] = osp_plotRaincloud(MRSCont,split_Selection{2},split_Selection{3},metab,'Raincloud plot',1,gui.controls.act_x,gui.controls.act_z);                
                end
            else
                if ~gui.controls.GM
                    [temp] = osp_plotRaincloud(MRSCont,split_Selection{2},split_Selection{3},metab,'Raincloud plot',0,1,gui.controls.act_z,gui.controls.act_x); 
                else
                    [temp] = osp_plotRaincloud(MRSCont,split_Selection{2},split_Selection{3},metab,'Raincloud plot',1,1,gui.controls.act_z,gui.controls.act_x);                
                end
            end
            set(gui.controls.pop_distrOvMetab, 'String', MRSCont.quantify.names.(split_Selection{1}){gui.controls.act_z,ind});
            split_Selection{4}=['basis ' num2str(gui.controls.act_z)];
        else
            quality = {'SNR','FWHM','freqShift'};
            if ~gui.controls.GM
                [temp] = osp_plotRaincloud(MRSCont,'Quality','Quality',quality{gui.overview.Selected.Metab},'Raincloud plot');  
            else
                [temp] = osp_plotRaincloud(MRSCont,'Quality','Quality',quality{gui.overview.Selected.Metab},'Raincloud plot',1);  
            end
            split_Selection{2}='Quality';
            split_Selection{1}='Spectral';
            split_Selection{3}=quality{gui.overview.Selected.Metab};
            split_Selection{4}='';
            metab = '';
            set(gui.controls.pop_distrOvMetab, 'String', quality);
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
            set(gui.Plot.distrOv.Children(3).Title, 'String', ['Raincloud plot voxel ' num2str(gui.controls.act_x) ' : ' split_Selection{1} ' ' split_Selection{2} ' ' split_Selection{3}  ' ' split_Selection{4} ' ' metab]) %Update title
        setappdata(gui.figure,'MRSCont',MRSCont); %Write  MRSCont into hidden container in gui class            
end
