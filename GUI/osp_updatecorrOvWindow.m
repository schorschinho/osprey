function osp_updatecorrOvWindow(gui) 
%% osp_updatecorrOvWindow
%   This function updates the correlation overview tab.
%
%
%   USAGE:
%       osp_updatecorrOvWindow(gui);
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
            delete(gui.Plot.corrOv.Children(3).Children)
            Selection = gui.quant.popMenuNames{gui.quant.Selected.Quant};
            split_Selection = strsplit(Selection,'-');  
            ind = find(strcmp(MRSCont.overview.FitSpecNamesStruct.(split_Selection{1})(1,:),split_Selection{2}));
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
           if strcmp(split_Selection{3},'AlphaCorrWaterScaled') || strcmp(split_Selection{3},'AlphaCorrWaterScaledGroupNormed')
                metab = 'GABA';
            else
                metab = MRSCont.quantify.names.(split_Selection{1}){gui.controls.act_z,ind}{gui.overview.Selected.Metab};
           end 
           all_subs_ind = 1:MRSCont.nDatasets(1);
           if isfield(MRSCont,'exclude')
                if~isempty(MRSCont.exclude)
                    all_subs_ind(MRSCont.exclude)=[]; 
                end
            end
           if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                if gui.overview.Selected.CorrChoice == 1
                    switch gui.overview.Selected.Corr
                        case 1
                        temp = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,MRSCont.QM.SNR.metab(1,all_subs_ind,ind)',gui.overview.Names.QM{gui.overview.Selected.Corr},1,gui.controls.act_z);
                        case 2
                        temp = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,MRSCont.QM.FWHM.metab(1,all_subs_ind,ind)',gui.overview.Names.QM{gui.overview.Selected.Corr},1,gui.controls.act_z);
                    end                                
                else if gui.overview.Selected.CorrChoice == 2
                    temp = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,MRSCont.quantify.names.(split_Selection{1}){1,ind}{gui.overview.Selected.Corr},MRSCont.quantify.names.(split_Selection{1}){1,ind}{gui.overview.Selected.Corr},1,gui.controls.act_z);
                    else
                        temp = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,gui.overview.CorrMeas{gui.overview.Selected.Corr},gui.overview.Names.Corr{gui.overview.Selected.Corr},1,gui.controls.act_z);
                    end
                end
           else
               if gui.overview.Selected.CorrChoice == 1
                    switch gui.overview.Selected.Corr
                        case 1
                        temp = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,MRSCont.QM.SNR.metab(1,all_subs_ind,ind)',gui.overview.Names.QM{gui.overview.Selected.Corr},gui.controls.act_x,gui.controls.act_z);
                        case 2
                        temp = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,MRSCont.QM.FWHM.metab(1,all_subs_ind,ind)',gui.overview.Names.QM{gui.overview.Selected.Corr},gui.controls.act_x,gui.controls.act_z);
                    end                                
                else if gui.overview.Selected.CorrChoice == 2
                    temp = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,MRSCont.quantify.names.(split_Selection{1}){1,ind}{gui.overview.Selected.Corr},MRSCont.quantify.names.(split_Selection{1}){1,ind}{gui.overview.Selected.Corr},gui.controls.act_x,gui.controls.act_z);
                    else
                        temp = osp_plotScatter(MRSCont,split_Selection{2},split_Selection{3},metab,gui.overview.CorrMeas{gui.overview.Selected.Corr},gui.overview.Names.Corr{gui.overview.Selected.Corr},gui.controls.act_x,gui.controls.act_z);
                    end
                end              
           end
            set(gui.Plot.corrOv.Children(3), 'XLim', temp.Children(2).XLim);
            set(gui.Plot.corrOv.Children(3), 'YLim', temp.Children(2).YLim);
            set(temp.Children(2).Children, 'Parent', gui.Plot.corrOv.Children(3) );
            set(gui.Plot.corrOv.Children(3), 'XLabel', temp.Children(2).XLabel);
            set(gui.Plot.corrOv.Children(3), 'YLabel', temp.Children(2).YLabel);
            close(temp);
            set(gui.Plot.corrOv,'Heights', [-0.1 -0.87 -0.03]);
            gui.Plot.corrOv.Children(3).Legend.Location = 'North'; %Update legend
            if gui.overview.Selected.CorrChoice == 1
                switch gui.overview.Selected.Corr
                    case 1
                    set(gui.Plot.corrOv.Children(3).Title, 'String', ['Voxel ' num2str(gui.controls.act_x) ' ' metab ' vs SNR']) %Update title
                    case 2
                    set(gui.Plot.corrOv.Children(3).Title, 'String', ['Voxel ' num2str(gui.controls.act_x) ' ' metab ' vs FHWM (ppm)']) %Update title
                end                                
            else if gui.overview.Selected.CorrChoice == 2
                set(gui.Plot.corrOv.Children(3).Title, 'String', ['Voxel ' num2str(gui.controls.act_x) ' ' MRSCont.quantify.names.(split_Selection{1}){gui.controls.act_z,ind}{gui.overview.Selected.Corr} ' vs ' metab]) %Update title
                else
                    set(gui.Plot.corrOv.Children(3).Title, 'String', ['Voxel ' num2str(gui.controls.act_x) ' ' gui.overview.Names.Corr{gui.overview.Selected.Corr} ' vs ' metab]) %Update title
                end
            end
        setappdata(gui.figure,'MRSCont',MRSCont); %Write  MRSCont into hidden container in gui class            
end