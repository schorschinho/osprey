function osp_pop_corrOvQuant_Call(src,~,gui) 
%% osp_pop_corrOvQuant_Call
%   This function is triggered when the correlation quantification popup menu is changed. It refreshes
%   the figure.
%
%
%   USAGE:
%       osp_pop_corrOvQuant_Call(src,~,gui);
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
        MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
        idx=(src.Value);
        src.Value=idx;
        gui.quant.Selected.Quant = idx;
        Selection = gui.quant.popMenuNames{gui.quant.Selected.Quant};
        selectedOvTab = get(gui.layout.overviewTab,'Selection');         
        set(gui.layout.(gui.layout.overviewTabhandels{selectedOvTab}).Children(1).Children(1).Children(3).Children(1).Children(1).Children(3),'String',gui.controls.act_z)
        set(gui.layout.(gui.layout.overviewTabhandels{selectedOvTab}).Children(1).Children(1).Children(3).Children(1).Children(1).Children(4),'String',gui.controls.act_x)             
        split_Selection = strsplit(Selection,'-');
        ind = find(strcmp(MRSCont.overview.FitSpecNamesStruct.(split_Selection{1})(1,:),split_Selection{2}));
        if strcmp(split_Selection{3},'AlphaCorrWaterScaled') || strcmp(split_Selection{3},'AlphaCorrWaterScaledGroupNormed')
               set(gui.controls.pop_corrOvMetab, 'String', MRSCont.quantify.names.(split_Selection{1}){gui.controls.act_z,ind});
               set(gui.controls.pop_corrOvMetab, 'String', {'GABA'});
               set(gui.controls.pop_corrOvMetab, 'Value', gui.quant.idx.GABA);
               set(gui.controls.pop_corrOvMetab, 'Enable', 'on');
            else
                metab_idx_old = get(gui.controls.pop_corrOvMetab, 'Value');
                name_list=get(gui.controls.pop_corrOvMetab, 'String');
                name_old = name_list{metab_idx_old}; 
                name_list_new=MRSCont.quantify.names.(split_Selection{1}){gui.controls.act_z,ind};
                metab_idx_new = find(strcmp(name_list_new,name_old));   
                if ~isempty(metab_idx_new)
                    set(gui.controls.pop_corrOvMetab, 'String', MRSCont.quantify.names.(split_Selection{1}){gui.controls.act_z,ind});
                else
                    metab_idx_new = 1;
                    set(gui.controls.pop_corrOvMetab, 'Value', 1);
                    set(gui.controls.pop_corrOvMetab, 'String', MRSCont.quantify.names.(split_Selection{1}){gui.controls.act_z,ind});
                    gui.overview.Selected.Metab = 1;
                end
               set(gui.controls.pop_corrOvMetab, 'Value', metab_idx_new);
               set(gui.controls.pop_corrOvMetab, 'Enable', 'on');                             
        end
        setappdata(gui.figure,'MRSCont',MRSCont); %Write  MRSCont into hidden container in gui class       
        osp_updatecorrOvWindow(gui);
end % pop_corrOvQuant_Call