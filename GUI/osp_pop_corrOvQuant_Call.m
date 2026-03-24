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
        split_Selection = strsplit(Selection,'-');
        if strcmp(split_Selection{2},'AlphaCorrWaterScaled') || strcmp(split_Selection{2},'AlphaCorrWaterScaledGroupNormed')
           set(gui.controls.pop_corrOvMetab, 'String', {'GABA'});
           set(gui.controls.pop_corrOvMetab, 'Value', gui.quant.idx.GABA);
           set(gui.controls.pop_corrOvMetab, 'Enable', 'off');
           set(gui.overview.Selected, 'Metab',1)
        else
           set(gui.controls.pop_corrOvMetab, 'String', MRSCont.quantify.metabs.(split_Selection{1}));
           set(gui.controls.pop_corrOvMetab, 'Value', gui.overview.Selected.Metab);
           set(gui.controls.pop_corrOvMetab, 'Enable', 'on');
        end
        setappdata(gui.figure,'MRSCont',MRSCont); %Write  MRSCont into hidden container in gui class       
        osp_updatecorrOvWindow(gui);
end % pop_corrOvQuant_Call