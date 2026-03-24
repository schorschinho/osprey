function osp_check_distrOv_Call(src,~,gui)
%% osp_check_distrOv_Call
%   This function is triggered when Grand Mean is checked for the raincloud
%   plots
%
%
%   USAGE:
%       osp_check_distrOv_Call(src,~,gui);
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
        gui.controls.GM = idx;
        metab = get(gui.controls.pop_distrOvMetab, 'Value');
        Selection = gui.quant.popMenuNames{gui.quant.Selected.Quant};
        if ~strcmp(Selection,'Quality')    
            split_Selection = strsplit(Selection,'-');        
            if strcmp(split_Selection{2},'AlphaCorrWaterScaled') || strcmp(split_Selection{2},'AlphaCorrWaterScaledGroupNormed')
               set(gui.controls.pop_distrOvMetab, 'String', {'GABA'});
               set(gui.controls.pop_distrOvMetab, 'Value', gui.quant.idx.GABA);
               set(gui.controls.pop_distrOvMetab, 'Enable', 'off');
            else
               set(gui.controls.pop_distrOvMetab, 'String', MRSCont.quantify.metabs.(split_Selection{1}));
               set(gui.controls.pop_distrOvMetab, 'Value', metab);
               set(gui.controls.pop_distrOvMetab, 'Enable', 'on');
            end
        else
               gui.overview.Selected.Metab = 1;
               set(gui.controls.pop_distrOvMetab, 'String', {'SNR','FWHM','freqShift'});
               set(gui.controls.pop_distrOvMetab, 'Value', 1);
               set(gui.controls.pop_distrOvMetab, 'Enable', 'on');            
        end
        osp_updatedistrOvWindow(gui);
end % pop_distrOvQuant_Call