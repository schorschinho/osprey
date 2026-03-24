function osp_pop_whichcorrOvCorr_Call(src,~,gui)
%% osp_pop_whichcorrOvCorr_Call
%   This function is triggered when the correlation overview choice popup menu is changed. It refreshes
%   the figure.
%
%
%   USAGE:
%       osp_pop_whichcorrOvCorr_Call(src,~,gui);
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
        gui.overview.Selected.CorrChoice = idx;
        gui.overview.Selected.Corr = 1;
        Selection = gui.quant.popMenuNames{gui.quant.Selected.Quant};
        split_Selection = strsplit(Selection,'-');  
        if idx == 1
            set(gui.controls.pop_corrOvCorr, 'String', gui.overview.Names.QM);
             set(gui.controls.pop_corrOvCorr, 'Value', gui.overview.Selected.Corr);

        else if idx == 2
            set(gui.controls.pop_corrOvCorr, 'String', MRSCont.quantify.metabs.(split_Selection{1}));
            set(gui.controls.pop_corrOvCorr, 'Value', gui.overview.Selected.Corr);
            else
                set(gui.controls.pop_corrOvCorr, 'String', gui.overview.Names.Corr);
                set(gui.controls.pop_corrOvCorr, 'Value', gui.overview.Selected.Corr);
            end
        end
        osp_updatecorrOvWindow(gui);
end % pop_corrOvCorr_Call