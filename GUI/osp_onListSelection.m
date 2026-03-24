function osp_onListSelection( src, ~,gui)
%% osp_onListSelection
%   This function is triggered when a file is selected in the listbox. It refreshes
%   the figure.
%
%
%   USAGE:
%       osp_onListSelection(src,~,gui);
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
%%% 1. UPDATE GUI %%%
    % User selected file in the list refreshs active tab
    if gui.controls.KeyPress == 0
        gui.controls.Selected = get( src, 'Value' );
        set(gui.layout.ListBox, 'value', gui.controls.Selected);
        switch gui.layout.tabs.Selection
            case 1
                osp_updateLoadWindow(gui);
            case 2
                gui.InfoText.pro = gui.layout.(gui.layout.proTabhandles{gui.load.Selected}).Children(2).Children;
                % Grid for Plot and Data control sliders
                gui.Plot.pro = gui.layout.(gui.layout.proTabhandles{gui.process.Selected});
                osp_updateProWindow(gui);
            case 3
                osp_updateFitWindow(gui);
            case 4
                osp_updateCoregWindow(gui);
            case 5
                osp_updateQuantifyWindow(gui);
        end
    else
        gui.controls.Selected = get(gui.layout.ListBox, 'value');
    end
end