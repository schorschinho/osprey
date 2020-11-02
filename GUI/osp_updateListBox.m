function osp_updateListBox(gui)
%% osp_updateListBox
%   This function updates the list box.
%
%
%   USAGE:
%       osp_updateListBox(gui);
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
%%% 1. UPDATE TABS %%%

    gui.controls.Selected = get(gui.layout.ListBox, 'value');
    switch gui.layout.tabs.Selection %Which tab?
        case 1 %Load tab?
            osp_updateLoadWindow(gui);
        case 2 %Process tab?
            gui.InfoText.pro = gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children;
            % Grid for Plot and Data control sliders
            gui.Plot.pro = gui.layout.(gui.layout.proTabhandles{gui.process.Selected});
            osp_updateProWindow(gui);
        case 3 %Fit tab?
            osp_updateFitWindow(gui);
        case 4 %Coreg tab?
            osp_updateCoregWindow(gui);
        case 5 %Quantify tab?
            osp_updateQuantifyWindow(gui);
    end
end