function osp_SelectionChangedFcn(src,~,gui)
%% osp_SelectionChangedFcn
%   This function is triggered when the main tab is changed. It refreshes
%   the figure.
%
%
%   USAGE:
%       osp_SelectionChangedFcn(src,~,gui);
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
%%% 1. GET HANDLES %%%
% Get MRSCont from hidden container in gui class
MRSCont = getappdata(gui.figure,'MRSCont');     
% User selected tab refreshs plot
   OldValue = src.Selection;
   switch OldValue
       case 1 %Load tab
            spec = gui.load.Selected;
            if ~(spec == 2 || spec ==3)
                gui.fit.Selected = 1;
            else
                if spec == 2
                    gui.fit.Selected = 2;
                else
                    gui.fit.Selected = 3;
                end
            end
       case 2 %Process tab
            spec = gui.process.Selected;

            if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                if ~(spec == 5 || spec == 6)
                    gui.load.Selected = 1;
                    gui.fit.Selected = 1;
                else
                    if spec == 5
                        gui.load.Selected = 2;
                        gui.fit.Selected = 4;
                    else
                        gui.load.Selected = 3;
                        gui.fit.Selected = 5;
                    end
                end
            elseif MRSCont.flags.isMEGA
                if ~(spec == 3 || spec == 4)
                    gui.load.Selected = 1;
                    gui.fit.Selected = 1;
                else
                    if spec == 3
                        gui.load.Selected = 2;
                        gui.fit.Selected = 3;
                    else
                        gui.load.Selected = 3;
                        gui.fit.Selected = 4;
                    end
                end
            elseif MRSCont.flags.isUnEdited
                        gui.load.Selected = spec;
                        gui.fit.Selected = spec;
            end

       case 4 %Fit Tab
            if (strcmp(gui.fit.Names{gui.fit.Selected},'off') || strcmp(gui.fit.Names{gui.fit.Selected},'diff1')||...
                strcmp(gui.fit.Names{gui.fit.Selected},'diff2') || strcmp(gui.fit.Names{gui.fit.Selected},'sum'))
                spec = 1;
            else
                if strcmp(gui.fit.Names{gui.fit.Selected},'w')
                    spec = 3;
                else
                    spec = 2;
                end
            end
            if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                if ~(spec == 2 || spec == 3)
                    gui.load.Selected = 1;
                    gui.process.Selected = 1;
                else
                    if spec == 2
                        gui.load.Selected = 2;
                        gui.process.Selected = 5;
                    else
                        gui.load.Selected = 3;
                        gui.process.Selected = 6;
                    end
                end
            elseif MRSCont.flags.isMEGA
                if ~(spec == 2 || spec == 3)
                    gui.load.Selected = 1;
                    gui.process.Selected = 1;
                else
                    if spec == 2
                        gui.load.Selected = 2;
                        gui.process.Selected = 3;
                    else
                        gui.load.Selected = 3;
                        gui.process.Selected = 4;
                    end
                end
            elseif MRSCont.flags.isUnEdited
                        gui.load.Selected = spec;
                        gui.process.Selected = spec;
            end
       case 5 %Quantify tab?
       otherwise
           idx = 1;
           spec = 1;
   end

%%% 2. UPDATE GUI %%%
   switch gui.layout.tabs.Selection
        case 1 %Load tab
        gui.layout.ListBox.Enable = 'on';
        gui.InfoText.data = gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children;
        gui.Plot.data = gui.layout.(gui.layout.rawTabhandles{gui.load.Selected});
        set(gui.layout.rawTab, 'selection', gui.load.Selected);
        osp_updateLoadWindow(gui);
        case 2 %Process tab
        gui.layout.ListBox.Enable = 'on';
        gui.InfoText.pro = gui.layout.(gui.layout.proTabhandles{gui.load.Selected}).Children(2).Children;
        gui.Plot.pro = gui.layout.(gui.layout.proTabhandles{gui.process.Selected});
        set(gui.layout.proTab, 'selection', gui.process.Selected);
        osp_updateProWindow(gui);
        case 3 %Fit tab
        gui.layout.ListBox.Enable = 'on';
        gui.InfoText.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected}).Children(2).Children;
        gui.Plot.fit = gui.layout.(gui.layout.fitTabhandles{gui.fit.Selected});
        set(gui.layout.fitTab, 'selection', gui.fit.Selected);
        osp_updateFitWindow(gui);
        case 4 %Coreg Tab
        gui.layout.ListBox.Enable = 'on';
        gui.InfoText.coreg = gui.layout.(gui.layout.proTabhandles{gui.load.Selected}).Children(2).Children;
        osp_updateCoregWindow(gui);
        case 5 % Quantify tab
        gui.layout.ListBox.Enable = 'on';
        osp_updateQuantifyWindow(gui);
       case 6 %Overview tab
        gui.layout.ListBox.Enable = 'off';
    end
end