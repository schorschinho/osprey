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
   OldValue = gui.controls.Tab;
   switch OldValue
       case 1 %Load tab
            tab = gui.layout.rawTab.TabTitles{gui.load.Selected};
            switch tab
                case {'metabolites','MM'}
                    gui.process.Selected = gui.load.Selected;
                    gui.fit.Selected = gui.load.Selected;
                case {'reference'}
                    gui.process.Selected = find(strcmp(gui.layout.proTab.TabTitles,'ref'));
                    gui.process.Selected = find(strcmp(gui.layout.fitTab.TabTitles,'ref'));
                case {'water'}
                    gui.process.Selected = find(strcmp(gui.layout.proTab.TabTitles,'w'));
                    gui.process.Selected = find(strcmp(gui.layout.fitTab.TabTitles,'w'));
                case {'MM reference'}   
                    gui.process.Selected = find(strcmp(gui.layout.proTab.TabTitles,'mm_ref'));
            end
       case 2 %Process tab
            tab = gui.layout.proTab.TabTitles{gui.process.Selected};
            switch tab
                case {'metab','mm'}
                    gui.load.Selected = gui.load.Selected;
                    gui.fit.Selected = gui.load.Selected;
                case {'ref'}
                    gui.load.Selected = find(strcmp(gui.layout.rawTab.TabTitles,'reference'));
                    gui.process.Selected = find(strcmp(gui.layout.fitTab.TabTitles,'ref'));
                case {'w'}
                    gui.load.Selected = find(strcmp(gui.layout.rawTab.TabTitles,'water'));
                    gui.process.Selected = find(strcmp(gui.layout.fitTab.TabTitles,'w'));
                case {'mm_ref'}   
                    gui.load.Selected = find(strcmp(gui.layout.rawTab.TabTitles,'MM reference'));
            end
            

       case 3 %Fit Tab
            tab = gui.layout.fitTab.TabTitles{gui.fit.Selected};
            switch tab
                case {'metab','mm'}
                    gui.load.Selected = gui.load.Selected;
                    gui.fit.Selected = gui.load.Selected;
                case {'ref'}
                    gui.load.Selected = find(strcmp(gui.layout.rawTab.TabTitles,'reference'));
                    gui.process.Selected = find(strcmp(gui.layout.fitTab.TabTitles,'ref'));
                case {'w'}
                    gui.load.Selected = find(strcmp(gui.layout.rawTab.TabTitles,'water'));
                    gui.process.Selected = find(strcmp(gui.layout.fitTab.TabTitles,'w'));
            end

       otherwise
           idx = 1;
           spec = 1;
   end

   Exp = 1;
   SubSpec = 1;
   Basis = 1;
   gui.controls.Tab = gui.layout.tabs.Selection;
%%% 2. UPDATE GUI %%%
   switch gui.layout.tabs.Selection
        case 1 %Load tab
            gui.layout.ListBox.Enable = 'on';
            set(gui.layout.rawTab, 'Selection', gui.load.Selected);            
            gui.info.nXvoxels = MRSCont.nDatasets(2);
            gui.controls.act_x = Exp;
            osp_updateLoadWindow(gui);
        case 2 %Process tab
            gui.layout.ListBox.Enable = 'on';
            gui.InfoText.pro = gui.upperBox.pro.Info{gui.process.Selected}.Children;
            gui.Plot.pro = gui.layout.(gui.layout.proTabhandles{gui.process.Selected});
            set(gui.layout.proTab, 'Selection', gui.process.Selected);
            if MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.dims.extras > 0
                gui.info.nXvoxels = MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.sz(MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.dims.extras);
            else
                gui.info.nXvoxels = 1;
            end            
            gui.controls.act_x = Exp;
            if MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.dims.subSpecs > 0
                gui.info.nYvoxels = MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.sz(MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.dims.subSpecs);
            else
                gui.info.nYvoxels = MRSCont.nDatasets(2);
            end
            gui.controls.act_y = SubSpec;
            osp_updateProWindow(gui);
        case 3 %Fit tab
            gui.layout.ListBox.Enable = 'on';
            set(gui.layout.fitTab, 'Selection', gui.fit.Selected);                                                        
            gui.info.nXvoxels = 1;
            gui.controls.act_x = Exp;
            gui.controls.b_left_y.Enable = 'off';
            gui.controls.b_left_z.Enable = 'off';
            gui.controls.b_right_y.Enable = 'off';
            gui.controls.b_right_z.Enable = 'off';
            if ~strcmp(MRSCont.opts.fit.method, 'Osprey_gLCM')
                gui.info.nYvoxels = size(MRSCont.fit.results.(gui.fit.Names{gui.fit.Selected}).fitParams,3);
                gui.info.nZvoxels = size(MRSCont.fit.results.(gui.fit.Names{gui.fit.Selected}).fitParams,1);
                if size(MRSCont.fit.results.(gui.fit.Names{gui.fit.Selected}).fitParams,3) > 1
                gui.controls.b_left_y.Enable = 'on';
                gui.controls.b_right_y.Enable = 'on';
                end
                if size(MRSCont.fit.results.(gui.fit.Names{gui.fit.Selected}).fitParams,1) > 1
                    gui.controls.b_left_z.Enable = 'on';
                    gui.controls.b_right_z.Enable = 'on';
                end
            else
                gui.info.nYvoxels = size(MRSCont.fit.results.metab,3);
                gui.info.nZvoxels = size(MRSCont.fit.results.metab,1);
                if size(MRSCont.fit.results.metab,3) > 1
                    gui.controls.b_left_y.Enable = 'on';
                    gui.controls.b_right_y.Enable = 'on';
                end
                if size(MRSCont.fit.results.metab,1) > 1
                    gui.controls.b_left_z.Enable = 'on';
                    gui.controls.b_right_z.Enable = 'on';
                end
            end
            gui.controls.act_y = 1;
            gui.controls.act_z = 1;  
            osp_updateFitWindow(gui);
        case 4 %Coreg Tab
            gui.layout.ListBox.Enable = 'on';
            gui.InfoText.coreg = gui.layout.(gui.layout.proTabhandles{gui.load.Selected}).Children(2).Children;
            osp_updateCoregWindow(gui);
        case 5 % Quantify tab
            gui.layout.ListBox.Enable = 'on';           
            gui.controls.act_x = Exp;
            gui.controls.b_left_y.Enable = 'off';
            gui.controls.b_left_z.Enable = 'off';
            gui.controls.b_right_y.Enable = 'off';
            gui.controls.b_right_z.Enable = 'off';
            if ~strcmp(MRSCont.opts.fit.method, 'Osprey_gLCM')
                gui.info.nYvoxels = size(MRSCont.fit.results.(gui.fit.Names{gui.fit.Selected}).fitParams,3);
                gui.info.nZvoxels = size(MRSCont.fit.results.(gui.fit.Names{gui.fit.Selected}).fitParams,1);
                if size(MRSCont.fit.results.metab.fitParams,3) > 1
                    gui.controls.b_left_y.Enable = 'on';
                    gui.controls.b_right_y.Enable = 'on';
                end
                if size(MRSCont.fit.results.metab.fitParams,1) > 1
                    gui.controls.b_left_z.Enable = 'on';
                    gui.controls.b_right_z.Enable = 'on';
                end
            else
                gui.info.nYvoxels = size(MRSCont.fit.results.metab,3);
                gui.info.nZvoxels = size(MRSCont.fit.results.metab,1);
                if size(MRSCont.fit.results.metab,3) > 1
                    gui.controls.b_left_y.Enable = 'on';
                    gui.controls.b_right_y.Enable = 'on';
                end
                if size(MRSCont.fit.results.metab,1) > 1
                    gui.controls.b_left_z.Enable = 'on';
                    gui.controls.b_right_z.Enable = 'on';
                end
            end
            gui.controls.act_y = 1;
            gui.controls.act_z = 1;
            osp_updateQuantifyWindow(gui);
       case 6 %Overview tab
            gui.layout.ListBox.Enable = 'off';
            gui.controls.act_x = Exp;
            gui.controls.b_left_x.Enable = 'off';
            gui.controls.b_left_z.Enable = 'off';
            gui.controls.b_right_x.Enable = 'off';
            gui.controls.b_right_z.Enable = 'off';
            gui.info.nXvoxels = 1;
            gui.controls.act_x = 1;
            if ~strcmp(MRSCont.opts.fit.method, 'Osprey_gLCM')
                gui.info.nZvoxels = size(MRSCont.fit.results.metab.fitParams,1);
            else
                gui.info.nZvoxels = size(MRSCont.fit.results.metab,1);
            end
            gui.controls.act_z = 1;
            if ~strcmp(MRSCont.opts.fit.method, 'Osprey_gLCM')
                if size(MRSCont.fit.results.metab.fitParams,1) > 1
                    gui.controls.b_left_z.Enable = 'on';
                    gui.controls.b_right_z.Enable = 'on';
                end
            else
                if size(MRSCont.fit.results.metab,1) > 1
                    gui.controls.b_left_z.Enable = 'on';
                    gui.controls.b_right_z.Enable = 'on';
                end
            end
    end
end