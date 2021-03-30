function osp_iniProcessWindow(gui)
%% osp_iniProcessWindow
%   This function creates the initial processed window in the gui.
%
%
%   USAGE:
%       osp_iniProcessWindow(gui);
%
%   INPUT:      gui      = gui class containing all handles and the MRSCont 
%
%   OUTPUT:     Changes in gui parameters and MRSCont are written into the
%               gui class
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
%This functions creates the initial process window    
        MRSCont = getappdata(gui.figure,'MRSCont');   % Get MRSCont from hidden container in gui class
        gui.layout.tabs.TabEnables{2} = 'on';
        gui.layout.tabs.Selection  = 2;
        gui.layout.EmptyProPlot = 0;
%%% 2. CREATING SUB TABS FOR THIS TAB %%% 
% In this case one tab fo each subspec (A,B,C,D,ref,water)
        gui.layout.AProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
        gui.layout.proTab.TabWidth   = 90;
        gui.layout.proTabhandles = {'AProTab'};
% Set up tabs with regard to the sequence type
        if MRSCont.flags.isUnEdited %Is UnEdited?
        if (MRSCont.flags.hasRef && MRSCont.flags.hasWater && MRSCont.flags.hasMM ) %Has all%re_mm
                gui.layout.mmProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);%re_mm
                gui.layout.refProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);%re_mm
                gui.layout.wProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);%re_mm
                gui.layout.proTab.TabTitles  = {'A', 'MM','ref','w'};%re_mm
                gui.layout.proTab.TabEnables = {'on', 'on','on','on'};%re_mm
                gui.layout.proTabhandles = {'AProTab','mmProTab', 'refProTab', 'wProTab'}; % Create 4 Tabs %re_mm
                gui.process.SNR = {'tNAA','MM09','water','water'};%re_mm
        elseif (~MRSCont.flags.hasRef && MRSCont.flags.hasWater && MRSCont.flags.hasMM ) %Has all%re_mm
                gui.layout.mmProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);%re_mm
                gui.layout.wProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);%re_mm
                gui.layout.proTab.TabTitles  = {'A', 'MM','w'};%re_mm
                gui.layout.proTab.TabEnables = {'on', 'on','on'};%re_mm
                gui.layout.proTabhandles = {'AProTab','mmProTab', 'wProTab'}; % Create 4 Tabs %re_mm
                gui.process.SNR = {'tNAA','MM09','water'};%re_mm
        elseif (MRSCont.flags.hasRef && ~MRSCont.flags.hasWater && MRSCont.flags.hasMM ) %Has all%re_mm
                gui.layout.mmProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);%re_mm
                gui.layout.refProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);%re_mm
                gui.layout.proTab.TabTitles  = {'A', 'MM','w'};%re_mm
                gui.layout.proTab.TabEnables = {'on', 'on','on'};%re_mm
                gui.layout.proTabhandles = {'AProTab','mmProTab', 'refProTab'}; % Create 4 Tabs %re_mm
                gui.process.SNR = {'tNAA','MM09','water'};%re_mm        
        elseif (MRSCont.flags.hasRef && MRSCont.flags.hasWater && ~MRSCont.flags.hasMM) %Has water and reference (This should actually not happen with UnEdited data...) 
                gui.layout.refProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.wProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.proTab.TabTitles  = {'A','ref','w'};
                gui.layout.proTab.TabEnables = {'on','on','on'};
                gui.layout.proTabhandles = {'AProTab', 'refProTab', 'wProTab'}; % Create 3 Tabs
                gui.process.SNR = {'tNAA','water','water'};
            elseif (~MRSCont.flags.hasRef && ~MRSCont.flags.hasWater && ~MRSCont.flags.hasMM) %Only metabolite data                
                gui.layout.proTab.TabTitles  = {'A'};
                gui.layout.proTab.TabEnables = {'on'};
                gui.layout.proTabhandles = {'AProTab'}; % Create 1 Tab
                gui.process.SNR = {'tNAA'};
            else
                if MRSCont.flags.hasRef %Has only reference?
                    gui.layout.refProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.proTab.TabTitles  = {'A','ref'};
                    gui.layout.proTab.TabEnables = {'on', 'on'};
                    gui.layout.proTabhandles = {'AProTab', 'refProTab'}; % Create 2 Tabs
                    gui.process.SNR = {'tNAA','water'};
                end
                if MRSCont.flags.hasWater %Has only water?
                    gui.layout.wProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.proTab.TabTitles  = {'A','w'};
                    gui.layout.proTab.TabEnables = {'on', 'on'};
                    gui.layout.proTabhandles = {'AProTab', 'wProTab'}; %Create 2 Tabs
                    gui.process.SNR = {'tNAA','water'};
                end
            end
        end
        if MRSCont.flags.isMEGA %Is MEGA?
            if (MRSCont.flags.hasRef && MRSCont.flags.hasWater) %Has water and reference?
                gui.layout.BProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.diff1ProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.sumProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.refProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.wProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.proTab.TabTitles  = {'A','B','diff1','sum','ref','w'};
                gui.layout.proTab.TabEnables = {'on', 'on','on', 'on', 'on', 'on'};
                gui.layout.proTabhandles = {'AProTab','BProTab', 'diff1ProTab','sumProTab', 'refProTab', 'wProTab'}; %Create 6 tabs
                gui.process.SNR = {'tNAA','tCr',MRSCont.processed.diff1{1,gui.process.Selected}.target,'tNAA','water','water'};
                if iscell(gui.process.SNR{3})
                        gui.process.SNR{3} = gui.process.SNR{3}{1};
                end
            elseif (~MRSCont.flags.hasRef && ~MRSCont.flags.hasWater) %Only metabolites?
                gui.layout.BProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.diff1ProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.sumProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.proTab.TabTitles  = {'A','B','diff1','sum'};
                gui.layout.proTab.TabEnables = {'on', 'on','on', 'on'};
                gui.layout.proTabhandles = {'AProTab','BProTab', 'diff1ProTab','sumProTab',}; %Create 4 tabs
                gui.process.SNR = {'tNAA','tCr',MRSCont.processed.diff1{1,gui.process.Selected}.target,'water'};
                if iscell(gui.process.SNR{3})
                        gui.process.SNR{3} = gui.process.SNR{3}{1};
                end
            else
                if MRSCont.flags.hasRef %Has only reference?
                    gui.layout.BProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.diff1ProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.sumProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.refProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.proTab.TabTitles  = {'A','B','diff1','sum','ref'};
                    gui.layout.proTab.TabEnables = {'on', 'on','on', 'on', 'on'};
                    gui.layout.proTabhandles = {'AProTab','BProTab', 'diff1ProTab','sumProTab', 'refProTab'}; %Create 5 tabs
                    gui.process.SNR = {'tNAA','tCr',MRSCont.processed.diff1{1,gui.process.Selected}.target,'tNAA','water'};
                    if iscell(gui.process.SNR{3})
                            gui.process.SNR{3} = gui.process.SNR{3}{1};
                    end
                end
                if MRSCont.flags.hasWater %Has only water?
                    gui.layout.BProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.diff1ProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.sumProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.wProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.proTab.TabTitles  = {'A','B','diff1','sum','w'};
                    gui.layout.proTab.TabEnables = {'on', 'on','on', 'on', 'on'};
                    gui.layout.proTabhandles = {'AProTab','BProTab', 'diff1ProTab','sumProTab', 'wProTab'}; %Create 5 tabs
                    gui.process.SNR = {'tNAA','tCr',MRSCont.processed.diff1{1,gui.process.Selected}.target,'tNAA','water'};
                    if iscell(gui.process.SNR{3})
                            gui.process.SNR{3} = gui.process.SNR{3}{1};
                    end
                end
            end
        end
        if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) %Is HERMES\HERCULES?
            if (MRSCont.flags.hasRef && MRSCont.flags.hasWater) %Has water and reference?
                gui.layout.BProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.CProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.DProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.diff1ProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.diff2ProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.sumProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.refProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.wProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.proTab.TabTitles  = {'A','B','C','D','diff1', 'diff2','sum','ref','w'};
                gui.layout.proTab.TabEnables = {'on', 'on','on', 'on', 'on', 'on', 'on', 'on', 'on'};
                gui.layout.proTabhandles = {'AProTab','BProTab','CProTab','DProTab', 'diff1ProTab', 'diff2ProTab', 'sumProTab', 'refProTab', 'wProTab'}; %Create 9 tabs
                try
                    gui.process.SNR = {'tNAA','tCr','tNAA', 'tCr', MRSCont.processed.diff1{1,gui.process.Selected}.target{1},MRSCont.processed.diff2{1,gui.process.Selected}.target{1},'tNAA','water','water'};
                catch
                   gui.process.SNR = {'tNAA','tCr','tNAA', 'tCr', MRSCont.processed.diff1{1,gui.process.Selected}.target,MRSCont.processed.diff2{1,gui.process.Selected}.target,'tNAA','water','water'}; 
                end
            elseif (~MRSCont.flags.hasRef && ~MRSCont.flags.hasWater) %Only metabs?
                gui.layout.BProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.CProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.DProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.diff1ProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.diff2ProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.sumProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                gui.layout.proTab.TabTitles  = {'A','B','C','D','diff1', 'diff2','sum'};
                gui.layout.proTab.TabEnables = {'on', 'on','on', 'on', 'on', 'on', 'on'};
                gui.layout.proTabhandles = {'AProTab','BProTab','CProTab','DProTab', 'diff1ProTab', 'diff2ProTab', 'sumProTab'}; %Create 7 tabs
                try
                    gui.process.SNR = {'tNAA','tCr','tNAA', 'tCr', MRSCont.processed.diff1{1,gui.process.Selected}.target{1},MRSCont.processed.diff2{1,gui.process.Selected}.target{1},'tNAA'};
                catch
                    gui.process.SNR = {'tNAA','tCr','tNAA', 'tCr', MRSCont.processed.diff1{1,gui.process.Selected}.target,MRSCont.processed.diff2{1,gui.process.Selected}.target,'tNAA'};
                end
            else
                if MRSCont.flags.hasRef %Has only reference?
                    gui.layout.BProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.CProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.DProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.diff1ProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.diff2ProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.sumProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);                    
                    gui.layout.refProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.proTab.TabTitles  = {'A','B','C','D','diff1', 'diff2','sum','ref'};
                    gui.layout.proTab.TabEnables = {'on', 'on','on','on', 'on', 'on', 'on', 'on'};
                    gui.layout.proTabhandles = {'AProTab','BProTab','CProTab','DProTab','diff1ProTab', 'diff2ProTab', 'sumProTab','refProTab'}; %Create 8 tabs
                    try
                        gui.process.SNR = {'tNAA','tCr','tNAA', 'tCr', MRSCont.processed.diff1{1,gui.process.Selected}.target{1},MRSCont.processed.diff2{1,gui.process.Selected}.target{1},'tNAA','water'};
                    catch
                        gui.process.SNR = {'tNAA','tCr','tNAA', 'tCr', MRSCont.processed.diff1{1,gui.process.Selected}.target,MRSCont.processed.diff2{1,gui.process.Selected}.target,'tNAA','water'};
                    end
                end
                if MRSCont.flags.hasWater %Has only water?
                    gui.layout.BProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.CProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.DProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.diff1ProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.diff2ProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.sumProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);                    
                    gui.layout.wProTab = uix.VBox('Parent', gui.layout.proTab, 'Padding', 5,'BackgroundColor',gui.colormap.Background,'Spacing',5);
                    gui.layout.proTab.TabTitles  = {'A','B','C','D','diff1', 'diff2','sum','w'};
                    gui.layout.proTab.TabEnables = {'on', 'on','on','on', 'on', 'on', 'on', 'on'};
                    gui.layout.proTabhandles = {'AProTab','BProTab','CProTab','DProTab','diff1ProTab', 'diff2ProTab', 'sumProTab','wProTab'}; %Create 8 tabs
                    try
                        gui.process.SNR = {'tNAA','tCr','tNAA', 'tCr', MRSCont.processed.diff1{1,gui.process.Selected}.target{1},MRSCont.processed.diff2{1,gui.process.Selected}.target{1},'tNAA','water'};
                    catch
                        gui.process.SNR = {'tNAA','tCr','tNAA', 'tCr', MRSCont.processed.diff1{1,gui.process.Selected}.target,MRSCont.processed.diff2{1,gui.process.Selected}.target,'tNAA','water'};
                    end
                end
            end
        end

%%% 3. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        for t = length(gui.layout.proTabhandles) : -1 : 1 %Loop over subspecs/tabs
            ind=find(ismember(gui.layout.proTabhandles,[gui.process.Names{t} 'ProTab']));
            gui.layout.proTab.Selection  = ind;
            gui.upperBox.pro.box = uix.HBox('Parent', gui.layout.(gui.layout.proTabhandles{ind}),'BackgroundColor',gui.colormap.Background,'Spacing',5);
            if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                gui.upperBox.data.upperLeftButtons = uix.Panel('Parent', gui.upperBox.pro.box, ...
                                         'Padding', 5, 'Title', ['Navigate voxel'],...
                                         'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                         'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
                gui.controls.Buttonbox = uix.HBox('Parent',gui.upperBox.data.upperLeftButtons, 'BackgroundColor',gui.colormap.Background);
                gui.controls.navigate_RawTab = uix.Grid('Parent',gui.controls.Buttonbox,'BackgroundColor',gui.colormap.Background);
                gui.controls.text_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','X:',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Y:',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Z:',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.b_left_x = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
                gui.controls.b_left_y = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
                gui.controls.b_left_z = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
                set(gui.controls.b_left_x,'Callback',{@osp_onLeftX,gui});
                set(gui.controls.b_left_y,'Callback',{@osp_onLeftY,gui});
                set(gui.controls.b_left_z,'Callback',{@osp_onLeftZ,gui});
                if gui.info.nXvoxels <= 1
                    gui.controls.b_left_x.Enable = 'off';
                end
                if gui.info.nYvoxels <= 1
                    gui.controls.b_left_y.Enable = 'off';
                end
                if gui.info.nZvoxels <= 1
                    gui.controls.b_left_z.Enable = 'off';
                end
                
                gui.controls.text_act_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_act_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_act_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
                    'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.b_right_x = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
                gui.controls.b_right_y = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
                gui.controls.b_right_z = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
                set(gui.controls.b_right_x,'Callback',{@osp_onRightX,gui});
                set(gui.controls.b_right_y,'Callback',{@osp_onRightY,gui});
                set(gui.controls.b_right_z,'Callback',{@osp_onRightZ,gui});
                if gui.info.nXvoxels <= 1
                    gui.controls.b_right_x.Enable = 'off';
                end
                if gui.info.nYvoxels <= 1
                    gui.controls.b_right_y.Enable = 'off';
                end
                if gui.info.nZvoxels <= 1
                    gui.controls.b_right_z.Enable = 'off';
                end   
                set( gui.controls.navigate_RawTab, 'Widths', [-20 -30 -20 -30], 'Heights', [-33 -33 -33] );
            end
            gui.upperBox.pro.Info = uix.Panel('Parent', gui.upperBox.pro.box, ...
                'Padding', 5, 'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}],...
                'HighlightColor', gui.colormap.Foreground,'FontName', 'Arial', 'BackgroundColor',...
                gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            gui.upperBox.pro.upperButtons = uix.Panel('Parent', gui.upperBox.pro.box, ...
                                     'Padding', 5, 'Title', ['Save'],...
                                     'FontName', 'Arial', 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                     'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            gui.controls.b_save_proTab = uicontrol('Parent',gui.upperBox.pro.upperButtons,'Style','PushButton');
            [img, ~, ~] = imread('Printer.png', 'BackgroundColor', gui.colormap.Background);
            [img2] = imresize(img, 0.1);
            set(gui.controls.b_save_proTab,'CData', img2, 'TooltipString', 'Create EPS figure from current file');
            set(gui.controls.b_save_proTab,'Callback',{@osp_onPrint,gui});
            if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                set(gui.upperBox.pro.box, 'Width', [-0.12 -0.78 -0.1]);
            else
                set(gui.upperBox.pro.box, 'Width', [-0.9 -0.1]);   
            end
            % Creates layout for plotting and data control
            gui.Plot.pro = uix.HBox('Parent', gui.layout.(gui.layout.proTabhandles{ind}), ...
                'Padding', 5,'BackgroundColor', gui.colormap.Background);
            set(gui.layout.(gui.layout.proTabhandles{ind}), 'Heights', [-0.1 -0.9]);
            % Get parameter from file to fill the info panel
            if ~(isfield(MRSCont.flags,'isPRIAM')|| isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                if (strcmp(gui.process.Names{t},'A') || strcmp(gui.process.Names{t},'B') || strcmp(gui.process.Names{t},'C') || strcmp(gui.process.Names{t},'D') || strcmp(gui.process.Names{t},'diff1') || strcmp(gui.process.Names{t},'diff2') || strcmp(gui.process.Names{t},'sum'))
                    StatText = ['Metabolite Data -> SNR(' gui.process.SNR{t} '): '  num2str(MRSCont.QM.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                                num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM.freqShift.(gui.process.Names{t})(gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.(gui.process.Names{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.(gui.process.Names{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq/1e6) ' Hz'];
                else if strcmp(gui.process.Names{t},'ref')
                StatText = ['Reference Data -> SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                            num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm'];
                    else
                        if ~strcmp(gui.process.Names{t},'mm') %re
                        StatText = ['Water Data -> SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                                    num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)) '/' (num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                                    ' Hz / ppm'];
                        else
                            StatText = ['MM Data -> SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                                    num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)) '/' (num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                                    ' Hz / ppm'];
                        end %re
                    end
                end
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                if (strcmp(gui.process.Names{t},'A') || strcmp(gui.process.Names{t},'B') || strcmp(gui.process.Names{t},'C') || strcmp(gui.process.Names{t},'D') || strcmp(gui.process.Names{t},'diff1') || strcmp(gui.process.Names{t},'diff2') || strcmp(gui.process.Names{t},'sum'))
                    StatText = ['Voxel ' num2str(gui.controls.act_x) ': Metabolite Data -> SNR(' gui.process.SNR{t} '): '  num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                                num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM{1,gui.controls.act_x}.freqShift.(gui.process.Names{t})(gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM{1,gui.controls.act_x}.drift.pre.AvgDeltaCr.(gui.process.Names{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM{1,gui.controls.act_x}.drift.post.AvgDeltaCr.(gui.process.Names{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq/1e6) ' Hz'];
                else if strcmp(gui.process.Names{t},'ref')
                StatText = ['Voxel ' num2str(gui.controls.act_x) ': Reference Data -> SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                            num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm'];
                    else
                        if ~strcmp(gui.process.Names{t},'mm') %re
                        StatText = ['Voxel ' num2str(gui.controls.act_x) ': Water Data -> SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                                    num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{t})(gui.controls.Selected)) '/' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                                    ' Hz / ppm'];
                        else
                            StatText = ['Voxel ' num2str(gui.controls.act_x) ': MM Data -> SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                                    num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{t})(gui.controls.Selected)) '/' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                                    ' Hz / ppm'];
                        end %re
                    end
                end
            else
                if (strcmp(gui.process.Names{t},'A') || strcmp(gui.process.Names{t},'B') || strcmp(gui.process.Names{t},'C') || strcmp(gui.process.Names{t},'D') || strcmp(gui.process.Names{t},'diff1') || strcmp(gui.process.Names{t},'diff2') || strcmp(gui.process.Names{t},'sum'))
                    StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': Metabolite Data -> SNR(' gui.process.SNR{t} '): '  num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                                num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.Names{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.freqShift.(gui.process.Names{t})(gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.drift.pre.AvgDeltaCr.(gui.process.Names{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.drift.post.AvgDeltaCr.(gui.process.Names{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq/1e6) ' Hz'];
                else if strcmp(gui.process.Names{t},'ref')
                StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': Reference Data -> SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                            num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.Names{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm'];
                    else
                        if ~strcmp(gui.process.Names{t},'mm') %re
                        StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': Water Data -> SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                                    num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.Names{t})(gui.controls.Selected)) '/' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                                    ' Hz / ppm'];
                        else
                            StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': MM Data -> SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                                    num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.Names{t})(gui.controls.Selected)) '/' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.Names{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq*1e6))...
                                    ' Hz / ppm'];
                        end %re
                    end
                end
            end
            gui.InfoText.pro  = uicontrol('Parent',gui.upperBox.pro.Info,'style','text',...
                                         'FontSize', 12, 'FontName', 'Arial',...
                                         'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                         'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

 %%% 4. VISUALIZATION PART OF THIS TAB %%%
 %osp_plotProcess is used to visualize the processed spectra
            temp = osp_plotProcess(MRSCont, gui.controls.Selected,gui.process.Names{t}); % Create figure
            %Subplots are distributed here
                gui.layout.proSpecs = uix.VBox('Parent', gui.Plot.pro, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.layout.proPre = uix.VBox('Parent', gui.layout.proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    gui.layout.proPost = uix.VBox('Parent', gui.layout.proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                gui.layout.proOut = uix.VBox('Parent', gui.Plot.pro,'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.layout.proDrift = uix.VBox('Parent', gui.layout.proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    gui.layout.proAlgn = uix.VBox('Parent', gui.layout.proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);

            set( temp.Children(1), 'Parent', gui.layout.proDrift );
            set( temp.Children(1), 'Parent', gui.layout.proAlgn );
            set( temp.Children(1), 'Parent', gui.layout.proPost );
            set( temp.Children(1), 'Parent', gui.layout.proPre );
            close( temp );
%%% 5. DATA CONTROLS FOR THIS TAB %%%
            set(gui.Plot.pro,'Widths', [-0.49 -0.49]);
            set(gui.layout.proPre.Children(1), 'Units', 'normalized')
            set(gui.layout.proPre.Children(1), 'OuterPosition', [0,0,1,1])
            set(gui.layout.proPost.Children(1), 'Units', 'normalized')
            set(gui.layout.proPost.Children(1), 'OuterPosition', [0,0,1,1])
            set(gui.layout.proDrift.Children(1), 'Units', 'normalized')
            set(gui.layout.proDrift.Children(1), 'OuterPosition', [0,0,1,1])
            set(gui.layout.proAlgn.Children(1), 'Units', 'normalized')
            set(gui.layout.proAlgn.Children(1), 'OuterPosition', [0,0,1,1])
        end
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end