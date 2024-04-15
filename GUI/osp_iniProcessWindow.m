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
        gui.layout.metabProTab = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
        gui.layout.proTab.TabWidth   = 90;
        gui.layout.proTabhandles = {'metabProTab'};
% Set up tabs with regard to the sequence type
         
        TabTitles = {'metab'};
        gui.process.metab.SNR = MRSCont.processed.metab{1,1}.QC_names;
        gui.process.metab.name = MRSCont.processed.metab{1,1}.names;
        if MRSCont.flags.hasMM
            TabTitles = horzcat(TabTitles, 'mm');
            gui.process.mm.SNR = MRSCont.processed.mm{1,1}.QC_names;
            gui.process.mm.name = MRSCont.processed.mm{1,1}.names;
        end
        if MRSCont.flags.hasRef
            TabTitles = horzcat(TabTitles, 'ref');
            gui.process.ref.SNR = MRSCont.processed.ref{1,1}.QC_names;
            gui.process.ref.name = MRSCont.processed.ref{1,1}.names;
        end
        if MRSCont.flags.hasMMRef
            TabTitles = horzcat(TabTitles, 'mm_ref');
            gui.process.mm_ref.SNR = MRSCont.processed.mm_ref{1,1}.QC_names;
            gui.process.mm_ref.name = MRSCont.processed.mm_ref{1,1}.names;
        end
        if MRSCont.flags.hasWater
            TabTitles = horzcat(TabTitles, 'w');
            gui.process.w.SNR = MRSCont.processed.w{1,1}.QC_names;
            gui.process.w.name = MRSCont.processed.w{1,1}.names;
        end
        for t = 2 : length(TabTitles)
            gui.layout.(strcat(TabTitles{t},'ProTab')) = uix.VBox('Parent', gui.layout.proTab,'BackgroundColor',gui.colormap.Background,'Spacing',5);
        end
        gui.layout.proTabhandles = strcat(TabTitles,'ProTab');
        gui.process.TabTitles = TabTitles;
        gui.layout.proTab.TabTitles  = TabTitles;
%%% 3. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        for t = length(gui.layout.proTabhandles) : -1 : 1 %Loop over subspecs/tabs
%             ind=find(ismember(gui.layout.proTabhandles,[gui.process.TabTitles{t} 'ProTab']));
            gui.layout.proTab.Selection  = t;
            gui.layout.proTab.TabEnables{t} = 'on';
            gui.upperBox.pro.box{t} = uix.HBox('Parent', gui.layout.(gui.layout.proTabhandles{t}),'BackgroundColor',gui.colormap.Background,'Spacing',5);
            if  (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI) 
                gui.upperBox.data.upperLeftButtons = uix.Panel('Parent', gui.upperBox.pro.box{t}, ...
                                         'Padding', 5, 'Title', ['Navigate voxel'],...
                                         'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                         'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
                gui.controls.Buttonbox = uix.HBox('Parent',gui.upperBox.data.upperLeftButtons, 'BackgroundColor',gui.colormap.Background);
                gui.controls.navigate_RawTab = uix.Grid('Parent',gui.controls.Buttonbox,'BackgroundColor',gui.colormap.Background);
                gui.controls.text_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','X:',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Y:',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Z:',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
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
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_act_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                gui.controls.text_act_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
                    'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
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
           
            
            
            gui.upperBox.pro.upperLeftButtons{t} = uix.Panel('Parent', gui.upperBox.pro.box{t}, ...
            'Padding', 5, 'Title', ['Navigate data'],...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
            'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
        gui.controls.Buttonbox{t} = uix.HBox('Parent',gui.upperBox.pro.upperLeftButtons{t}, 'BackgroundColor',gui.colormap.Background);
        gui.controls.navigate_RawTab{t} = uix.Grid('Parent',gui.controls.Buttonbox{t},'BackgroundColor',gui.colormap.Background);

        gui.controls.text_x = uicontrol(gui.controls.navigate_RawTab{t},'Style','text','String','Exp:',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.controls.text_y = uicontrol(gui.controls.navigate_RawTab{t},'Style','text','String','Spec:',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.controls.b_left_x = uicontrol(gui.controls.navigate_RawTab{t},'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
        gui.controls.b_left_y = uicontrol(gui.controls.navigate_RawTab{t},'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
        set(gui.controls.b_left_x,'Callback',{@osp_onLeftX,gui});
        set(gui.controls.b_left_y,'Callback',{@osp_onLeftY,gui});


        gui.controls.text_act_x = uicontrol(gui.controls.navigate_RawTab{t},'Style','text','String','1',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.controls.text_act_y = uicontrol(gui.controls.navigate_RawTab{t},'Style','text','String','1',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

        gui.controls.b_right_x = uicontrol(gui.controls.navigate_RawTab{t},'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
        gui.controls.b_right_y = uicontrol(gui.controls.navigate_RawTab{t},'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
        set(gui.controls.b_right_x,'Callback',{@osp_onRightX,gui});
        set(gui.controls.b_right_y,'Callback',{@osp_onRightY,gui});

        gui.controls.b_left_x.Enable = 'off';
        gui.controls.b_left_y.Enable = 'off';
        gui.controls.b_right_x.Enable = 'off';
        gui.controls.b_right_y.Enable = 'off';

        buttonString = [num2str(MRSCont.processed.(TabTitles{t}){1}.dims.extras > 1) num2str(MRSCont.processed.(TabTitles{t}){1}.subspecs>1)];
        switch buttonString
                case '01' 
                    gui.controls.b_left_y.Enable = 'on';
                    gui.controls.b_right_y.Enable = 'on';
                case '10' 
                    gui.controls.b_left_x.Enable = 'on';
                    gui.controls.b_right_x.Enable = 'on';
                case '11' 
                    gui.controls.b_left_x.Enable = 'on';
                    gui.controls.b_right_x.Enable = 'on';
                    gui.controls.b_left_y.Enable = 'on';
                    gui.controls.b_right_y.Enable = 'on';                                                 
        end
        set( gui.controls.navigate_RawTab{t}, 'Widths', [-20 -30 -20 -30], 'Heights', [-50 -50]);
    
            
            gui.upperBox.pro.Info{t} = uix.Panel('Parent', gui.upperBox.pro.box{t}, ...
                'Padding', 5, 'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}],...
                'HighlightColor', gui.colormap.Foreground,'FontName', gui.font, 'BackgroundColor',...
                gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            gui.upperBox.pro.upperButtons = uix.Panel('Parent', gui.upperBox.pro.box{t}, ...
                                     'Padding', 5, 'Title', ['Save'],...
                                     'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
                                     'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
            gui.controls.b_save_proTab{t} = uicontrol('Parent',gui.upperBox.pro.upperButtons,'Style','PushButton');
            [img, ~, ~] = imread('Printer.png', 'BackgroundColor', gui.colormap.Background);
            [img2] = imresize(img, 0.1);
            set(gui.controls.b_save_proTab{t},'CData', img2, 'TooltipString', 'Create EPS figure from current file');
            set(gui.controls.b_save_proTab{t},'Callback',{@osp_onPrint,gui});
            set(gui.upperBox.pro.box{t}, 'Width', [-0.16 -0.74 -0.1]);

            % Creates layout for plotting and data control
            gui.Plot.pro{t} = uix.HBox('Parent', gui.layout.(gui.layout.proTabhandles{t}), ...
                'Padding', 5,'BackgroundColor', gui.colormap.Background);
            set(gui.layout.(gui.layout.proTabhandles{t}), 'Heights', [-0.1 -0.9]);
            % Get parameter from file to fill the info panel
            Exp = 1;
            if ~(isfield(MRSCont.flags,'isPRIAM')|| isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                if strcmp(gui.process.TabTitles{t},'metab')
                    StatText = ['SNR(' gui.process.(gui.process.TabTitles{t}).SNR{1} '): '  num2str(MRSCont.QM.SNR.(gui.process.TabTitles{t})(Exp,gui.controls.Selected)) '; FWHM (' gui.process.(gui.process.TabTitles{t}).SNR{1} '): '...
                                num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{t})(Exp,gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{t})(Exp,gui.controls.Selected)/MRSCont.processed.(gui.process.TabTitles{t}){gui.controls.Selected}.txfrq(Exp)*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM.freqShift.(gui.process.TabTitles{t})(Exp,gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.A(gui.controls.Selected)*MRSCont.processed.(gui.process.TabTitles{t}){gui.controls.Selected}.txfrq(Exp)/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.A(gui.controls.Selected)*MRSCont.processed.(gui.process.TabTitles{t}){gui.controls.Selected}.txfrq(Exp)/1e6) ' Hz'];
                else
                StatText = ['SNR(' gui.process.(gui.process.TabTitles{t}).SNR{1} '): ' num2str(MRSCont.QM.SNR.(gui.process.TabTitles{t})(Exp,gui.controls.Selected)) '; FWHM (' gui.process.(gui.process.TabTitles{t}).SNR{1} '): '...
                            num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{t})(Exp,gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{t})(Exp,gui.controls.Selected)/MRSCont.processed.(gui.process.TabTitles{t}){gui.controls.Selected}.txfrq(Exp)*1e6))...
                            ' Hz / ppm'];
                end
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                if  (contains(gui.process.TabTitles{t},{'A','B','C','D','diff1','diff2','sum'}))
                    StatText = ['Voxel ' num2str(gui.controls.act_x) ': SNR(' gui.process.SNR{t} '): '  num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.TabTitles{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                                num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.TabTitles{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.TabTitles{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.StructFields{t}){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM{1,gui.controls.act_x}.freqShift.(gui.process.TabTitles{t})(gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM{1,gui.controls.act_x}.drift.pre.AvgDeltaCr.(gui.process.TabTitles{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.StructFields{t}){gui.controls.Selected}.txfrq/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM{1,gui.controls.act_x}.drift.post.AvgDeltaCr.(gui.process.TabTitles{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.StructFields{t}){gui.controls.Selected}.txfrq/1e6) ' Hz'];
                else
                StatText = ['Voxel ' num2str(gui.controls.act_x) ':SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.TabTitles{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                            num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.TabTitles{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.TabTitles{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.StructFields{t}){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm'];
                end
            else
                if  (contains(gui.process.TabTitles{t},{'A','B','C','D','diff1','diff2','sum'}))
                    StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': SNR(' gui.process.SNR{t} '): '  num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(gui.process.TabTitles{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                                num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.TabTitles{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.TabTitles{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.StructFields{t}){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.freqShift.(gui.process.TabTitles{t})(gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.drift.pre.AvgDeltaCr.(gui.process.TabTitles{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.StructFields{t}){gui.controls.Selected}.txfrq/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.drift.post.AvgDeltaCr.(gui.process.TabTitles{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.StructFields{t}){gui.controls.Selected}.txfrq/1e6) ' Hz'];
                else 
                StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(gui.process.TabTitles{t})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{t} '): '...
                            num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.TabTitles{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.TabTitles{t})(gui.controls.Selected)/MRSCont.processed.(gui.process.StructFields{t}){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm'];
                end
            end
            gui.InfoText.pro{t}  = uicontrol('Parent',gui.upperBox.pro.Info{t},'style','text',...
                                         'FontSize', 12, 'FontName', gui.font,...
                                         'HorizontalAlignment', 'left', 'String', sprintf(StatText),...
                                         'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

 %%% 4. VISUALIZATION PART OF THIS TAB %%%
 %osp_plotProcess is used to visualize the processed spectra
            temp = osp_plotProcess(MRSCont, gui.controls.Selected,gui.process.TabTitles{t},1,Exp); % Create figure
            %Subplots are distributed here
                gui.layout.proSpecs = uix.VBox('Parent', gui.Plot.pro{t}, 'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.layout.proPre = uix.VBox('Parent', gui.layout.proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    gui.layout.proPost = uix.VBox('Parent', gui.layout.proSpecs,'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                gui.layout.proOut = uix.VBox('Parent', gui.Plot.pro{t},'Padding', 5, 'BackgroundColor',gui.colormap.Background);
                    gui.layout.proDrift = uix.VBox('Parent', gui.layout.proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);
                    gui.layout.proAlgn = uix.VBox('Parent', gui.layout.proOut, 'Padding', 5,'Units', 'Normalized', 'BackgroundColor',gui.colormap.Background);

            set( temp.Children(1), 'Parent', gui.layout.proDrift );
            set( temp.Children(1), 'Parent', gui.layout.proAlgn );
            set( temp.Children(1), 'Parent', gui.layout.proPost );
            set( temp.Children(1), 'Parent', gui.layout.proPre );
            close( temp );
%%% 5. DATA CONTROLS FOR THIS TAB %%%
            set(gui.Plot.pro{t},'Widths', [-0.49 -0.49]);
            set(gui.layout.proPre.Children(1), 'Units', 'normalized')
            set(gui.layout.proPre.Children(1), 'OuterPosition', [0,0,1,1])
            set(gui.layout.proPost.Children(1), 'Units', 'normalized')
            set(gui.layout.proPost.Children(1), 'OuterPosition', [0,0,1,1])
            set(gui.layout.proDrift.Children(1), 'Units', 'normalized')
            set(gui.layout.proDrift.Children(1), 'OuterPosition', [0,0,1,1])
            set(gui.layout.proAlgn.Children(1), 'Units', 'normalized')
            set(gui.layout.proAlgn.Children(1), 'OuterPosition', [0,0,1,1])
        end
    h = findall(groot,'Type','figure');
    for ff = 1 : length(h)
        if ~(strcmp(h(ff).Tag, 'Osprey') ||  strcmp(h(ff).Tag, 'TMWWaitbar'))
            close(h(ff))
        end
    end
    gui.controls.YLimits = [];
    gui.info.nYvoxels = MRSCont.processed.metab{1, 1}.subspecs;
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end