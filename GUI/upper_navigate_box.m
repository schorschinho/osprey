function gui = upper_navigate_box(gui, name,Selection,Overview)
        MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class
    
        gui.upperBox.(name).upperLeftButtons = uix.Panel('Parent', gui.upperBox.(name).box, ...
            'Padding', 5, 'Title', ['Navigate model'],...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground,...
            'HighlightColor', gui.colormap.Foreground, 'ShadowColor', gui.colormap.Foreground);
        gui.controls.Buttonbox = uix.HBox('Parent',gui.upperBox.(name).upperLeftButtons, 'BackgroundColor',gui.colormap.Background);
        gui.controls.navigate_RawTab = uix.Grid('Parent',gui.controls.Buttonbox,'BackgroundColor',gui.colormap.Background);

        gui.controls.text_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Exp',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        if ~Overview
            gui.controls.text_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Spec',...
                'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        end
        gui.controls.text_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','Basis',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        gui.controls.b_left_x = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
        if ~Overview
            gui.controls.b_left_y = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
        end
        gui.controls.b_left_z = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','<');
        set(gui.controls.b_left_x,'Callback',{@osp_onLeftX,gui});
        if ~Overview
            set(gui.controls.b_left_y,'Callback',{@osp_onLeftY,gui});
        end
        set(gui.controls.b_left_z,'Callback',{@osp_onLeftZ,gui});


        gui.controls.text_act_x = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        if ~Overview
            gui.controls.text_act_y = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
                'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        end
        gui.controls.text_act_z = uicontrol(gui.controls.navigate_RawTab,'Style','text','String','1',...
            'FontName', gui.font, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);

        gui.controls.b_right_x = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
        if ~Overview
            gui.controls.b_right_y = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
        end
        gui.controls.b_right_z = uicontrol(gui.controls.navigate_RawTab,'Style','PushButton', 'BackgroundColor',gui.colormap.Background,'String','>');
        set(gui.controls.b_right_x,'Callback',{@osp_onRightX,gui});
        if ~Overview
            set(gui.controls.b_right_y,'Callback',{@osp_onRightY,gui});
        end
        set(gui.controls.b_right_z,'Callback',{@osp_onRightZ,gui});

        gui.controls.b_left_x.Enable = 'off';
        if ~Overview
            gui.controls.b_left_y.Enable = 'off';
        end
        gui.controls.b_left_z.Enable = 'off';
        gui.controls.b_right_x.Enable = 'off';
        if ~Overview
            gui.controls.b_right_y.Enable = 'off';
        end
        gui.controls.b_right_z.Enable = 'off';

        if ~Overview
            if ~strcmp(MRSCont.opts.fit.method, 'Osprey_gLCM')
                buttonString = [num2str(MRSCont.nDatasets(2) > 1) num2str(size(MRSCont.fit.results.(Selection).fitParams,3)>1) num2str(size(MRSCont.fit.results.(Selection).fitParams,1)>1)];
            else
                buttonString = [num2str(MRSCont.nDatasets(2) > 1) num2str(size(MRSCont.fit.results.(Selection),3)>1) num2str(size(MRSCont.fit.results.(Selection),2)>1)];
            end
            switch buttonString
                    case '001' 
                        gui.controls.b_left_z.Enable = 'on';
                        gui.controls.b_right_z.Enable = 'on';
                    case '010' 
                        gui.controls.b_left_y.Enable = 'on';
                        gui.controls.b_right_y.Enable = 'on';
                    case '100' 
                        gui.controls.b_left_x.Enable = 'on';
                        gui.controls.b_right_x.Enable = 'on';
                    case '011' 
                        gui.controls.b_left_y.Enable = 'on';
                        gui.controls.b_right_y.Enable = 'on';
                        gui.controls.b_left_z.Enable = 'on';
                        gui.controls.b_right_z.Enable = 'on';  
                   case '101' 
                        gui.controls.b_left_x.Enable = 'on';
                        gui.controls.b_right_x.Enable = 'on';
                        gui.controls.b_left_z.Enable = 'on';
                        gui.controls.b_right_z.Enable = 'on';
                   case '110' 
                        gui.controls.b_left_x.Enable = 'on';
                        gui.controls.b_right_x.Enable = 'on';
                        gui.controls.b_left_y.Enable = 'on';
                        gui.controls.b_right_y.Enable = 'on'; 
                    case '111'
                        gui.controls.b_left_x.Enable = 'on';
                        gui.controls.b_left_y.Enable = 'on';
                        gui.controls.b_left_z.Enable = 'on';
                        gui.controls.b_right_x.Enable = 'on';
                        gui.controls.b_right_y.Enable = 'on';
                        gui.controls.b_right_z.Enable = 'on';                                
            end
            set( gui.controls.navigate_RawTab, 'Widths', [-20 -30 -20 -30], 'Heights', [-33 -33 -33]);
        else
            if ~strcmp(MRSCont.opts.fit.method, 'Osprey_gLCM')
                buttonString = [num2str(MRSCont.nDatasets(2) > 1) num2str(size(MRSCont.fit.results.(Selection).fitParams,1)>1) ];
            else
                buttonString = [num2str(MRSCont.nDatasets(2) > 1) num2str(size(MRSCont.fit.results.(Selection),1)>1) ];
            end
            switch buttonString
                    case '01' 
                        gui.controls.b_left_z.Enable = 'on';
                        gui.controls.b_right_z.Enable = 'on';
                    case '10' 
                        gui.controls.b_left_x.Enable = 'on';
                        gui.controls.b_right_x.Enable = 'on';
                    case '11' 
                        gui.controls.b_left_x.Enable = 'on';
                        gui.controls.b_right_x.Enable = 'on';
                        gui.controls.b_left_z.Enable = 'on';
                        gui.controls.b_right_z.Enable = 'on';  
            end
            set( gui.controls.navigate_RawTab, 'Widths', [-20 -30 -20 -30], 'Heights', [-50 -50]);
        end
        
end