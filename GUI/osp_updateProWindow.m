function osp_updateProWindow(gui)
%% osp_updateProWindow
%   This function updates the process tab.
%
%
%   USAGE:
%       osp_updateProWindow(gui);
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
%%% 1. INITIALIZE %%%
        MRSCont = getappdata(gui.figure,'MRSCont');  % Get MRSCont from hidden container in gui class 
         if (isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) && (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            set(gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children(3).Children(1).Children.Children(4),'String',gui.controls.act_z);
            set(gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children(3).Children(1).Children.Children(5),'String',gui.controls.act_y);
            set(gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children(3).Children(1).Children.Children(6),'String',gui.controls.act_x);
         else
            set(gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children(3).Children(1).Children.Children(3),'String',gui.controls.act_y);
            set(gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children(3).Children(1).Children.Children(4),'String',gui.controls.act_x);
        end
        gui.layout.EmptyProPlot = 0;
        Selection = gui.process.TabTitles{gui.process.Selected};
        Exp = gui.controls.act_x;
        SubSpec = gui.controls.act_y;
        
        if (MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.dims.extras == 0) || (Exp > MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.sz(MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.dims.extras))
            Exp = 1;
            gui.controls.act_x = 1;
        end
        if (MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.dims.subSpecs == 0) || (SubSpec > MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.sz(MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.dims.subSpecs))
            SubSpec = 1;
            gui.controls.act_y = 1;
        end
                        
%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        if ~(isfield(MRSCont.flags,'isPRIAM')|| isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                 if strcmp(gui.process.TabTitles{gui.process.Selected},'metab')
                    StatText = ['SNR(' gui.process.(gui.process.TabTitles{gui.process.Selected}).SNR{SubSpec} '): '  num2str(MRSCont.QM.SNR.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected,SubSpec)) '; FWHM (' gui.process.(gui.process.TabTitles{gui.process.Selected}).SNR{1} '): '...
                                num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected,SubSpec)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected,SubSpec)/MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM.freqShift.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected,SubSpec)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.A(gui.controls.Selected)*MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.A(gui.controls.Selected)*MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)/1e6) ' Hz'];
                else
                StatText = ['SNR(' gui.process.(gui.process.TabTitles{gui.process.Selected}).SNR{1} '): ' num2str(MRSCont.QM.SNR.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) '; FWHM (' gui.process.(gui.process.TabTitles{gui.process.Selected}).SNR{1} '): '...
                            num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)/MRSCont.processed.(gui.process.TabTitles{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)*1e6))...
                            ' Hz / ppm'];
                end
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                if  (contains(gui.process.TabTitles{gui.process.Selected},{'A','B','C','D','diff1','diff2','sum'}))
                    StatText = ['Voxel ' num2str(gui.controls.act_x) ': SNR(' gui.process.SNR{gui.process.Selected} '): '  num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                                num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)/MRSCont.processed.(gui.process.StructFields{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM{1,gui.controls.act_x}.freqShift.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM{1,gui.controls.act_x}.drift.pre.AvgDeltaCr.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)*MRSCont.processed.(gui.process.StructFields{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM{1,gui.controls.act_x}.drift.post.AvgDeltaCr.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)*MRSCont.processed.(gui.process.StructFields{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)/1e6) ' Hz'];
                else
                StatText = ['Voxel ' num2str(gui.controls.act_x) ':SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                            num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)/MRSCont.processed.(gui.process.StructFields{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm'];
                end
            else
                if  (contains(gui.process.TabTitles{gui.process.Selected},{'A','B','C','D','diff1','diff2','sum'}))
                    StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': SNR(' gui.process.SNR{gui.process.Selected} '): '  num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                                num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)/MRSCont.processed.(gui.process.StructFields{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)*1e6))...
                                ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.freqShift.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.drift.pre.AvgDeltaCr.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)*MRSCont.processed.(gui.process.StructFields{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)/1e6)...
                                ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.drift.post.AvgDeltaCr.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)*MRSCont.processed.(gui.process.StructFields{gui.process.Selected}){gui.controls.Selected}.txfrq(Exp)/1e6) ' Hz'];
                else 
                StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                            num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(gui.process.TabTitles{gui.process.Selected})(Exp,gui.controls.Selected)/MRSCont.processed.(gui.process.StructFields{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm'];
                end
        end                

        set(gui.InfoText.pro, 'String',sprintf(StatText));
        
%%% 3. VISUALIZATION PART OF THIS TAB %%%
        if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            temp = osp_plotProcess(MRSCont, gui.controls.Selected,Selection,SubSpec,Exp); %Create figure
        elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
            temp = osp_plotProcess(MRSCont, gui.controls.Selected,Selection,SubSpec,Exp,gui.controls.act_x); %Create figure
        else
            temp = osp_plotProcess(MRSCont, gui.controls.Selected,Selection,SubSpec,Exp,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]); %Create figure
        end
        %Delete old content
        delete(gui.layout.proDrift.Children.Children)
        delete(gui.layout.proAlgn.Children.Children)
        delete(gui.layout.proPost.Children.Children)
        delete(gui.layout.proPre.Children.Children)
        %Fill window with new content        
        set( temp.Children(1).Children, 'Parent', gui.layout.proDrift.Children ); % Update drift plot
        set( gui.layout.proDrift.Children,'Children',flipud(gui.layout.proDrift.Children.Children));
        set(  gui.layout.proDrift.Children, 'YLim', temp.Children(1).YLim);
        set( flipud(temp.Children(2).Children), 'Parent', gui.layout.proAlgn.Children ); % Update aligned and averaged plot
        set(  gui.layout.proAlgn.Children, 'XLim', temp.Children(2).XLim);
        set(  gui.layout.proAlgn.Children, 'YLim', temp.Children(2).YLim);
        set(  gui.layout.proAlgn.Children, 'XTick', temp.Children(2).XTick);
        set( temp.Children(3).Children, 'Parent', gui.layout.proPost.Children ); % Update post alignment plot
        set(  gui.layout.proPost.Children, 'XLim', temp.Children(3).XLim);
        set(  gui.layout.proPost.Children, 'YLim', temp.Children(3).YLim);
        set(  gui.layout.proPost.Children, 'XTick', temp.Children(3).XTick);
        set( temp.Children(4).Children, 'Parent', gui.layout.proPre.Children ); % Update pre alignment plot
        set(  gui.layout.proPre.Children, 'XLim', temp.Children(4).XLim);
        set(  gui.layout.proPre.Children, 'YLim', temp.Children(4).YLim);
        set(  gui.layout.proPre.Children, 'XTick', temp.Children(4).XTick);
        close( temp );

        % If it is Multivoxel data we have to update the Voxel Position
        % window
        if MRSCont.flags.isMRSI 
            temp = osp_plotRawMRSIpos(MRSCont, 1, [gui.controls.act_y gui.controls.act_x gui.controls.act_z]);
            ViewAxes = gca();
            drawnow
            set( gui.layout.LocPanel.Children,'ColorData', ViewAxes.ColorData );
            close(temp)
        end
        set(gui.upperBox.pro.Info{gui.process.Selected},'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}]);
        set(gui.controls.b_save_proTab{gui.process.Selected},'Callback',{@osp_onPrint,gui});
        setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class 
end