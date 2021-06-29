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
        gui.controls.b_save_proTab = gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children(1).Children;
        gui.upperBox.pro.Info = gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children(2);
        if (isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) && (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            set(gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children(3).Children(1).Children.Children(4),'String',gui.controls.act_z);
            set(gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children(3).Children(1).Children.Children(5),'String',gui.controls.act_y);
            set(gui.layout.(gui.layout.proTabhandles{gui.process.Selected}).Children(2).Children(3).Children(1).Children.Children(6),'String',gui.controls.act_x);
        end
        gui.layout.EmptyProPlot = 0;
        Selection = gui.process.Names{gui.process.Selected};
%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            
            if ~strcmp (Selection, 'ref') && ~strcmp (Selection, 'w') && ~strcmp (Selection, 'mm') %Metabolite data?
                StatText = ['Metabolite Data -> SNR(' gui.process.SNR{gui.process.Selected} '): '  num2str(MRSCont.QM.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                            num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM.freqShift.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq/1e6)...
                            ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq/1e6) ' Hz'];
            else if strcmp (Selection, 'ref') %Reference data?
            StatText = ['Reference Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                        num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                        ' Hz / ppm'];
                else if ~strcmp (Selection, 'mm') % re_mm
                    StatText = ['Water Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                                num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '/' (num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm'];
                    else
                    StatText = ['MM Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                                num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '/' (num2str(MRSCont.QM.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm'];
                    end
                end
            end
        elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
          if ~strcmp (Selection, 'ref') && ~strcmp (Selection, 'w') && ~strcmp (Selection, 'mm') %Metabolite data?
                StatText = ['Voxel ' num2str(gui.controls.act_x) ': Metabolite Data -> SNR(' gui.process.SNR{gui.process.Selected} '): '  num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                            num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM{1,gui.controls.act_x}.freqShift.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM{1,gui.controls.act_x}.drift.pre.AvgDeltaCr.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq/1e6)...
                            ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM{1,gui.controls.act_x}.drift.post.AvgDeltaCr.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq/1e6) ' Hz'];
            else if strcmp (Selection, 'ref') %Reference data?
            StatText = ['Voxel ' num2str(gui.controls.act_x) ': Reference Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                        num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                        ' Hz / ppm'];
                else if ~strcmp (Selection, 'mm') % re_mm
                    StatText = ['Voxel ' num2str(gui.controls.act_x) ': Water Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                                num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '/' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm'];
                    else
                    StatText = ['Voxel ' num2str(gui.controls.act_x) ': MM Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM{1,gui.controls.act_x}.SNR.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                                num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)) '/' (num2str(MRSCont.QM{1,gui.controls.act_x}.FWHM.(gui.process.Names{gui.process.Selected})(gui.controls.Selected)/MRSCont.processed.(gui.process.Names{gui.process.Selected}){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm'];
                    end
                end
            end
        else
            if (strcmp(Selection,'A') || strcmp(Selection,'B') || strcmp(Selection,'C') || strcmp(Selection,'D') || strcmp(Selection,'diff1') || strcmp(Selection,'diff2') || strcmp(Selection,'sum'))
                StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': Metabolite Data -> SNR(' gui.process.SNR{gui.process.Selected} '): '  num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(Selection)(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                            num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(Selection)(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(Selection)(gui.controls.Selected)/MRSCont.processed.(Selection){gui.controls.Selected}.txfrq*1e6))...
                            ' Hz / ppm \nReference shift: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.freqShift.(Selection)(gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.drift.pre.AvgDeltaCr.(Selection)(gui.controls.Selected)*MRSCont.processed.(Selection){gui.controls.Selected}.txfrq/1e6)...
                            ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.drift.post.AvgDeltaCr.(Selection)(gui.controls.Selected)*MRSCont.processed.(Selection){gui.controls.Selected}.txfrq/1e6) ' Hz'];
            else if strcmp(Selection,'ref')
            StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': Reference Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(Selection)(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                        num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(Selection)(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(Selection)(gui.controls.Selected)/MRSCont.processed.(Selection){gui.controls.Selected}.txfrq*1e6))...
                        ' Hz / ppm'];
                else
                    if ~strcmp(Selection,'mm') %re
                    StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': Water Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(Selection)(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                                num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(Selection)(gui.controls.Selected)) '/' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(Selection)(gui.controls.Selected)/MRSCont.processed.(Selection){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm'];
                    else
                        StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ': MM Data -> SNR(' gui.process.SNR{gui.process.Selected} '): ' num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.SNR.(Selection)(gui.controls.Selected)) '; FWHM (' gui.process.SNR{gui.process.Selected} '): '...
                                num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(Selection)(gui.controls.Selected)) '/' (num2str(MRSCont.QM{gui.controls.act_x,gui.controls.act_y}.FWHM.(Selection)(gui.controls.Selected)/MRSCont.processed.(Selection){gui.controls.Selected}.txfrq*1e6))...
                                ' Hz / ppm'];
                    end %re
                end
            end      
        end

        set(gui.upperBox.pro.Info.Children, 'String',sprintf(StatText));
        
%%% 3. VISUALIZATION PART OF THIS TAB %%%
        if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            temp = osp_plotProcess(MRSCont, gui.controls.Selected,Selection); %Create figure
        elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
            temp = osp_plotProcess(MRSCont, gui.controls.Selected,Selection,gui.controls.act_x); %Create figure
        else
            temp = osp_plotProcess(MRSCont, gui.controls.Selected,Selection,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]); %Create figure
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
        set( temp.Children(2).Children, 'Parent', gui.layout.proAlgn.Children ); % Update aligned and averaged plot
        set(  gui.layout.proAlgn.Children, 'XLim', temp.Children(2).XLim);
        set(  gui.layout.proAlgn.Children, 'YLim', temp.Children(2).YLim);
        set( temp.Children(3).Children, 'Parent', gui.layout.proPost.Children ); % Update post alignment plot
        set(  gui.layout.proPost.Children, 'XLim', temp.Children(3).XLim);
        set(  gui.layout.proPost.Children, 'YLim', temp.Children(3).YLim);
        set( temp.Children(4).Children, 'Parent', gui.layout.proPre.Children ); % Update pre alignment plot
        set(  gui.layout.proPre.Children, 'XLim', temp.Children(4).XLim);
        set(  gui.layout.proPre.Children, 'YLim', temp.Children(4).YLim);
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
        set(gui.upperBox.pro.Info,'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}]);
        set(gui.controls.b_save_proTab,'Callback',{@osp_onPrint,gui});
        setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class 
end