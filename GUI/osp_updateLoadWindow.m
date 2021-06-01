function osp_updateLoadWindow(gui)
%% osp_updateLoadWindow
%   This function updates the load tab.
%
%
%   USAGE:
%       osp_updateLoadWindow(gui);
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
        gui.upperBox.data.Info = gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children(2);
        if (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            set(gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children(3).Children(1).Children.Children(4),'String',gui.controls.act_z)
            set(gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children(3).Children(1).Children.Children(5),'String',gui.controls.act_y)
            set(gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children(3).Children(1).Children.Children(6),'String',gui.controls.act_x)
        end
        gui.InfoText.data = gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children(2).Children;
        % Grid for Plot and Data control sliders
        gui.Plot.data = gui.layout.(gui.layout.rawTabhandles{gui.load.Selected});
        gui.controls.b_save_RawTab = gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children(1).Children;
        gui.layout.EmptyPlot.data = 0;
%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        if gui.load.Selected == 1 %Is metabolite data?
            StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) '\naverages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
        else
            if MRSCont.flags.hasMM   %re_mm
               if gui.load.Selected == 2 %re_mm
                    StatText = ['MM Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.spectralwidth) ' Hz'...   %re_mm
                         '\nraw subspecs: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...   %re_mm
                         num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];   %re_mm
            end  %re_mm
           if gui.load.Selected == 3 %Is water or ref data?    %re_mm
            StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.spectralwidth) ' Hz'...   %re_mm
                         '\nraw subspecs: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...   %re_mm
                         num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];   %re_mm
            end   %re_mm
            else    %re_mm
            if gui.load.Selected == 2 %Is water or ref data?
                if MRSCont.flags.hasRef
                    StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                                 '\nraw subspecs: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.averages)...
                                 '; Sz: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                                 num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
                elseif MRSCont.flags.hasWater
                    StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '\nraw subspecs: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
                end
            else
                StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '\nraw subspecs: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
            end
            end    %re_mm

        end
        if ~isfield(MRSCont.flags,'isPRIAM') && ~isfield(MRSCont.flags,'isMRSI') && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
            set(gui.InfoText.data, 'String',sprintf(StatText))
        elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
            StatText = ['Voxel ' num2str(gui.controls.act_x) ': ' StatText];
            set(gui.InfoText.data, 'String',sprintf(StatText))
        else
            StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ' ' num2str(gui.controls.act_z) ': ' StatText];
            set(gui.InfoText.data, 'String',sprintf(StatText))
        end


%%% 3. VISUALIZATION PART OF THIS TAB %%%
        temp = figure( 'Visible', 'off' );
        if gui.load.Selected == 1 %Is Metabolite data/tab?
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mets');
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mets',gui.controls.act_x);
               else
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mets',[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
            end
            if MRSCont.flags.isUnEdited %Is UnEdited?
                ViewAxes = gca();
                delete(gui.Plot.data.Children(1).Children(1).Children)
                set(gui.Plot.data.Children(1).Children(1), 'XLim', ViewAxes.XLim)
                set(gui.Plot.data.Children(1).Children(1), 'YLim', ViewAxes.YLim)
                set(ViewAxes.Children, 'Parent', gui.Plot.data.Children(1).Children(1));
                set(gui.Plot.data.Children(1).Children(1).Title, 'String', ViewAxes.Title.String)
            end
            if MRSCont.flags.isMEGA %Is MEGA?
                delete(gui.Plot.data.Children(1).Children(1).Children)
                delete(gui.Plot.data.Children(1).Children(2).Children)
                set(gui.Plot.data.Children(1).Children(2), 'XLim', temp.Children(2).XLim)
                set(gui.Plot.data.Children(1).Children(1), 'XLim', temp.Children(1).XLim)
                set(gui.Plot.data.Children(1).Children(2), 'YLim', temp.Children(2).YLim)
                set(gui.Plot.data.Children(1).Children(1), 'YLim', temp.Children(1).YLim)
                set(temp.Children(2).Children, 'Parent', gui.Plot.data.Children(1).Children(2));
                set(temp.Children(1).Children, 'Parent', gui.Plot.data.Children(1).Children(1));
                set(gui.Plot.data.Children(1).Children(2).Title, 'String', temp.Children(2).Title.String)
            end
            if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) % Is HERMES/HERCULES
               try
                    delete(gui.layout.multiAload.Children.Children)
                    delete(gui.layout.multiBload.Children.Children)
                    delete(gui.layout.multiCload.Children.Children)
                    delete(gui.layout.multiDload.Children.Children)
               catch
               end
                %Fill window with new content
                set(  gui.layout.multiDload.Children, 'XLim', temp.Children(1).XLim);
                set(  gui.layout.multiDload.Children, 'YLim', temp.Children(1).YLim);
                set( temp.Children(1).Children, 'Parent', gui.layout.multiDload.Children ); % Update drift plot
                set(  gui.layout.multiCload.Children, 'XLim', temp.Children(2).XLim);
                set(  gui.layout.multiCload.Children, 'YLim', temp.Children(2).YLim);
                set( temp.Children(2).Children, 'Parent', gui.layout.multiCload.Children ); % Update aligned and averaged plot
                set(  gui.layout.multiBload.Children, 'XLim', temp.Children(3).XLim);
                set(  gui.layout.multiBload.Children, 'YLim', temp.Children(3).YLim);
                set( temp.Children(3).Children, 'Parent', gui.layout.multiBload.Children ); % Update post alignment plot
                set(  gui.layout.multiAload.Children, 'XLim', temp.Children(4).XLim);
                set(  gui.layout.multiAload.Children, 'YLim', temp.Children(4).YLim);
                set( temp.Children(4).Children, 'Parent', gui.layout.multiAload.Children ); % Update pre alignment plot
            end
        else
            if MRSCont.flags.hasMM %re_mm
                if gui.load.Selected == 2 %ref data/tab  %re_mm
                    if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm');
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm',gui.controls.act_x);
                          else
                            temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm',[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                    end
                ViewAxes = gca();
                delete(gui.Plot.data.Children(1).Children(1).Children)
                set(  gui.Plot.data.Children(1).Children(1), 'XLim',ViewAxes.XLim)
                set(ViewAxes.Children, 'Parent', gui.Plot.data.Children(1).Children(1));
                set(  gui.Plot.data.Children(1).Children(1).Title, 'String',ViewAxes.Title.String)
                end %re_mm
                if gui.load.Selected == 3 %ref data/tab %re_mm
                    if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref');
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref',gui.controls.act_x);
                      else
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref',[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                    end
                ViewAxes = gca();
                delete(gui.Plot.data.Children(1).Children(1).Children)
                set(ViewAxes.Children, 'Parent', gui.Plot.data.Children(1).Children(1));
                set(gui.Plot.data.Children(1).Children(1).Title, 'String', ViewAxes.Title.String);
                set(gui.Plot.data.Children(1).Children(1), 'XLim',ViewAxes.XLim);
                end %re_mm
                if gui.load.Selected == 4 %ref data/tab %re_mm
                    if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w');
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',gui.controls.act_x);
                    else
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                    end
                ViewAxes = gca();
                delete(gui.Plot.data.Children(1).Children(1).Children)
                set(ViewAxes.Children, 'Parent', gui.Plot.data.Children(1).Children(1));
                set(  gui.Plot.data.Children(1).Children(1).Title, 'String',ViewAxes.Title.String)
                set(  gui.Plot.data.Children(1).Children(1), 'XLim',ViewAxes.XLim)
                end
                        else %re_mm
            if gui.load.Selected == 2 %Is Ref data/tab?
                if MRSCont.flags.hasRef
                    if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref');
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref',gui.controls.act_x);
                    end
                    ViewAxes = gca();
                    delete(gui.Plot.data.Children(1).Children(1).Children)
                    set(ViewAxes.Children, 'Parent', gui.Plot.data.Children(1).Children(1));
                    set(gui.Plot.data.Children(1).Children(1).Title, 'String', ViewAxes.Title.String)
                    set(gui.Plot.data.Children(1).Children(1), 'XLim',ViewAxes.XLim)
                elseif MRSCont.flags.hasWater
                    if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w');
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',gui.controls.act_x);
                    else
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                    end
                    ViewAxes = gca();
                    delete(gui.Plot.data.Children(1).Children(1).Children)
                    set(ViewAxes.Children, 'Parent', gui.Plot.data.Children(1).Children(1));
                    set(  gui.Plot.data.Children(1).Children(1).Title, 'String',ViewAxes.Title.String)
                    set(  gui.Plot.data.Children(1).Children(1), 'XLim',ViewAxes.XLim)
                end
            else %Is water data/tab?
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w');
                elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',gui.controls.act_x);
              else
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                end
                ViewAxes = gca();
                delete(gui.Plot.data.Children(1).Children(1).Children)
                set(ViewAxes.Children, 'Parent', gui.Plot.data.Children(1).Children(1));
                set(  gui.Plot.data.Children(1).Children(1).Title, 'String',ViewAxes.Title.String)
                set(  gui.Plot.data.Children(1).Children(1), 'XLim',ViewAxes.XLim)
            end
                        end
        end
        % Get rid of the Load figure
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
        set(gui.upperBox.data.Info,'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}] );
        set(gui.controls.b_save_RawTab,'Callback',{@osp_onPrint,gui});
        setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end
