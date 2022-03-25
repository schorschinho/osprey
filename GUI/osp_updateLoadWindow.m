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
        if (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            set(gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children(3).Children(1).Children.Children(4),'String',gui.controls.act_z)
            set(gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children(3).Children(1).Children.Children(5),'String',gui.controls.act_y)
            set(gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children(3).Children(1).Children.Children(6),'String',gui.controls.act_x)
        end
        % Grid for Plot and Data control sliders
         gui.layout.EmptyPlot.data = 0;
%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        switch gui.load.Names.Spec{gui.load.Selected}
            case 'metabolites'
            StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) '\naverages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
            case 'MM'
                    StatText = ['MM Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.spectralwidth) ' Hz'...   %re_mm
                         '\nraw subspecs: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...   %re_mm
                         num2str(MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_mm{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];   %re_mm
            case 'reference'
            StatText = ['Reference Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.spectralwidth) ' Hz'...   %re_mm
                         '\nraw subspecs: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...   %re_mm
                         num2str(MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];   %re_mm
           case 'MM reference' 
                    StatText = ['MM reference Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '\nraw subspecs: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_mm_ref{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
            case 'water'
                StatText = ['Water Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '\nraw subspecs: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.sz) '; dimensions: ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw_w{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
        end

        if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
            set(gui.InfoText.data{gui.load.Selected}, 'String',sprintf(StatText))
        elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
            StatText = ['Voxel ' num2str(gui.controls.act_x) ': ' StatText];
            set(gui.InfoText.data{gui.load.Selected}, 'String',sprintf(StatText))
        else
            StatText = ['Voxel ' num2str(gui.controls.act_x) ' ' num2str(gui.controls.act_y) ' ' num2str(gui.controls.act_z) ': ' StatText];
            set(gui.InfoText.data{gui.load.Selected}, 'String',sprintf(StatText))
        end


%%% 3. VISUALIZATION PART OF THIS TAB %%%
        temp = figure( 'Visible', 'off' );
        Exp = gui.controls.act_Exp;
        switch gui.load.Names.Spec{gui.load.Selected}
            case 'metabolites'
                if Exp > max(MRSCont.opts.MultipleSpectra.metab)
                    Exp = 1;
                    gui.controls.act_Exp = 1;
                end
            if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mets',Exp);
            elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mets',Exp,gui.controls.act_x);
               else
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mets',Exp,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
            end
            if MRSCont.flags.isUnEdited %Is UnEdited?
                ViewAxes = gca();
                delete(gui.Plot.data{gui.load.Selected}.Children(1).Children)
                set(gui.Plot.data{gui.load.Selected}.Children, 'XLim', ViewAxes.XLim)
                set(gui.Plot.data{gui.load.Selected}.Children, 'YLim', ViewAxes.YLim)
                set(ViewAxes.Children, 'Parent', gui.Plot.data{gui.load.Selected}.Children);
                set(gui.Plot.data{gui.load.Selected}.Children.Title, 'String', ViewAxes.Title.String)
            end
            if MRSCont.flags.isMEGA %Is MEGA?
                delete(gui.Plot.data{gui.load.Selected}.Children(1).Children)
                delete(gui.Plot.data{gui.load.Selected}.Children(2).Children)
                set(gui.Plot.data{gui.load.Selected}.Children(2), 'XLim', temp.Children(2).XLim)
                set(gui.Plot.data{gui.load.Selected}.Children(1), 'XLim', temp.Children(1).XLim)
                set(gui.Plot.data{gui.load.Selected}.Children(2), 'YLim', temp.Children(2).YLim)
                set(gui.Plot.data{gui.load.Selected}.Children(1), 'YLim', temp.Children(1).YLim)
                set(temp.Children(2).Children, 'Parent', gui.Plot.data{gui.load.Selected}.Children(2));
                set(temp.Children(1).Children, 'Parent', gui.Plot.data{gui.load.Selected}.Children(1));
                set(gui.Plot.data{gui.load.Selected}.Children(2).Title, 'String', temp.Children(2).Title.String)
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
            set(gui.upperBox.data.Info{gui.load.Selected},'Title', ['Actual file: ' MRSCont.files{Exp,gui.controls.Selected}] );
            case 'MM'
                if Exp > max(MRSCont.opts.MultipleSpectra.mm)
                    Exp = 1;
                    gui.controls.act_Exp = 1;
                end
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm',Exp);
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm',Exp,gui.controls.act_x);
                          else
                            temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm',Exp,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                end
                if MRSCont.flags.isUnEdited %Is UnEdited?  
                    ViewAxes = gca();
                    delete(gui.Plot.data{gui.load.Selected}.Children(1).Children)
                    set(gui.Plot.data{gui.load.Selected}.Children, 'XLim', ViewAxes.XLim)
                    set(gui.Plot.data{gui.load.Selected}.Children, 'YLim', ViewAxes.YLim)
                    set(ViewAxes.Children, 'Parent', gui.Plot.data{gui.load.Selected}.Children);
                    set(gui.Plot.data{gui.load.Selected}.Children.Title, 'String', ViewAxes.Title.String)
                end
                if MRSCont.flags.isMEGA %Is MEGA?
                    delete(gui.Plot.data{gui.load.Selected}.Children(1).Children)
                    delete(gui.Plot.data{gui.load.Selected}.Children(2).Children)
                    set(gui.Plot.data{gui.load.Selected}.Children(2), 'XLim', temp.Children(2).XLim)
                    set(gui.Plot.data{gui.load.Selected}.Children(1), 'XLim', temp.Children(1).XLim)
                    set(gui.Plot.data{gui.load.Selected}.Children(2), 'YLim', temp.Children(2).YLim)
                    set(gui.Plot.data{gui.load.Selected}.Children(1), 'YLim', temp.Children(1).YLim)
                    set(temp.Children(2).Children, 'Parent', gui.Plot.data{gui.load.Selected}.Children(2));
                    set(temp.Children(1).Children, 'Parent', gui.Plot.data{gui.load.Selected}.Children(1));
                    set(gui.Plot.data{gui.load.Selected}.Children(2).Title, 'String', temp.Children(2).Title.String)
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
            set(gui.upperBox.data.Info{gui.load.Selected},'Title', ['Actual file: ' MRSCont.files_mm{Exp,gui.controls.Selected}] );
            case 'reference'
                if Exp > max(MRSCont.opts.MultipleSpectra.ref)
                    Exp = 1;
                    gui.controls.act_Exp = 1;
                end
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref',Exp);
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref',Exp,gui.controls.act_x);
                      else
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'ref',Exp,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                end
                ViewAxes = gca();
                delete(gui.Plot.data{gui.load.Selected}.Children.Children)
                set(ViewAxes.Children, 'Parent', gui.Plot.data{gui.load.Selected}.Children);
                set(gui.Plot.data{gui.load.Selected}.Children.Title, 'String', ViewAxes.Title.String);
                set(gui.Plot.data{gui.load.Selected}.Children, 'XLim',ViewAxes.XLim);
                set(gui.upperBox.data.Info{gui.load.Selected},'Title', ['Actual file: ' MRSCont.files_ref{Exp,gui.controls.Selected}] );
            case 'MM reference'
                if Exp > max(MRSCont.opts.MultipleSpectra.mm_ref)
                    Exp = 1;
                    gui.controls.act_Exp = 1;
                end
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm_ref',Exp);
                    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm_ref',Exp,gui.controls.act_x);
                      else
                        temp = osp_plotLoad(MRSCont, gui.controls.Selected,'mm_ref',Exp,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                end
                ViewAxes = gca();
                delete(gui.Plot.data{gui.load.Selected}.Children.Children)
                set(ViewAxes.Children, 'Parent', gui.Plot.data{gui.load.Selected}.Children);
                set(gui.Plot.data{gui.load.Selected}.Children.Title, 'String', ViewAxes.Title.String);
                set(gui.Plot.data{gui.load.Selected}.Children, 'XLim',ViewAxes.XLim);
                set(gui.upperBox.data.Info{gui.load.Selected},'Title', ['Actual file: ' MRSCont.files_mm_ref{Exp,gui.controls.Selected}] );
            case 'water'
                if Exp > max(MRSCont.opts.MultipleSpectra.w)
                    Exp = 1;
                    gui.controls.act_Exp = 1;
                end
                if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',Exp);
                elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',Exp,gui.controls.act_x);
              else
                temp = osp_plotLoad(MRSCont, gui.controls.Selected,'w',Exp,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]);
                end
                ViewAxes = gca();
                delete(gui.Plot.data{gui.load.Selected}.Children(1).Children)
                set(ViewAxes.Children, 'Parent', gui.Plot.data{gui.load.Selected}.Children(1));
                set(  gui.Plot.data{gui.load.Selected}.Children(1).Title, 'String',ViewAxes.Title.String)
                set(  gui.Plot.data{gui.load.Selected}.Children(1), 'XLim',ViewAxes.XLim)   
                set(gui.upperBox.data.Info{gui.load.Selected},'Title', ['Actual file: ' MRSCont.files_w{Exp,gui.controls.Selected}] );
        end

            

        % Get rid of the Load figure
        close( temp );
        h = findall(groot,'Type','figure');
        for ff = 1 : length(h)
            if ~(strcmp(h(ff).Tag, 'Osprey') ||  strcmp(h(ff).Tag, 'TMWWaitbar'))
                close(h(ff))
            end
        end

        % If it is Multivoxel data we have to update the Voxel Position
        % window
        if MRSCont.flags.isMRSI
            temp = osp_plotRawMRSIpos(MRSCont, 1, [gui.controls.act_y gui.controls.act_x gui.controls.act_z]);
            ViewAxes = gca();
            drawnow
            set( gui.layout.LocPanel.Children,'ColorData', ViewAxes.ColorData );
            close(temp)
        end        
        set(gui.controls.b_save_RawTab{gui.load.Selected},'Callback',{@osp_onPrint,gui});
        if MRSCont.nDatasets(2) > 1
             set(gui.layout.(gui.layout.rawTabhandles{gui.load.Selected}).Children(2).Children(3).Children(1).Children.Children(2),'String',gui.controls.act_Exp)
         end
        setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end
