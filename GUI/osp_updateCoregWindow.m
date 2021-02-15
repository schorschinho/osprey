function osp_updateCoregWindow(gui)
%% osp_updateCoregWindow
%   This function updates the coreg tab.
%
%
%   USAGE:
%       osp_updateCoregWindow(gui);
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
        addpath(genpath([gui.folder.spmversion filesep])); % Add SPM  path
        gui.controls.b_save_coregTab = gui.layout.coregTab.Children(2).Children(1).Children;
         if (isfield(MRSCont.flags, 'isPRIAM') || isfield(MRSCont.flags, 'isMRSI')) &&  (MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            set(gui.layout.coregTab.Children(2).Children(4).Children(1).Children.Children(4),'String',gui.controls.act_z)
            set(gui.layout.coregTab.Children(2).Children(4).Children(1).Children.Children(5),'String',gui.controls.act_y)
            set(gui.layout.coregTab.Children(2).Children(4).Children(1).Children.Children(6),'String',gui.controls.act_x)
         end
        gui.layout.EmptyProPlot = 0;
%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '\nraw subspecs: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) ';  dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
       set(gui.upperBox.coreg.Info.Children, 'String',sprintf(StatText))
%%% 3. VISUALIZATION PART OF THIS TAB %%%
        if MRSCont.flags.didSeg && length(gui.Results.coreg.Children) == 2 %Did seg & has been visualized already
            delete( gui.Results.coreg.Children(1) ); %delete seg contents
            delete( gui.Results.coreg.Children(1) ); %delete coreg contents
            temp = figure( 'Visible', 'off' );
            if ~isfield(MRSCont.flags,'isPRIAM')  && ~MRSCont.flags.isPRIAM
                temp = osp_plotCoreg(MRSCont, gui.controls.Selected);
            else
                temp = osp_plotCoreg(MRSCont, gui.controls.Selected,gui.controls.act_x);
            end               
            set( temp.Children, 'Parent', gui.Results.coreg);
            colormap(gui.Results.coreg.Children(1),'gray');
            close( temp );
            temp = figure( 'Visible', 'off' );
            if gui.controls.Selected <= length(MRSCont.seg.tissue.fGM)
                if ~isfield(MRSCont.flags,'isPRIAM')  && ~MRSCont.flags.isPRIAM
                    temp = osp_plotSegment(MRSCont, gui.controls.Selected);
                else
                    temp = osp_plotSegment(MRSCont, gui.controls.Selected,gui.controls.act_x);
                end
                set( temp.Children, 'Parent', gui.Results.coreg );
                colormap(gui.Results.coreg.Children(1),'gray');
                close( temp );
            end
        else if MRSCont.flags.didSeg && length(gui.Results.coreg.Children) == 1 %Did seg but has not been visualized yet (Initial run of segement button)
            delete( gui.Results.coreg.Children(1) );
            temp = figure( 'Visible', 'off' );
            if ~isfield(MRSCont.flags,'isPRIAM')  && ~MRSCont.flags.isPRIAM
                temp = osp_plotCoreg(MRSCont, gui.controls.Selected);
            else
                temp = osp_plotCoreg(MRSCont, gui.controls.Selected,gui.controls.act_x);
            end
            set( temp.Children, 'Parent', gui.Results.coreg);
            colormap(gui.Results.coreg.Children(1),'gray');
            close( temp );
            temp = figure( 'Visible', 'off' );
            if ~isfield(MRSCont.flags,'isPRIAM')  && ~MRSCont.flags.isPRIAM
                temp = osp_plotSegment(MRSCont, gui.controls.Selected);
            else
                temp = osp_plotSegment(MRSCont, gui.controls.Selected,gui.controls.act_x);
            end
            set( temp.Children, 'Parent', gui.Results.coreg );
            colormap(gui.Results.coreg.Children(1),'gray');
            close( temp );
            else if length(gui.Results.coreg.Children) == 1 %Only coreg has been performed
                    temp = figure( 'Visible', 'off' );
                    if ~isfield(MRSCont.flags,'isPRIAM')  && ~MRSCont.flags.isPRIAM
                        temp = osp_plotCoreg(MRSCont, gui.controls.Selected);
                    else
                        temp = osp_plotCoreg(MRSCont, gui.controls.Selected,gui.controls.act_x);
                    end
                    delete( gui.Results.coreg.Children(1) );
                    set( temp.Children, 'Parent', gui.Results.coreg );
                    colormap(gui.Results.coreg.Children,'gray')
                    close( temp );
                else % Neither coreg nor segment has been performed
                    temp = figure( 'Visible', 'off' );
                    if ~isfield(MRSCont.flags,'isPRIAM')  && ~MRSCont.flags.isPRIAM
                        temp = osp_plotCoreg(MRSCont, gui.controls.Selected);
                    else
                        temp = osp_plotCoreg(MRSCont, gui.controls.Selected,gui.controls.act_x);
                    end
                    set( temp.Children, 'Parent', gui.Results.coreg );
                    colormap(gui.Results.coreg.Children,'gray')
                    close( temp ); 
                end
            end
        end

        rmpath(genpath([gui.folder.spmversion filesep])); %Remove SPM path
        set(gui.controls.b_save_coregTab,'Callback',{@osp_onPrint,gui});
        setappdata(gui.figure,'MRSCont',MRSCont); %Write  MRSCont into hidden container in gui class
end