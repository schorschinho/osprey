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
        gui.upperBox.coreg.Info = gui.layout.(gui.layout.coregTabhandles{gui.coreg.Selected}).Children(2).Children(2);
        gui.InfoText.coreg = gui.layout.(gui.layout.coregTabhandles{gui.coreg.Selected}).Children(2).Children(2).Children;
        % Grid for Plot and Data control sliders
        gui.Plot.coreg = gui.layout.(gui.layout.coregTabhandles{gui.coreg.Selected});
        gui.controls.b_save_coregTab = gui.layout.(gui.layout.coregTabhandles{gui.coreg.Selected}).Children(2).Children(1).Children;
        gui.layout.EmptyProPlot = 0;
%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        StatText = ['Metabolite Data -> Sequence: ' gui.load.Names.Seq '; B0: ' num2str(MRSCont.raw{1,gui.controls.Selected}.Bo) '; TE / TR: ' num2str(MRSCont.raw{1,gui.controls.Selected}.te) ' / ' num2str(MRSCont.raw{1,gui.controls.Selected}.tr) ' ms ' '; spectral bandwidth: ' num2str(MRSCont.raw{1,gui.controls.Selected}.spectralwidth) ' Hz'...
                         '\nraw subspecs: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawSubspecs) '; raw averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.rawAverages) '; averages: ' num2str(MRSCont.raw{1,gui.controls.Selected}.averages)...
                         '; Sz: ' num2str(MRSCont.raw{1,gui.controls.Selected}.sz) ';  dimensions: ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2})) ' x ' num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})) ' mm = '...
                         num2str(MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{1}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{2}) * MRSCont.raw{1,gui.controls.Selected}.geometry.size.(gui.load.Names.Geom{3})/1000) ' ml'];
       set(gui.upperBox.coreg.Info.Children, 'String',sprintf(StatText))
%%% 3. VISUALIZATION PART OF THIS TAB %%%
if gui.coreg.Selected == 1 %Is first structural tab?
        if MRSCont.flags.didSeg && length(gui.Results.coreg.Children) == 2 %Did seg & has been visualized already
            delete( gui.Results.coreg.Children(1) ); %delete seg contents
            delete( gui.Results.coreg.Children(1) ); %delete coreg contents
            temp = figure( 'Visible', 'off' );
            temp = osp_plotCoreg(MRSCont, gui.controls.Selected);
            set( temp.Children, 'Parent', gui.Results.coreg);
            colormap(gui.Results.coreg.Children(1),'gray');
            close( temp );
            temp = figure( 'Visible', 'off' );
            if gui.controls.Selected <= length(MRSCont.seg.tissue.fGM)
                temp = osp_plotSegment(MRSCont, gui.controls.Selected);
                set( temp.Children, 'Parent', gui.Results.coreg );
                colormap(gui.Results.coreg.Children(1),'gray');
                close( temp );
            end
        else
            if MRSCont.flags.didSeg && length(gui.Results.coreg.Children) == 1 %Did seg but has not been visualized yet (Initial run of segement button)
            delete( gui.Results.coreg.Children(1) );
            temp = figure( 'Visible', 'off' );
            temp = osp_plotCoreg(MRSCont, gui.controls.Selected);
            set( temp.Children, 'Parent', gui.Results.coreg);
            colormap(gui.Results.coreg.Children(1),'gray');
            close( temp );
            temp = figure( 'Visible', 'off' );
            temp = osp_plotSegment(MRSCont, gui.controls.Selected);
            set( temp.Children, 'Parent', gui.Results.coreg );
            colormap(gui.Results.coreg.Children(1),'gray');
            close( temp );
            else
                if length(gui.Results.coreg.Children) == 2
                    temp = figure( 'Visible', 'off' );
                    temp = osp_plotCoreg(MRSCont, gui.controls.Selected);
                    delete( gui.Results.coreg.Children(1) );
                    delete( gui.Results.coreg.Children(1) );
                    set( temp.Children, 'Parent', gui.Results.coreg );
                    colormap(gui.Results.coreg.Children,'gray')
                    close( temp );
                elseif length(gui.Results.coreg.Children) == 1 %Only coreg has been performed
                    temp = figure( 'Visible', 'off' );
                    temp = osp_plotCoreg(MRSCont, gui.controls.Selected);
                    delete( gui.Results.coreg.Children(1) );
                    set( temp.Children, 'Parent', gui.Results.coreg );
                    colormap(gui.Results.coreg.Children,'gray')
                    close( temp );
                else % Neither coreg nor segment has been performed
                    temp = figure( 'Visible', 'off' );
                    temp = osp_plotCoreg(MRSCont, gui.controls.Selected);
                    set( temp.Children, 'Parent', gui.Results.coreg );
                    colormap(gui.Results.coreg.Children,'gray')
                    close( temp );
                end
            end
        end
elseif gui.coreg.Selected == 2 %Is second structural tab?
    delete( gui.Results.coreg.Children(1) );
    temp = figure( 'Visible', 'off' );
    temp = osp_plotCoregSecond(MRSCont, gui.controls.Selected);
    set( temp.Children, 'Parent', gui.Results.coreg);
    colormap(gui.Results.coreg.Children(1),'gray');
    close( temp );
    
elseif gui.coreg.Selected == 3 %Is PET tab ?
    
    if length(gui.Results.coreg.Children) > 1
        delete( gui.Results.coreg.Children(1) );
    end
    if length(gui.Results.coreg.Children) > 1
        delete( gui.Results.coreg.Children(1) );
    end
    if length(gui.Results.coreg.Children) > 1
        delete( gui.Results.coreg.Children(1) );
    end
    delete( gui.Results.coreg.Children(1) );
    temp = figure( 'Visible', 'off' );
    temp = osp_plotCoregPET(MRSCont, gui.controls.Selected);
    set( temp.Children(2), 'Parent', gui.Results.coreg );
    set( temp.Children(1), 'Parent', gui.Results.coreg );
    set(gui.Results.coreg,'Heights', [-0.49 -0.49]);
    set(gui.Results.coreg.Children(2), 'Units', 'normalized')
    set(gui.Results.coreg.Children(2), 'OuterPosition', [0,0.5,1,0.5])
    set(gui.Results.coreg.Children(1), 'Units', 'normalized')
    set(gui.Results.coreg.Children(1), 'OuterPosition', [0,0,1,0.5])
    colormap(gui.Results.coreg.Children(2),'gray');
    close( temp );
end
        rmpath(genpath([gui.folder.spmversion filesep])); %Remove SPM path
        set(gui.controls.b_save_coregTab,'Callback',{@osp_onPrint,gui});
        setappdata(gui.figure,'MRSCont',MRSCont); %Write  MRSCont into hidden container in gui class
end