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
        gui.layout.EmptyProPlot = 0;
        if (MRSCont.flags.isUnEdited)
            t = gui.process.Selected;
        end
        if (MRSCont.flags.isMEGA) %Is Edited? Pick the right tab
            if (gui.process.Selected == 1 || gui.process.Selected == 2 || gui.process.Selected == 3)
                t = gui.process.Selected;                   
            else if gui.process.Selected == 4
                    t = 5;                 
                else if MRSCont.flags.hasWater
                    t = 6;
                    else
                    t = 4;
                    end
                end
            end
        end
        if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) %Is Hadamard? Pick the right tab
            if (gui.process.Selected == 1 || gui.process.Selected == 2 || gui.process.Selected == 3 || gui.process.Selected == 4 || gui.process.Selected == 5 || gui.process.Selected == 6)
                t = gui.process.Selected;                   
            else if gui.process.Selected == 7
                    t = 8;                 
                else if MRSCont.flags.hasWater
                    t = 9;
                    else
                    t = 7;
                    end
                end
            end
        end
%%% 2. FILLING INFO PANEL FOR THIS TAB %%%
% All the information from the Raw data is read out here
        if (strcmp(gui.process.Names{t},'A') || strcmp(gui.process.Names{t},'B') || strcmp(gui.process.Names{t},'C') || strcmp(gui.process.Names{t},'D') || strcmp(gui.process.Names{t},'diff1') || strcmp(gui.process.Names{t},'diff2') || strcmp(gui.process.Names{t},'sum'))
            StatText = ['Metabolite Data -> SNR(' gui.process.SNR{t} '): '  num2str(MRSCont.QM.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM: '...
                        num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq/1e6))...
                        ' ppm / Hz \nReference shift: ' num2str(MRSCont.QM.freqShift.(gui.process.Names{t})(gui.controls.Selected)) ' Hz \nAverage Delta F0 Pre Registration: ' num2str(MRSCont.QM.drift.pre.AvgDeltaCr.(gui.process.Names{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq/1e6)...
                        ' Hz; Average Delta F0 Post Registration: ' num2str(MRSCont.QM.drift.post.AvgDeltaCr.(gui.process.Names{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq/1e6) ' Hz'];
        else if strcmp(gui.process.Names{t},'ref')
        StatText = ['Reference Data -> SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM: '...
                    num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)) ' / ' (num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq/1e6))...
                    ' ppm / Hz'];
            else
                StatText = ['Water Data -> SNR(' gui.process.SNR{t} '): ' num2str(MRSCont.QM.SNR.(gui.process.Names{t})(gui.controls.Selected)) '; FWHM: '...
                            num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)) '/' (num2str(MRSCont.QM.FWHM.(gui.process.Names{t})(gui.controls.Selected)*MRSCont.processed.(gui.process.Names{t}){gui.controls.Selected}.txfrq/1e6))...
                            ' ppm / Hz'];
            end
        end
        set(gui.InfoText.pro, 'String',sprintf(StatText))
%%% 3. VISUALIZATION PART OF THIS TAB %%%
        temp = figure( 'Visible', 'off' );
        temp = osp_plotProcess(MRSCont, gui.controls.Selected,gui.process.Names{t},1 ); %Create figure
        %Delete old content
        delete(gui.layout.proDrift.Children.Children)
        delete(gui.layout.proAlgn.Children.Children)
        delete(gui.layout.proPost.Children.Children)
        delete(gui.layout.proPre.Children.Children)
        %Fill window with new content
        set( temp.Children(1).Children, 'Parent', gui.layout.proDrift.Children ); % Update drift plot
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
        set(gui.Info.pro,'Title', ['Actual file: ' MRSCont.files{gui.controls.Selected}] )
        setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class 
end