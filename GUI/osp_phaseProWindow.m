function osp_phaseProWindow(gui)
%% osp_phaseProWindow
%   This function updates the process tab for manual phasing.
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
  
%%% 2. PHASING AND FREQUENCY SHIFT HAPPENS HERE %%%   
    if isfield(MRSCont.processed.metab{gui.controls.Selected}, 'extras')
        procData = op_takeextra(MRSCont.processed.metab{gui.controls.Selected},Exp);
        if gui.controls.posPhase0Shift
            procData=op_addphase(procData ,5);
            if isfield(procData, 'manual') && isfield(procData.manual,'ph0') && Exp <= length(procData.manual.ph0)
                procData.manual.ph0(Exp) = procData.manual.ph0(Exp) + 5;
            else
                procData.manual.ph0(Exp) =  5;
            end
        end
        if gui.controls.negPhase0Shift
            procData=op_addphase(procData ,-5);
            if isfield(procData, 'manual') && isfield(procData.manual,'ph0') && Exp <= length(procData.manual.ph0)
                procData.manual.ph0(Exp) = procData.manual.ph0(Exp) - 5;
            else
                procData.manual.ph0(Exp) =  -5;
            end
        end
        if gui.controls.posPhase1Shift
            procData=op_addphase(procData ,0,0.00001);
            if isfield(procData, 'manual') && isfield(procData.manual,'ph1') && Exp <= length(procData.manual.ph1)
                procData.manual.ph1(Exp) = procData.manual.ph1(Exp) + 0.00001;
            else
                procData.manual.ph1(Exp) =  0.00001;
            end
        end
        if gui.controls.negPhase1Shift
            procData=op_addphase(procData ,0,-0.00001);
            if isfield(procData, 'manual') && isfield(procData.manual,'ph1') && Exp <= length(procData.manual.ph1)
                procData.manual.ph1(Exp) = procData.manual.ph1(Exp) - 0.00001;
            else
                procData.manual.ph1(Exp) =  -0.00001;
            end
        end
        if gui.controls.flipPhase
            procData=op_addphase(procData ,180);
            if isfield(procData, 'manual') && isfield(procData.manual,'ph0') && Exp <= length(procData.manual.ph0)
                procData.manual.ph0(Exp) = procData.manual.ph0(Exp) + 180;
            else
                procData.manual.ph0(Exp) =  180;
            end
        end
        if gui.controls.posFreqShift
            procData=op_freqshift(procData ,1);
            if isfield(procData, 'manual') && isfield(procData.manual,'f') && Exp <= length(procData.manual.f)
                procData.manual.f(Exp) = procData.manual.f(Exp) + 1;
            else
                procData.manual.f(Exp) =  1;
            end
        end
        if gui.controls.negFreqShift
            procData=op_freqshift(procData ,-1);
            if isfield(procData, 'manual') && isfield(procData.manual,'f') && Exp <= length(procData.manual.f)
                procData.manual.f(Exp) = procData.manual.f(Exp) - 1;
            else
                procData.manual.f(Exp) =  -1;
            end
        end
        if isfield(procData, 'manual') && isfield(procData.manual,'ph0')
            if abs(procData.manual.ph0(Exp)) > 360
                procData.manual.ph0(Exp) = procData.manual.ph0(Exp) + sign(procData.manual.ph0(Exp))*(-1)*360;
            end
        end
        MRSCont.processed.metab{gui.controls.Selected}=op_addextra(MRSCont.processed.metab{gui.controls.Selected},procData,Exp); 
    
        if isfield(MRSCont, 'processed_no_align')
            procData = op_takeextra(MRSCont.processed_no_align.metab{gui.controls.Selected},Exp);
            if gui.controls.posPhase0Shift
                procData=op_addphase(procData ,5);
                if isfield(procData, 'manual') && isfield(procData.manual,'ph0') && Exp <= length(procData.manual.ph0)
                    procData.manual.ph0(Exp) = procData.manual.ph0(Exp) + 5;
                else
                    procData.manual.ph0(Exp) =  5;
                end
            end
            if gui.controls.negPhase0Shift
                procData=op_addphase(procData ,-5);
                if isfield(procData, 'manual') && isfield(procData.manual,'ph0') && Exp <= length(procData.manual.ph0)
                    procData.manual.ph0(Exp) = procData.manual.ph0(Exp) - 5;
                else
                    procData.manual.ph0(Exp) =  -5;
                end
            end
            if gui.controls.posPhase1Shift
                procData=op_addphase(procData ,0,0.00001);
                if isfield(procData, 'manual') && isfield(procData.manual,'ph1') && Exp <= length(procData.manual.ph1)
                    procData.manual.ph1(Exp) = procData.manual.ph1(Exp) + 0.00001;
                else
                    procData.manual.ph1(Exp) =  0.00001;
                end
            end
            if gui.controls.negPhase1Shift
                procData=op_addphase(procData ,0,-0.00001);
                if isfield(procData, 'manual') && isfield(procData.manual,'ph1') && Exp <= length(procData.manual.ph1)
                    procData.manual.ph1(Exp) = procData.manual.ph1(Exp) - 0.00001;
                else
                    procData.manual.ph1(Exp) =  -0.00001;
                end
            end
            if gui.controls.flipPhase
                procData=op_addphase(procData ,180);
                if isfield(procData, 'manual') && isfield(procData.manual,'ph0') && Exp <= length(procData.manual.ph0)
                    procData.manual.ph0(Exp) = procData.manual.ph0(Exp) + 180;
                else
                    procData.manual.ph0(Exp) =  180;
                end
            end
            if gui.controls.posFreqShift
                procData=op_freqshift(procData ,1);
                if isfield(procData, 'manual') && isfield(procData.manual,'f') && Exp <= length(procData.manual.f)
                    procData.manual.f(Exp) = procData.manual.f(Exp) + 1;
                else
                    procData.manual.f(Exp) =  1;
                end
            end
            if gui.controls.negFreqShift
                procData=op_freqshift(procData ,-1);
                if isfield(procData, 'manual') && isfield(procData.manual,'f') && Exp <= length(procData.manual.f)
                    procData.manual.f(Exp) = procData.manual.f(Exp) - 1;
                else
                    procData.manual.f(Exp) =  -1;
                end
            end
            if isfield(procData, 'manual') && isfield(procData.manual,'ph0')
                if abs(procData.manual.ph0(Exp)) > 360
                    procData.manual.ph0(Exp) = procData.manual.ph0(Exp) + sign(procData.manual.ph0(Exp))*(-1)*360;
                end
            end
            MRSCont.processed_no_align.metab{gui.controls.Selected}=op_addextra(MRSCont.processed_no_align,procData,Exp);
        end
    else
        if gui.controls.posPhase0Shift
            MRSCont.processed.metab{gui.controls.Selected}=op_addphase(MRSCont.processed.metab{gui.controls.Selected} ,5);
            if isfield(MRSCont.processed.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed.metab{gui.controls.Selected}.manual,'ph0')
                MRSCont.processed.metab{gui.controls.Selected}.manual.ph0 = MRSCont.processed.metab{gui.controls.Selected}.manual.ph0 + 5;
            else
                MRSCont.processed.metab{gui.controls.Selected}.manual.ph0 =  5;
            end
        end
        if gui.controls.negPhase0Shift
            MRSCont.processed.metab{gui.controls.Selected}=op_addphase(MRSCont.processed.metab{gui.controls.Selected} ,-5);
            if isfield(MRSCont.processed.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed.metab{gui.controls.Selected}.manual,'ph0')
                MRSCont.processed.metab{gui.controls.Selected}.manual.ph0 = MRSCont.processed.metab{gui.controls.Selected}.manual.ph0 - 5;
            else
                MRSCont.processed.metab{gui.controls.Selected}.manual.ph0 =  -5;
            end
        end
        if gui.controls.posPhase1Shift
            MRSCont.processed.metab{gui.controls.Selected}=op_addphase(MRSCont.processed.metab{gui.controls.Selected} ,0,0.00001);
            if isfield(MRSCont.processed.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed.metab{gui.controls.Selected}.manual,'ph1')
                MRSCont.processed.metab{gui.controls.Selected}.manual.ph1 = MRSCont.processed.metab{gui.controls.Selected}.manual.ph1 + 0.00001;
            else
                MRSCont.processed.metab{gui.controls.Selected}.manual.ph1 =  0.00001;
            end
        end
        if gui.controls.negPhase1Shift
            MRSCont.processed.metab{gui.controls.Selected}=op_addphase(MRSCont.processed.metab{gui.controls.Selected} ,0,-0.00001);
            if isfield(MRSCont.processed.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed.metab{gui.controls.Selected}.manual,'ph1')
                MRSCont.processed.metab{gui.controls.Selected}.manual.ph1 = MRSCont.processed.metab{gui.controls.Selected}.manual.ph1 -0.00001;
            else
                MRSCont.processed.metab{gui.controls.Selected}.manual.ph1 =  -0.00001;
            end
        end
        if gui.controls.flipPhase
            MRSCont.processed.metab{gui.controls.Selected}=op_addphase(MRSCont.processed.metab{gui.controls.Selected} ,180);
            if isfield(MRSCont.processed.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed.metab{gui.controls.Selected}.manual,'ph0')
                MRSCont.processed.metab{gui.controls.Selected}.manual.ph0 = MRSCont.processed.metab{gui.controls.Selected}.manual.ph0 + 180;
            else
                MRSCont.processed.metab{gui.controls.Selected}.manual.ph0 =  180;
            end
        end
        if gui.controls.posFreqShift
            MRSCont.processed.metab{gui.controls.Selected}=op_freqshift(MRSCont.processed.metab{gui.controls.Selected} ,1);
            if isfield(MRSCont.processed.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed.metab{gui.controls.Selected}.manual,'f')
                MRSCont.processed.metab{gui.controls.Selected}.manual.f = MRSCont.processed.metab{gui.controls.Selected}.manual.f + 1;
            else
                MRSCont.processed.metab{gui.controls.Selected}.manual.f =  1;
            end
        end
        if gui.controls.negFreqShift
            MRSCont.processed.metab{gui.controls.Selected}=op_freqshift(MRSCont.processed.metab{gui.controls.Selected} ,-1);
            if isfield(MRSCont.processed.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed.metab{gui.controls.Selected}.manual,'f')
                MRSCont.processed.metab{gui.controls.Selected}.manual.f = MRSCont.processed.metab{gui.controls.Selected}.manual.f - 1;
            else
                MRSCont.processed.metab{gui.controls.Selected}.manual.f =  -1;
            end
        end
        if isfield(MRSCont.processed.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed.metab{gui.controls.Selected}.manual,'ph0')
            if abs(MRSCont.processed.metab{gui.controls.Selected}.manual.ph0) > 360
                MRSCont.processed.metab{gui.controls.Selected}.manual.ph0 = MRSCont.processed.metab{gui.controls.Selected}.manual.ph0 + sign(MRSCont.processed.metab{gui.controls.Selected}.manual.ph0)*(-1)*360;
            end
        end
        
        if isfield(MRSCont, 'processed_no_align')
            if gui.controls.posPhase0Shift
                MRSCont.processed_no_align.metab{gui.controls.Selected}=op_addphase(MRSCont.processed_no_align.metab{gui.controls.Selected} ,5);
                if isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}.manual,'ph0')
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0 = MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0 + 5;
                else
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0 =  5;
                end
            end
            if gui.controls.negPhase0Shift
                MRSCont.processed_no_align.metab{gui.controls.Selected}=op_addphase(MRSCont.processed_no_align.metab{gui.controls.Selected} ,-5);
                if isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}.manual,'ph0')
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0 = MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0 - 5;
                else
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0 =  -5;
                end
            end
            if gui.controls.posPhase1Shift
                MRSCont.processed_no_align.metab{gui.controls.Selected}=op_addphase(MRSCont.processed_no_align.metab{gui.controls.Selected} ,0,0.00001);
                if isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}.manual,'ph1')
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph1 = MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph1 + 0.00001;
                else
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph1 =  0.00001;
                end
            end
            if gui.controls.negPhase1Shift
                MRSCont.processed_no_align.metab{gui.controls.Selected}=op_addphase(MRSCont.processed_no_align.metab{gui.controls.Selected} ,0,-0.00001);
                if isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}.manual,'ph1')
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph1 = MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph1 -0.00001;
                else
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph1 =  -0.00001;
                end
            end
            if gui.controls.flipPhase
                MRSCont.processed_no_align.metab{gui.controls.Selected}=op_addphase(MRSCont.processed_no_align.metab{gui.controls.Selected} ,180);
                if isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}.manual,'ph0')
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0 = MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0 + 180;
                else
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0 =  180;
                end
            end
            if gui.controls.posFreqShift
                MRSCont.processed_no_align.metab{gui.controls.Selected}=op_freqshift(MRSCont.processed_no_align.metab{gui.controls.Selected} ,1);
                if isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}.manual,'f')
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.f = MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.f + 1;
                else
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.f =  1;
                end
            end
            if gui.controls.negFreqShift
                MRSCont.processed_no_align.metab{gui.controls.Selected}=op_freqshift(MRSCont.processed_no_align.metab{gui.controls.Selected} ,-1);
                if isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}.manual,'f')
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.f = MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.f - 1;
                else
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.f =  -1;
                end
            end
            if isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}, 'manual') && isfield(MRSCont.processed_no_align.metab{gui.controls.Selected}.manual,'ph0')
                if abs(MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0) > 360
                    MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0 = MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0 + sign(MRSCont.processed_no_align.metab{gui.controls.Selected}.manual.ph0)*(-1)*360;
                end
            end

        end
    end
        
%%% 3. VISUALIZATION PART OF THIS TAB %%%
    temp_match = MRSCont.plot.processed.match;
    MRSCont.plot.processed.match = 2;
    if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
        temp = osp_plotProcess(MRSCont, gui.controls.Selected,Selection,SubSpec,Exp); %Create figure
    elseif isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
        temp = osp_plotProcess(MRSCont, gui.controls.Selected,Selection,SubSpec,Exp,gui.controls.act_x); %Create figure
    else
        temp = osp_plotProcess(MRSCont, gui.controls.Selected,Selection,SubSpec,Exp,[gui.controls.act_x gui.controls.act_y gui.controls.act_z]); %Create figure
    end
    MRSCont.plot.processed.match = temp_match;
    %Delete old content
    delete(gui.layout.proAlgn.Children.Children)
    %Fill window with new content        
    set( flipud(temp.Children(2).Children), 'Parent', gui.layout.proAlgn.Children ); % Update aligned and averaged plot
    set(  gui.layout.proAlgn.Children, 'XLim', temp.Children(2).XLim);
    set(  gui.layout.proAlgn.Children, 'YLimMode', 'auto');
    if isempty(gui.controls.YLimits)
        gui.controls.YLimits = gui.layout.proAlgn.Children.YLim;
    end

    if gui.controls.resetAxis
        set(  gui.layout.proAlgn.Children, 'YLimMode', 'auto');
        gui.controls.YLimits = gui.layout.proAlgn.Children.YLim;
    end

    range = abs(gui.controls.YLimits(1)) + abs(gui.controls.YLimits(2));
    if gui.controls.zoomOut
        gui.controls.YLimits = [gui.controls.YLimits(1)+sign(gui.controls.YLimits(1))*range*0.05 gui.controls.YLimits(2)+sign(gui.controls.YLimits(2))*range*0.05];
    end
    if gui.controls.zoomIn
        gui.controls.YLimits = [gui.controls.YLimits(1)-sign(gui.controls.YLimits(1))*range*0.05 gui.controls.YLimits(2)-sign(gui.controls.YLimits(2))*range*0.05];
    end
    if gui.controls.MoveUp
        gui.controls.YLimits = [gui.controls.YLimits(1)+range*0.01 gui.controls.YLimits(2)+range*0.01];
    end
    if gui.controls.MoveDown
        gui.controls.YLimits = [gui.controls.YLimits(1)-range*0.01 gui.controls.YLimits(2)-range*0.01];
    end

    % For debugging purposes
    % if isfield(MRSCont.processed.metab{gui.controls.Selected},'manual')
    %     MRSCont.processed.metab{gui.controls.Selected}.manual
    % end
    % gui.controls.YLimits
    
    set(gui.layout.proAlgn.Children, 'YLim', gui.controls.YLimits);
    set(  gui.layout.proAlgn.Children, 'XTick', temp.Children(2).XTick);
    close( temp );

    set(gui.controls.b_save_proTab{gui.process.Selected},'Callback',{@osp_onPrint,gui});
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class 
end