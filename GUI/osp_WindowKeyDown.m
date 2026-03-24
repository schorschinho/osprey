function osp_WindowKeyDown(~,EventData,gui) 
%% osp_WindowKeyDown
%   This function is triggered when key up and down is pressed. It refreshes
%   the listbox.
%
%
%   USAGE:
%       osp_WindowKeyDown(,EventData,gui);
%
%   INPUT:  EventData      = Key press event
%           gui            = gui class containing all handles and the MRSCont             
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
    % 28 leftarrow
    % 29 rightarrow
    % 30 uparrow
    % 31 downarrow
% Get MRSCont from hidden class container
    MRSCont = getappdata(gui.figure,'MRSCont');    
    if strcmp(EventData.Key, 'uparrow') % scrolling up
        OldValue = get( gui.layout.ListBox,'value');
        gui.controls.KeyPress = 1;
        if OldValue == 1
            set(gui.layout.ListBox, 'value', MRSCont.nDatasets );
        else
            set(gui.layout.ListBox, 'value', OldValue-1 );
        end
    end
    if strcmp(EventData.Key, 'downarrow') % scrolling down
        OldValue = get( gui.layout.ListBox,'value');
        gui.controls.KeyPress = 1;
        if OldValue == MRSCont.nDatasets
            set(gui.layout.ListBox, 'value', 1 );
        else
            set(gui.layout.ListBox, 'value', OldValue+1 );
        end
    end
    if strcmp(EventData.Key, 'leftarrow') % excluding spectrum from OspreyQuantify and OspreyOverview
        Idx = get( gui.layout.ListBox,'value');
        gui.controls.KeyPress = 1;
        if ~strcmp(gui.layout.RedFileList{Idx}(1),'<')
            gui.layout.RedFileList{Idx} = ['<HTML><s>' gui.layout.RedFileList{Idx} '</HTML></s>'];
            set(gui.layout.ListBox,'String',gui.layout.RedFileList)
        end
        if ~isfield(MRSCont,'exclude')
            MRSCont.exclude = Idx;
            if MRSCont.flags.didQuantify
                gui.controls.waitbar = waitbar(0,'Start','Name','Exclude spectra');
                waitbar(0,gui.controls.waitbar,'Exclude spectra')
                MRSCont.flags.didQuantify =0;
                MRSCont.flags.didOverview =0;
                MRSCont = rmfield(MRSCont,'quantify');
                MRSCont = rmfield(MRSCont,'overview');
                waitbar(0.33,gui.controls.waitbar,'Call OspreyQuantify')
                MRSCont = OspreyQuantify(MRSCont);
                waitbar(0.66,gui.controls.waitbar,'Call OspreyOverview')
                MRSCont = OspreyOverview(MRSCont);
                waitbar(1,gui.controls.waitbar,'Finished');
                pause(1);
                close(gui.controls.waitbar);
            end
        else if ~ismember(Idx, MRSCont.exclude)
            MRSCont.exclude = [MRSCont.exclude Idx];
            if MRSCont.flags.didQuantify
                gui.controls.waitbar = waitbar(0,'Start','Name','Exclude spectra');
                waitbar(0,gui.controls.waitbar,'Exclude spectra')
                MRSCont.flags.didQuantify =0;
                MRSCont.flags.didOverview =0;
                MRSCont = rmfield(MRSCont,'quantify');
                MRSCont = rmfield(MRSCont,'overview');
                waitbar(0.33,gui.controls.waitbar,'Call OspreyQuantify')
                MRSCont = OspreyQuantify(MRSCont);
                waitbar(0.66,gui.controls.waitbar,'Call OspreyOverview')
                MRSCont = OspreyOverview(MRSCont);
                waitbar(1,gui.controls.waitbar,'Finished');
                pause(1);
                close(gui.controls.waitbar);                
            end
            end

        end
    end
    if strcmp(EventData.Key, 'rightarrow') % Adding an excluded spectrum from OspreyQuantify and OspreyOverview
        Idx = get( gui.layout.ListBox,'value');
        gui.controls.KeyPress = 1;
        if strcmp(gui.layout.RedFileList{Idx}(1),'<')
            tempTxt = gui.layout.RedFileList{Idx};
            tempTxt = strrep(tempTxt,'<HTML><s>','');
            tempTxt = strrep(tempTxt,'</HTML></s>','');
            gui.layout.RedFileList{Idx} = tempTxt;
            set(gui.layout.ListBox,'String',gui.layout.RedFileList);
            MRSCont.exclude(MRSCont.exclude == Idx) = [];
            if MRSCont.flags.didQuantify
                gui.controls.waitbar = waitbar(0,'Start','Name','Include spectra');
                waitbar(0,gui.controls.waitbar,'Include spectra')
                MRSCont.flags.didQuantify =0;
                MRSCont.flags.didOverview =0;
                MRSCont = rmfield(MRSCont,'quantify');
                MRSCont = rmfield(MRSCont,'overview');
                waitbar(0.33,gui.controls.waitbar,'Call OspreyQuantify')
                MRSCont = OspreyQuantify(MRSCont);
                waitbar(0.66,gui.controls.waitbar,'Call OspreyOverview')
                MRSCont = OspreyOverview(MRSCont);
                waitbar(1,gui.controls.waitbar,'Finished');
                pause(1);
                close(gui.controls.waitbar);                
            end
        end
    end
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
end