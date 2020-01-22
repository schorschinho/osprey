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
    if strcmp(EventData.Key, 'uparrow')
        OldValue = get( gui.layout.ListBox,'value');
        gui.controls.KeyPress = 1;
        if OldValue == 1
            set(gui.layout.ListBox, 'value', MRSCont.nDatasets );
        else
            set(gui.layout.ListBox, 'value', OldValue-1 );
        end
    end
    if strcmp(EventData.Key, 'downarrow')
        OldValue = get( gui.layout.ListBox,'value');
        gui.controls.KeyPress = 1;
        if OldValue == MRSCont.nDatasets
            set(gui.layout.ListBox, 'value', 1 );
        else
            set(gui.layout.ListBox, 'value', OldValue+1 );
        end
    end
end