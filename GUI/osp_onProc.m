function osp_onProc( ~, ~ ,gui)
%% osp_onProc
%   Callback function on process button click.
%
%
%   USAGE:
%       osp_onProc( ~, ~ ,gui);
%
%   INPUT:      gui      = gui class containing all handles and the MRSCont 
%
%   OUTPUT:     Changes in gui parameters and MRSCont are written into the
%               gui class
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
    set(gui.figure,'HandleVisibility','off');
    gui.layout.tabs.TabEnables{2} = 'on';
    gui.layout.tabs.Selection  = 2;
    [gui,MRSCont] = osp_processingWindow(gui,MRSCont);
    % User wants to process the data
%%% 2. CALL OPSREYPROCESS %%%    
    MRSCont = OspreyProcess(MRSCont);
    delete(gui.layout.dummy);
    gui.process.Number = length(fieldnames(MRSCont.processed));
    gui.process.Names = fieldnames(MRSCont.processed);
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
    set(gui.figure,'HandleVisibility','on');
%%% 3. INITIALIZE OUTPUT WINDOW %%%    
    osp_iniProcessWindow(gui);
    gui.layout.b_proc.Enable = 'off';
    gui.layout.b_fit.Enable = 'on';
end % onProc