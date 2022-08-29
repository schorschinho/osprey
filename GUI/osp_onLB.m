function osp_onLB( source, ~ ,gui)
%% osp_onLB
%   Callback function on process click on left x voxel direction button.
%
%
%   USAGE:
%       osp_onLeftX( ~, ~ ,gui);
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
%       2021-01-21: First version of the code.
%%% 1. INITIALIZE %%%
    MRSCont = getappdata(gui.figure,'MRSCont');  % Get MRSCont from hidden container in gui class 
    selectedTab = get(gui.layout.tabs, 'Selection');

    % User wants to process the data
%%% 2. UPDATEWINDOW %%%    

    selectedLB = get(source, 'Value');
    MRSCont.opts.cosmetics.LB = selectedLB;
    set(source, 'Tooltip', ['Exponential Linebroadening (Hz) ' num2str(selectedLB)]);
    setappdata(gui.figure,'MRSCont',MRSCont);   % Write MRSCont into hidden container in gui class 
    osp_updateFitWindow(gui);

end % osp_onZo