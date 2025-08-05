function osp_onModelPickq( source, ~ ,gui)
%% osp_onModelPickq
%   Callback function on process click on model pick slider.
%
%
%   USAGE:
%       osp_onModelPickq( source, ~ ,gui);
%
%   INPUT:      gui      = gui class containing all handles and the MRSCont 
%
%   OUTPUT:     Changes in gui parameters and MRSCont are written into the
%               gui class
%
%
%   AUTHORS:
%       Dr. Chris Davies-Jenkins(Johns Hopkins University, 2025-07-03)
%       cdavies9@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%%% 1. INITIALIZE %%%
    MRSCont = getappdata(gui.figure,'MRSCont');  % Get MRSCont from hidden container in gui class 

    % User wants to process the data
%%% 2. UPDATEWINDOW %%%    
    selectedModel = round(get(source, 'Value'));
    set(source, 'Value', selectedModel);
    set(source, 'Tooltip', ['Model Chosen ' num2str(selectedModel)]);
    osp_updateQuantifyWindow(gui);


end % osp_onModelStep