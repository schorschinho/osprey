function osp_onModelStep( source, ~ ,gui)
%% osp_onModelStep
%   Callback function on process click on model step slider.
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

    % User wants to process the data
%%% 2. UPDATEWINDOW %%%    
    selectedModel = round(get(source, 'Value'));
    set(source, 'Value', selectedModel);
    set(source, 'Tooltip', ['Model Step ' num2str(selectedModel)]);
    osp_updateFitWindow(gui);


end % osp_onModelStep