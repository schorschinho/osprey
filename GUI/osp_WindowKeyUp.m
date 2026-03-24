function osp_WindowKeyUp(~,~,gui)
%% osp_WindowKeyUp
%   This function is triggered when a key is released. It refreshes
%   the the listbox.
%
%
%   USAGE:
%       osp_WindowKeyUp(~,~,gui);
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
%%% 1. GET HANDLES %%%
    gui.controls.KeyPress = 0;
    osp_updateListBox(gui);
end