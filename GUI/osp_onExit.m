function osp_onExit( ~, ~,gui)
%% osp_onExit
%   Callback function on exit button click. Closes gui.
%
%
%   USAGE:
%       osp_onExit( ~, ~ ,gui);
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
%%% 1. CLOSE %%%
    % User wants to quit out of the application
    
    
    answer = questdlg('Do you want to quit Osprey or to start a new analysis?','Exit Osprey','Quit','New analysis','Cancel','Cancel');
    
    switch answer
        case 'Quit'
            [ospFFolder,~,~] = fileparts(which('Osprey.m'));
            curdir = cd(ospFFolder);
            load('startpath.mat')
            path(startpath);
            delete('startpath.mat')
            cd(curdir)
    
            delete( gui.figure );
        case 'New analysis'
            [ospFFolder,~,~] = fileparts(which('Osprey.m'));
            curdir = cd(ospFFolder);
            load('startpath.mat')
            path(startpath);
            delete('startpath.mat')
            cd(curdir)
            
            delete( gui.figure );
            
            Osprey;
        case 'Cancel'
            return
    end
end % onExit