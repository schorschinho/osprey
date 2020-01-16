function osp_onExit( ~, ~,gui)
%% osp_onExit
%   Callback function on exit button click. Saves MRSCont and closes gui.
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
%%% 1. SAVE AND CLOSE %%%
    % Save the output structure to the output folder
    % Determine output folder
    MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class    
    outputFolder    = MRSCont.outputFolder;
    outputFile      = MRSCont.outputFile;
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end
    save(fullfile(outputFolder, outputFile), 'MRSCont');
    % User wants to quit out of the application
    delete( gui.figure );
end % onExit