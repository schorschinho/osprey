function osp_onNewAnalysis( ~, ~,gui)
%% osp_onNewAnalysis
%   Callback function on new analysis button click. Closes gui and start a new MRS analysis.
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
%       Dr. Peter Van Schuerbeek, UZ Brussel (VUB)
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-09-16: First version of the code.
%%% 1. CLOSE %%%
    % User wants to quit out of the application
            
delete( gui.figure );
            
CreateOspreyJob_app;
end