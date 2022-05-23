function [hasSPM,OspreyVersion] = osp_CheckRunPreviousModule(MRSCont, module)
%% osp_CheckRunPreviousModule(MRSCont, module)
%   This function checks that all required modules have been run prior to
%   executing the module name that is passed as an argument to this
%   function.
%
%   There is no output, the function simply throws an error highlighting
%   which module needs to be executed before.
%
%   USAGE:
%      osp_CheckRunPreviousModule(module)
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%       module      = String with the module name
%                     Options: OspreyProcess
%                              OspreyFit
%                              OspreyCoreg
%                              OspreySeg
%
%   AUTHOR:
%       Georg Oeltzschner (Johns Hopkins University, 2021-07-07)
%       goeltzs1@jhmi.edu
%

switch module
    case 'OspreyLoad'
        requiredModules     = {'OspreyJob'};
        requiredModuleVerbs = {'defined (job)'};
        currentModuleVerb   = 'load';
    case 'OspreyProcess'
        requiredModules     = {'OspreyLoad'};
        requiredModuleVerbs = {'defined (job)','loaded'};
        currentModuleVerb   = 'process';
        
    case 'OspreyFit'
        requiredModules     = {'OspreyLoad', 'OspreyProcess'};
        requiredModuleVerbs = {'defined (job)','loaded', 'processed'};
        currentModuleVerb   = 'fit';

    case 'OspreyCoreg'
        requiredModules     = {'OspreyLoad'};
        requiredModuleVerbs = {'defined (job)','loaded'};
        currentModuleVerb   = 'coregister';
        
    case 'OspreySeg'
        requiredModules     = {'OspreyLoad', 'OspreyCoreg'};
        requiredModuleVerbs = {'defined (job)','loaded', 'coregistered'};
        currentModuleVerb   = 'segment';
        
    case 'OspreyQuantify'
        requiredModules     = {'OspreyLoad', 'OspreyProcess', 'OspreyFit'};
        requiredModuleVerbs = {'defined (job)','loaded', 'processed', 'fit'};
        currentModuleVerb   = 'quantify';
        
    case 'OspreyOverview'
        requiredModules     = {'OspreyLoad', 'OspreyProcess'};
        requiredModuleVerbs = {'defined (job)','loaded', 'OspreyProcess'};
        currentModuleVerb   = 'create overview of';
        
    otherwise
        error('Unknown module name: %s', module);
end

% Loop over required module names
for rr = 1:length(requiredModules)
    % Flag names are the name of the module minus the 'Osprey' string.
    flagName = strrep(requiredModules{rr}, 'Osprey', 'did');
       
    if ~MRSCont.flags.(flagName)
        error('Trying to %s data, but data has not been %s yet. Run %s first.', currentModuleVerb, requiredModuleVerbs{rr}, requiredModules{rr});
    end
        
end

%Do the toolbox check here
OspreyVersion = 'Osprey 2.1.0';
hasSPM = 1;
[hasSPM,OspreyVersion ] = osp_Toolbox_Check (module,MRSCont.flags.isGUI);

end