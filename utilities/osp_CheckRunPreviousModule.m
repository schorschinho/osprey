function osp_CheckRunPreviousModule(MRSCont, module)
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
    case 'OspreyProcess'
        requiredModules     = {'OspreyLoad'};
        requiredModuleVerbs = {'loaded'};
        currentModuleVerb   = 'process';
        
    case 'OspreyFit'
        requiredModules     = {'OspreyLoad', 'OspreyProcess'};
        requiredModuleVerbs = {'loaded', 'processed'};
        currentModuleVerb   = 'fit';

    case 'OspreyCoreg'
        requiredModules     = {'OspreyLoad'};
        requiredModuleVerbs = {'loaded'};
        currentModuleVerb   = 'coregister';
        
    case 'OspreySeg'
        requiredModules     = {'OspreyLoad', 'OspreyCoreg'};
        requiredModuleVerbs = {'loaded', 'coregistered'};
        currentModuleVerb   = 'segment';
        
    case 'OspreyQuantify'
        requiredModules     = {'OspreyLoad', 'OspreyProcess', 'OspreyFit'};
        requiredModuleVerbs = {'loaded', 'processed', 'fit'};
        currentModuleVerb   = 'quantify';
        
    case 'OspreyOverview'
        requiredModules     = {'OspreyLoad', 'OspreyProcess'};
        requiredModuleVerbs = {'loaded', 'OspreyProcess'};
        currentModuleVerb   = 'create overview of';
        
    otherwise
        error('Unknown module name: %s', module);
end

% Loop over required module names
for rr = 1:length(requiredModules)
    % Flag names are the name of the module minus the 'Osprey' string.
    flagName = strrep(requiredModules{rr}, 'Osprey', 'did');
    
    % For backwards compatibility:
    if ~isfield(MRSCont.flags, 'didLoad') 
        flagName = 'didLoadData';
    end
    
    if ~MRSCont.flags.(flagName)
        error('Trying to %s data, but data has not been %s yet. Run %s first.', currentModuleVerb, requiredModuleVerbs{rr}, requiredModules{rr});
    end
        
end


end