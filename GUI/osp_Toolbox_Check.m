function [hasSPM,OspreyVersion] = osp_Toolbox_Check (Module,ToolChecked)
%% [hasSPM] = osp_Toolbox_Check (Module,ToolChecked)
%   This function checks the availabilty of the required MATLAB toolboxes
%   and SPM versions. Adds the version number of Osprey.
%
%   USAGE:
%      [hasSPM] = osp_Toolbox_Check (Module,ToolChecked)
%
%   INPUTS:
%       Module     = String with the Module name
%                    Options:  OspreyGUI
%                              OspreyProcess
%                              OspreyFit
%                              OspreyCoreg
%                              OspreySeg
%      ToolChecked = Flag whether Toolboxes have been checked before. 
%
%   OUTPUTS:
%       hasSPM     = SPM flag.
%
%   AUTHOR:
%       Helge Zoellner (Johns Hopkins University, 2020-05-15)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-05-15: First version of the code.
%%
%%% 1. GET SPMPATH AND TOOLBOXES%%%
OspreyVersion = 'Osprey 1.0.2';
fprintf(['Timestamp %s ' OspreyVersion '  ' Module '\n'], datestr(now,'mmmm dd, yyyy HH:MM:SS'));
addons = matlab.addons.installedAddons;
available = cellstr(table2cell(addons(:,1)));
lic = strcmp({'Enabled'}, addons.Properties.VariableNames);
if ~isempty(lic)
    enabled = table2cell(addons(:,lic==1));
else
    enabled = table2cell(addons(:,1));
end

[settingsFolder,~,~] = fileparts(which('OspreySettings.m'));
allFolders      = strsplit(settingsFolder, filesep);
ospFolder       = strjoin(allFolders(1:end-1), filesep); % parent folder (= Osprey folder)
 
% SPM
if isfile(fullfile(ospFolder,'GUI','SPMpath.mat')) % Load path to SPM
    load(fullfile(ospFolder,'GUI','SPMpath.mat'),'SPMpath')
    if exist(SPMpath, 'dir')
        spmversion = SPMpath;
    else  %This isn't the right SPM path (maybe it was copied from another machine)
        spmversion = uipickfiles('FilterSpec',ospFolder,'REFilter', '\','NumFiles',1,'Prompt','Select your SPM-folder (Will be saved in SPMpath.mat file in the GUI folder)');
        spmversion = spmversion{1};
        SPMpath = spmversion;
        save(fullfile(ospFolder,'GUI','SPMpath.mat'),'SPMpath');
    end
else %Set path to SPM
    spmversion = uipickfiles('FilterSpec',ospFolder,'REFilter', '\','NumFiles',1,'Prompt','Select your SPM-folder (Will be saved in SPMpath.mat file in the GUI folder)');
    spmversion = spmversion{1};
    SPMpath = spmversion;
    save(fullfile(ospFolder,'GUI','SPMpath.mat'),'SPMpath');
end
% Check if SPM12 is installed
if isempty(spmversion)
    hasSPM = 0;
elseif strcmpi(spmversion(end-3:end),'spm8')
    available{end+1} = 'SPM8';
    enabled{end+1} = false;
    hasSPM = 0;
else
    available{end+1} = 'SPM12';
    enabled{end+1} = true;
    hasSPM = 1;
    rmpath(genpath([spmversion filesep 'external' filesep 'fieldtrip']));
end 

try
    if ~isempty(lic)
        available(find(cellfun(@(a)~isempty(a)&&a<1,enabled)), :) = [];
    end    

    %%% 2. CHECK AVAILABILTY %%%
    switch Module
        case 'OspreyGUI'
            ModuleString = 'fully run \bfOsprey';
            neededGlobal = {'Optimization Toolbox', 'Statistics and Machine Learning Toolbox','SPM12'};
            neededSpecific = {'Widgets Toolbox', 'GUI Layout Toolbox'};
        case 'OspreyProcess'
            ModuleString = 'run \bfOspreyProcess';
            neededGlobal = {'Optimization Toolbox', 'Statistics and Machine Learning Toolbox','SPM12'};
            neededSpecific = {'Optimization Toolbox', 'Statistics and Machine Learning Toolbox'}; 
        case 'OspreyFit'
            ModuleString = 'run \bfOspreyFit';
            neededGlobal = {'Optimization Toolbox', 'Statistics and Machine Learning Toolbox','SPM12'};
            neededSpecific = {'Optimization Toolbox', 'Statistics and Machine Learning Toolbox'};
        case 'OspreyCoreg'
            ModuleString = 'run \bfOspreyCoreg';
            neededGlobal = {'Optimization Toolbox', 'Statistics and Machine Learning Toolbox','SPM12'};
            neededSpecific = {'SPM12'};
        case 'OspreySeg'
            ModuleString = 'run \bfOspreySeg';
            neededGlobal = {'Optimization Toolbox', 'Statistics and Machine Learning Toolbox','SPM12'};
            neededSpecific = {'SPM12'};        
        otherwise
            ModuleString = ['run \bf' Module];
            neededGlobal = {'Optimization Toolbox', 'Statistics and Machine Learning Toolbox','SPM12'};
            neededSpecific = cellstr({});
    end
    
    %To account for the re-naming of new downloads of the Widget Toolbox
    %for Matlab versions earlier than 2020b, while maintaining
    %functionality for older downloads, we need to check for all naming
    %conventions of the Widgets Toolbox HZ
    for tb = 1 : length(available)
        if contains(available{tb},'Widgets Toolbox')
            available{tb} = 'Widgets Toolbox';
        end
    end
    missingSpecific = setdiff(neededSpecific,available);
    missing = setdiff(neededGlobal,available); 

    %%% 3. CREATE WARNING MESSAGES %%%
    if ~ToolChecked
        warning = cellstr({});
        warning_count = 1;
        if ~isempty(missing) || ~isempty(missingSpecific)
            opts.Interpreter = 'tex';
            opts.WindowStyle = 'modal';
            warning{1} = ['The following toolboxes are missing to ' ModuleString '\rm:'];
            for i = 1 : length(missing)
                warning{i+1} = ['\bf' missing{i} '\rm'];
            end
            warning_count = warning_count +length(missing) + 1;
            warning{warning_count} = ['Please install them to ' ModuleString '\rm'];
            warning_count = warning_count + 1;
            if ~isempty(missingSpecific)
                warning{warning_count} = ['The following toolboxes are missing to run ' Module ':']; 
                warningc = ['Please install and include the following toolboxes to use ' Module ':'];
                for i = 1 : length(missingSpecific)
                    warning{warning_count + i} = ['\bf' missingSpecific{i} '\rm'];
                    warningc = [warningc ' ' missingSpecific{i}];
                end
                warning{warning_count + length(missingSpecific) + 1} = ['Please install them to use \bf' Module '\rm'];
                warndlg(warning,'Missing Toolboxes',opts);
                error(warningc);
            end    
            warndlg(warning,'Missing Toolboxes',opts);
        end
    end

catch %If the MATLAB version pre-dates the inmplementation of matlab.addons.installedAddons
    warning = cellstr({});
    warning{1} = 'Your current MATLAB version does not allow the automated toolbox check. We assume that all required toolboxes are available.';
    warndlg(warning,'Automated toolbox check not working.',opts);
end   

end