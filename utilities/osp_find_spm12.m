function [spmPath] = osp_find_spm12(ospFolder)
%% [spmPath] = osp_find_spm12(ospFolder)
%   This function attempts to locate a valid SPM12 installation on the
%   current machine without requiring the user to manually call `addpath`
%   in every MATLAB session. It checks, in order:
%       1. Whether SPM is already on the MATLAB path.
%       2. A path previously registered with `osp_set_spm_path`, stored as
%          a persistent MATLAB preference ('Osprey', 'SPM12Path').
%       3. A short list of common installation locations, e.g. a `spm12`
%          folder next to the Osprey installation, or inside `userpath`.
%
%   This function does not modify the MATLAB path itself; it only reports
%   a candidate folder. Callers (e.g. OspreySettings.m) are responsible
%   for calling `addpath` if a path is returned.
%
%   USAGE:
%       [spmPath] = osp_find_spm12(ospFolder);
%
%   INPUTS:
%       ospFolder   = Full path to the Osprey root folder (optional). Used
%                     to probe for a sibling `spm12`/`spm` folder.
%
%   OUTPUTS:
%       spmPath     = Full path to a valid SPM12 installation (folder
%                     containing spm.m), or '' if none could be found.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University)
%       goeltzs1@jhmi.edu

if nargin < 1
    ospFolder = '';
end

spmPath = '';

% 1. Already on the MATLAB path?
onPathFolder = fileparts(which('spm'));
if isValidSPM12Folder(onPathFolder)
    spmPath = onPathFolder;
    return;
end

% 2. Previously registered via osp_set_spm_path?
if ispref('Osprey','SPM12Path')
    candidate = getpref('Osprey','SPM12Path');
    if isValidSPM12Folder(candidate)
        spmPath = candidate;
        return;
    end
end

% 3. Common installation locations.
candidates = {};
if ~isempty(ospFolder)
    parentFolder = fileparts(ospFolder);
    candidates = [candidates, {fullfile(parentFolder,'spm12'), fullfile(parentFolder,'spm')}];
end
upFolder = userpath;
if ~isempty(upFolder)
    candidates = [candidates, {fullfile(upFolder,'spm12'), fullfile(upFolder,'spm')}];
end

for rr = 1:length(candidates)
    if isValidSPM12Folder(candidates{rr})
        spmPath = candidates{rr};
        return;
    end
end

end


function tf = isValidSPM12Folder(folder)
% A folder is considered a valid SPM12 installation if it contains
% spm.m and is not an SPM8 installation.
tf = ~isempty(folder) && exist(fullfile(folder,'spm.m'),'file') == 2 && ...
     ~(length(folder) >= 4 && strcmpi(folder(end-3:end),'spm8'));
end
