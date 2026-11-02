function osp_set_spm_path(spmPath)
%% osp_set_spm_path(spmPath)
%   Registers an SPM12 installation with Osprey permanently, so that it is
%   automatically added to the MATLAB path in all future Osprey sessions
%   (via OspreySettings.m) without requiring a manual `addpath` call every
%   time MATLAB is restarted.
%
%   Call this function once after installing SPM12, pointing it to the
%   folder that contains spm.m. The path is stored as a persistent MATLAB
%   preference ('Osprey', 'SPM12Path') and is also added to the path of
%   the current MATLAB session immediately.
%
%   USAGE:
%       osp_set_spm_path('C:\Users\me\Documents\MATLAB\spm12');
%       osp_set_spm_path('/Users/me/Documents/MATLAB/spm12');
%
%   INPUTS:
%       spmPath     = Full path to the SPM12 installation folder (i.e.
%                     the folder that contains spm.m).
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University)
%       goeltzs1@jhmi.edu

if nargin < 1 || isempty(spmPath)
    error(['Please provide the full path to your SPM12 installation, e.g.:' newline ...
           '    osp_set_spm_path(''C:\path\to\spm12'')']);
end

if exist(fullfile(spmPath,'spm.m'),'file') ~= 2
    error(['No spm.m found in ''%s''. Please point osp_set_spm_path to the ' ...
           'folder that contains spm.m (the root of your SPM12 installation).'], spmPath);
end

if length(spmPath) >= 4 && strcmpi(spmPath(end-3:end),'spm8')
    error('SPM8 was detected at ''%s'', but only SPM12 is supported by Osprey.', spmPath);
end

addpath(spmPath);
setpref('Osprey','SPM12Path',spmPath);
fprintf('SPM12 path set to ''%s''.\nOsprey will automatically add it to the path in all future MATLAB sessions.\n', spmPath);

end
