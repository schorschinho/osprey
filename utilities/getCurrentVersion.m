function VersionStruct = getCurrentVersion()
%% VersionStruct = getCurrentVersion()
%   This function reads the current version and release date from the
%   version.json file in the settings folder
%
%   USAGE:
%       VersionStruct = getCurrentVersion()
%
%   INPUTS: none
%       
%
%   OUTPUTS:
%       VersionStruct     = Struct with version number and release date.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2024-07-26)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2024-07-26: First version of the code.
    %% Get version
    fid = fopen(which(fullfile('settings','version.json')));  % Open and read
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
       
    if strcmp('win',osp_platform('filesys'))                % Correct single backslashes in Windows paths
        str = strrep(str,'\','\\');
    end
    
    % Following line replaces white space characters
    pattern = '[ \t\n]*"'; % Match zero or more spaces, tabs, or newlines, followed by a double quote
    replacement = '"'; % Replace the matched string with just a double quote
    str = regexprep(str, pattern, replacement);
    
    % Return struct
    VersionStruct  = jsondecode(str);
end