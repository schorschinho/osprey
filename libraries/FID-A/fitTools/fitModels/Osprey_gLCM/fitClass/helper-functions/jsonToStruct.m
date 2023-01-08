function jsonStruct = jsonToStruct(jsonfile)
% Reads a JSON-encoded file and returns the content as struct.

% Intercept if file doesn't exist
if ~isfile(jsonfile)
    error('JSON file %s does not exist.', jsonfile);
end

% Open and read
fid = fopen(jsonfile); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 

% Correct single backslashes in Windows paths
if strcmp('win',osp_platform('filesys'))
    str = strrep(str,'\','\\');
end

% Following line replaces white space characters, but only introduced in
% MATLAB R2020b - commented out for now and seemed to work, might need a
% version-agnostic workaround?
% str = replace(str, whitespacePattern + '"', '"');

% Return struct
jsonStruct  = jsondecode(str);

end
