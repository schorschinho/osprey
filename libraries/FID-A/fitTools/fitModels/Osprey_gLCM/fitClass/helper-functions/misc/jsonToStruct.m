function jsonStruct = jsonToStruct(jsonfile)
%% jsonStruct = jsonToStruct(jsonfile)
%  Reads a JSON-encoded file and returns the content as struct.
%
%   USAGE:
%       jsonStruct = jsonToStruct(jsonfile)
%
%   INPUTS:
%       jsonfile       =  json file to convert
%
%   OUTPUTS:
%       jsonStruct     = struct generated from json
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu

%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017) 
%%  Find gyromagnetic ratio
    
    if (isstring(jsonfile) && strcmp(jsonfile{1}(1:5),'which'))  ||...
        (ischar(jsonfile) && strcmp(jsonfile(1:5),'which'))             % Full path is not given we have to eval
        strdef = append('jsonfile = ', jsonfile, ';');
        eval(strdef);
    end

    if ~isfile(jsonfile)                                     % Intercept if file doesn't exist
        error('JSON file %s does not exist.', jsonfile);
    end
    
    
    fid = fopen(jsonfile);                                   % Open and read
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
    jsonStruct  = jsondecode(str);

    % We want all parameters as row vectors so we have to change this in
    % the model json
    if isfield(jsonStruct, 'Steps')
        if iscell(jsonStruct.Steps)
            for ss = 1 : length(jsonStruct.Steps)
                if isfield(jsonStruct.Steps{ss}, 'parametrizations')
                    params = fieldnames(jsonStruct.Steps{ss}.parametrizations);
                    if ~isempty(params)
                        for pp = 1 : length(params)
                            jsonStruct.Steps{ss}.parametrizations.(params{pp}) = structfun(@transpose,jsonStruct.Steps{ss}.parametrizations.(params{pp}),'UniformOutput',false);
                            if isfield(jsonStruct.Steps{ss}.parametrizations.(params{pp}), 'type')
                                jsonStruct.Steps{ss}.parametrizations.(params{pp}).type = transpose(jsonStruct.Steps{ss}.parametrizations.(params{pp}).type);
                            end
                            if isfield(jsonStruct.Steps{ss}.parametrizations.(params{pp}), 'RegFun')
                                jsonStruct.Steps{ss}.parametrizations.(params{pp}).RegFun = transpose(jsonStruct.Steps{ss}.parametrizations.(params{pp}).RegFun);
                            end
                        end
                    end
                end
            end
        else
            if isfield(jsonStruct.Steps, 'parametrizations')
                    params = fieldnames(jsonStruct.Steps.parametrizations);
                    if ~isempty(params)
                        for pp = 1 : length(params)
                            jsonStruct.Steps.parametrizations.(params{pp}) = structfun(@transpose,jsonStruct.Steps.parametrizations.(params{pp}),'UniformOutput',false);
                            if isfield(jsonStruct.Steps.parametrizations.(params{pp}), 'type')
                                jsonStruct.Steps.parametrizations.(params{pp}).type = transpose(jsonStruct.Steps.parametrizations.(params{pp}).type);
                            end
                            if isfield(jsonStruct.Steps.parametrizations.(params{pp}), 'RegFun')
                                jsonStruct.Steps.parametrizations.(params{pp}).RegFun = transpose(jsonStruct.Steps.parametrizations.(params{pp}).RegFun);
                            end
                        end
                    end
                end
        end
    end

    % We have to do the same for the indirect parametrization file
    if isfield(jsonStruct, 'parameters')
        params = fieldnames(jsonStruct.parameters);
        if ~isempty(params)
            for pp = 1 : length(params)
                jsonStruct.parameters.(params{pp}) = structfun(@transpose,jsonStruct.parameters.(params{pp}),'UniformOutput',false);
                if isfield(jsonStruct.parameters.(params{pp}), 'type')
                    jsonStruct.parameters.(params{pp}).type = transpose(jsonStruct.parameters.(params{pp}).type);
                end
            end
        end
    end
end
