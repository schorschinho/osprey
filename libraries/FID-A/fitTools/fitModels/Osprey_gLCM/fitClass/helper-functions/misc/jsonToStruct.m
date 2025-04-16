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
                            if isfield(jsonStruct.Steps{ss}.parametrizations.(params{pp}), 'sd')
                                if ischar(jsonStruct.Steps{ss}.parametrizations.(params{pp}).sd) && strcmp(convertCharsToStrings(jsonStruct.Steps{ss}.parametrizations.(params{pp}).sd),'Inf')
                                    jsonStruct.Steps{ss}.parametrizations.(params{pp}).sd = Inf;
                                end
                            end
                            if isfield(jsonStruct.Steps{ss}.parametrizations.(params{pp}),'sc')
                                if ~iscell(jsonStruct.Steps{ss}.parametrizations.(params{pp}).sc.fix_factor)
                                    jsonStruct.Steps{ss}.parametrizations.(params{pp}).sc.fix_factor = {{jsonStruct.Steps{ss}.parametrizations.(params{pp}).sc.fix_factor'}};
                                end
                                if ~iscell(jsonStruct.Steps{ss}.parametrizations.(params{pp}).sc.ex)
                                    jsonStruct.Steps{ss}.parametrizations.(params{pp}).sc.ex = {{jsonStruct.Steps{ss}.parametrizations.(params{pp}).sc.ex'}};
                                end
                                if ~iscell(jsonStruct.Steps{ss}.parametrizations.(params{pp}).sc.sd)
                                    jsonStruct.Steps{ss}.parametrizations.(params{pp}).sc.sd = {{jsonStruct.Steps{ss}.parametrizations.(params{pp}).sc.sd'}};
                                end
                            end
                        end
                    end
                end
            end
        else
            for ss = 1 : length(jsonStruct.Steps)
                if isfield(jsonStruct.Steps(ss), 'parametrizations')
                        params = fieldnames(jsonStruct.Steps(ss).parametrizations);
                        if ~isempty(params)
                            for pp = 1 : length(params)
                                jsonStruct.Steps(ss).parametrizations.(params{pp}) = structfun(@transpose,jsonStruct.Steps(ss).parametrizations.(params{pp}),'UniformOutput',false);
                                if isfield(jsonStruct.Steps(ss).parametrizations.(params{pp}), 'type')
                                    jsonStruct.Steps(ss).parametrizations.(params{pp}).type = transpose(jsonStruct.Steps(ss).parametrizations.(params{pp}).type);
                                end
                                if isfield(jsonStruct.Steps(ss).parametrizations.(params{pp}), 'RegFun')
                                    jsonStruct.Steps(ss).parametrizations.(params{pp}).RegFun = transpose(jsonStruct.Steps(ss).parametrizations.(params{pp}).RegFun);
                                end
                                if isfield(jsonStruct.Steps(ss).parametrizations.(params{pp}), 'sd')
                                    if ischar(jsonStruct.Steps(ss).parametrizations.(params{pp}).sd) && strcmp(convertCharsToStrings(jsonStruct.Steps(ss).parametrizations.(params{pp}).sd),'Inf')
                                        jsonStruct.Steps(ss).parametrizations.(params{pp}).sd = Inf;
                                    end
                                end
                                if isfield(jsonStruct.Steps(ss).parametrizations.(params{pp}),'sc')
                                    if ~iscell(jsonStruct.Steps(ss).parametrizations.(params{pp}).sc.fix_factor)
                                        jsonStruct.Steps(ss).parametrizations.(params{pp}).sc.fix_factor = {{jsonStruct.Steps(ss).parametrizations.(params{pp}).sc.fix_factor'}};
                                    end
                                    if ~iscell(jsonStruct.Steps(ss).parametrizations.(params{pp}).sc.ex)
                                        jsonStruct.Steps(ss).parametrizations.(params{pp}).sc.ex = {{jsonStruct.Steps(ss).parametrizations.(params{pp}).sc.ex'}};
                                    end
                                    if ~iscell(jsonStruct.Steps(ss).parametrizations.(params{pp}).sc.sd)
                                        jsonStruct.Steps(ss).parametrizations.(params{pp}).sc.sd = {{jsonStruct.Steps(ss).parametrizations.(params{pp}).sc.sd'}};
                                    end
                                end
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
                if isfield(jsonStruct.parameters.(params{pp}), 'RegFun')
                    jsonStruct.parameters.(params{pp}).RegFun = transpose(jsonStruct.parameters.(params{pp}).RegFun);
                end
                if isfield(jsonStruct.parameters.(params{pp}),'parametrizations') && isfield(jsonStruct.parameters.(params{pp}).parametrizations, 'sc')
                    if ~iscell(jsonStruct.parameters.(params{pp}).parametrizations.sc.fix_factor)
                        jsonStruct.parameters.(params{pp}).parametrizations.sc.fix_factor = num2cell(jsonStruct.parameters.(params{pp}).parametrizations.sc.fix_factor);
                    end
                    if ~iscell(jsonStruct.parameters.(params{pp}).parametrizations.sc.ex)
                        % jsonStruct.parameters.(params{pp}).parametrizations.sc.ex = num2cell(jsonStruct.parameters.(params{pp}).parametrizations.sc.ex);
                        jsonStruct.parameters.(params{pp}).parametrizations.sc.ex = {{jsonStruct.parameters.(params{pp}).parametrizations.sc.ex'}};
                    end
                    if ~iscell(jsonStruct.parameters.(params{pp}).parametrizations.sc.sd)
                        % jsonStruct.parameters.(params{pp}).parametrizations.sc.sd = num2cell(jsonStruct.parameters.(params{pp}).parametrizations.sc.sd);
                        jsonStruct.parameters.(params{pp}).parametrizations.sc.sd = {{jsonStruct.parameters.(params{pp}).parametrizations.sc.sd'}};
                    end
                    for qq = 1 : size(jsonStruct.parameters.(params{pp}).parametrizations.sc.fix_factor,1)
                        if ~iscell(jsonStruct.parameters.(params{pp}).parametrizations.sc.fix_factor{qq})
                            jsonStruct.parameters.(params{pp}).parametrizations.sc.fix_factor{qq} = {jsonStruct.parameters.(params{pp}).parametrizations.sc.fix_factor{qq}};
                        end
                        if ~iscell(jsonStruct.parameters.(params{pp}).parametrizations.sc.ex{qq})
                            jsonStruct.parameters.(params{pp}).parametrizations.sc.ex{qq} = {jsonStruct.parameters.(params{pp}).parametrizations.sc.ex{qq}};
                        end                        
                        if ~iscell(jsonStruct.parameters.(params{pp}).parametrizations.sc.sd{qq})
                            jsonStruct.parameters.(params{pp}).parametrizations.sc.sd{qq} = {jsonStruct.parameters.(params{pp}).parametrizations.sc.sd{qq}};
                        end
                    end
                    if isfield(jsonStruct.parameters.(params{pp}).parametrizations, 'sd')
                        temp = zeros(length(jsonStruct.parameters.(params{pp}).parametrizations.sd),1);
                        for ss = 1 : length(jsonStruct.parameters.(params{pp}).parametrizations.sd)
                            if ischar(jsonStruct.parameters.(params{pp}).parametrizations.sd{ss}) && strcmp(convertCharsToStrings(jsonStruct.parameters.(params{pp}).parametrizations.sd{ss}),'Inf')
                                temp(ss) = Inf;
                            else
                                temp(ss) = cell2mat(jsonStruct.parameters.(params{pp}).parametrizations.sd(ss));
                            end
                        end
                        jsonStruct.parameters.(params{pp}).parametrizations.sd = temp;
                    end
                    if isfield(jsonStruct.parameters.(params{pp}).parametrizations, 'ub')
                        temp = zeros(length(jsonStruct.parameters.(params{pp}).parametrizations.ub),1);
                        for ss = 1 : length(jsonStruct.parameters.(params{pp}).parametrizations.ub)
                            if ischar(jsonStruct.parameters.(params{pp}).parametrizations.ub{ss}) && strcmp(convertCharsToStrings(jsonStruct.parameters.(params{pp}).parametrizations.ub{ss}),'Inf')
                                temp(ss) = Inf;
                            else
                                temp(ss) = cell2mat(jsonStruct.parameters.(params{pp}).parametrizations.ub(ss));
                            end
                        end
                        jsonStruct.parameters.(params{pp}).parametrizations.ub = temp;
                    end
                    if isfield(jsonStruct.parameters.(params{pp}).parametrizations, 'lb') && iscell(jsonStruct.parameters.(params{pp}).parametrizations.lb)
                        temp = zeros(length(jsonStruct.parameters.(params{pp}).parametrizations.lb),1);
                        for ss = 1 : length(jsonStruct.parameters.(params{pp}).parametrizations.lb)
                            if ischar(jsonStruct.parameters.(params{pp}).parametrizations.lb{ss}) && strcmp(convertCharsToStrings(jsonStruct.parameters.(params{pp}).parametrizations.lb{ss}),'Inf')
                                temp(ss) = Inf;
                            else
                                temp(ss) = cell2mat(jsonStruct.parameters.(params{pp}).parametrizations.lb(ss));
                            end
                        end
                        jsonStruct.parameters.(params{pp}).parametrizations.lb = temp;
                    end
                end
            end
        end
    end
end
