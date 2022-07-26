function LCMparam = osp_readlcm_control(controlFile)
%% LCMparam = osp_readlcm_control(controlFile)
%   This function reads a LCModel compatible .control file and saves the
%   parameters as a struct that can be modified.
%
%   USAGE:
%       LCMparam = osp_readlcm_control(controlFile)
%
%   INPUTS:
%      controlFile  = LCModel .control file
%
%   OUTPUTS:
%      LCMparam     = struct with LCModel control parameters
%
%   AUTHORS:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2021-08-03)
%       goeltzs1@jhmi.edu
%   

% Throw error if no input is provided
if nargin < 1
    error('No control file provided to read.');
end

% Open the provided file, create empty cell to store the parsed entries
fid = fopen(controlFile);
C = {};

% The first line is always '$LCMODL', the last line always '$END'
firstLine = '$LCMODL';
lastLine  = '$END';

% Get lines and continue until we get to the last line indicator
tline = fgets(fid);
while isempty(strfind(tline, lastLine))
    tline = fgets(fid);
    
    % Parse every line between the start and the end line
    if (isempty(strfind(tline ,firstLine)) + isempty(strfind(tline, lastLine)) == 2)
        C{end+1}= strtrim(strsplit(tline, {'=',','}));
    end
    
end

% Loop over everything that has been parsed
LCMparam = struct;
for rr = 1:length(C)
    % Most entries should be pairs, but there are cases where control
    % parameters have been entered on the same line and separated by
    % commas, or parameters like 'title' where the value may contain the
    % delimiters '=' and ','
    if length(C{rr}) == 2
        LCMparam = parseControlFileLine(LCMparam, C{rr});
    else
        if strcmpi(C{rr}{1}, 'title')% If the title field is split, join back together
            title = strjoin(C{rr}(2:end));
            LCMparam.title = title;
        elseif strlength(strrep(C{rr}{1},"'",''))<=1
            % Ignore line
        elseif strcmp(C{rr}{1}(1:6),'CHSIMU') % simulation parameters contains '='
            P{1} = C{rr}{1};
            for pp=2:length(C{rr})
                if length(C{rr}{pp})>4 && (strcmp(C{rr}{pp}(end-3:end),'FWHM') || strcmp(C{rr}{pp}(end-2:end),'AMP'))
                    C{rr}{pp} = [C{rr}{pp},'=']; %Re-add the "="
                end
            end
            P(2) = join(C{rr}(2:end));
            LCMparam = parseControlFileLine(LCMparam,P);
        elseif length(C{rr}) == 3% If length=3, then assume 2-value input
            P{1} = C{rr}{1};
            P{2} = [C{rr}{2},',',C{rr}{3}];
            LCMparam = parseControlFileLine(LCMparam,P);
        else
            % If the number of entries in this line is even, proceed, otherwise
            % throw an error
            if mod(length(C{rr}),2) == 0
                % Evaluate pairwise
                for pp = 1:length(C{rr})/2
                    P{1} = C{rr}{2*pp-1};
                    P{2} = C{rr}{2*pp};
                    LCMparam = parseControlFileLine(LCMparam, P);
                end
            else
                Out = join(C{rr});
                error('Invalid control file line: %s', Out{1});
            end
        
        end
        
    end   

end

fclose(fid);

end

function LCMparam = parseControlFileLine(LCMparam, keyValuePair)
    % Test whether the input key has parentheses in it, for example chomit
    key = keyValuePair{1};
    expr = '(\w+)';
    [tokens,matches] = regexp(key,expr,'tokens','match');
    % If this regular expression search returns two tokens, it means that
    % the control parameter name is a vector. We'll save the entries as a
    % cell.
    if length(tokens) == 1
        LCMparam.(key) = keyValuePair{2};
    elseif length(tokens) == 2
        LCMparam.(tokens{1}{1}){str2num(tokens{2}{1})} = keyValuePair{2};
    end
    
end
