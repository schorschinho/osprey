function GenerateBIDSStructure(nSubs, nSes, dataTypes)
%% GenerateBIDSStructure(nSubs, nSes, dataTypes)
%   This function generates a coarse empty BIDS structure that can be used to
%   pre-sort imaging and MRS data into a structure than can be reasonably
%   well BIDS-ified with tools like BIDScoin.
%
%   USAGE:
%       GenerateBIDSStructure(nSubs, nSes, dataTypes)
%
%   INPUTS:
%       nSubs     = (integer) number of subjects
%       nSes      = (integer) number of sessions
%       dataTypes = (cell) names of BIDS datatypes
%
%   OUTPUTS:
%       BIDS folder structure with nSubs subjects, nSes sessions, and data
%       type folders.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2023-10-06)
%       goeltzs1@jhmi.edu


if nargin<3
    dataTypes = {'anat'};
    if nargin<2
        nSes = 1;
        if nargin<1
            nSubs = 1;
        end
    end
end

% Check validity of input data types
dataTypes = checkValidDataTypes(dataTypes);

% Create subject layer
lengthCell = length(dataTypes) * nSubs * nSes;
foldersToCreate = cell(lengthCell, 1);
folderCreateIndex = 1;
for subIdx = 1:nSubs

    % Create subject string
    subStr = sprintf('sub-%03d', subIdx);
    
    % Create session layer
    for sesIdx = 1:nSes

        % Create session string
        sesStr = sprintf('ses-%02d', sesIdx);

            % Create data type layer
            for dataTypeIdx = 1:length(dataTypes)

                % Create data type string
                dataTypeStr = dataTypes{dataTypeIdx};
                
                % Assemble the string
                if nSes == 1
                    % Omit session layer if only one session
                    stringToCreate = fullfile(subStr, dataTypeStr);
                else
                    stringToCreate = fullfile(subStr, sesStr, dataTypeStr);
                end
                
                % Save and advance counter
                foldersToCreate{folderCreateIndex} = stringToCreate;
                folderCreateIndex = folderCreateIndex + 1;
            end
    end
end

% Create folders
for ff = 1:numel(foldersToCreate)
    if 2~=exist(foldersToCreate{ff},'dir')
        mkdir(foldersToCreate{ff});
    end
end

    function dataTypes = checkValidDataTypes(dataTypes)

        % Catch invalid input
        if ~iscell(dataTypes)
            error('dataTypes input must be a cell array of valid BIDS datatypes.')
        end

        % List of valid BIDS data types
        validDataTypes = {'anat', 'func', 'fmap', 'dwi', 'perf', 'eeg', ...
                          'meg', 'ieeg', 'beh', 'pet', 'micr', 'nirs', ...
                          'motion', 'mrs'};

        % Determine match
        match = ismember(dataTypes, validDataTypes);
        if ~all(match) % Checks if all elements of 'match' are TRUE
            % Find indices that are 0
            invalidElements = find(~match);
            for rr = 1:length(invalidElements)
                fprintf('Invalid BIDS datatype: %s \n', dataTypes{invalidElements(rr)});
            end
            error('Invalid BIDS datatypes detected. Abort!')
        end

    end

end