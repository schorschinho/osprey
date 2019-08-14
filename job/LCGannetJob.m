function [MRSCont] = LCGannetJob(jobFile)
%% [MRSCont] = LCGannetJob(jobFile)
%   This function loads the LCGannet job defined in jobFile.
%   The job can be submitted in the following formats:
%       - .m
%
%   A valid LCGannet job contains four distinct classes of items:
%       1. basic information on the MRS sequence used
%       2. several settings for data handling and modeling
%       3. a list of MRS (and, optionally, structural imaging) data files 
%          to be loaded
%       4. an output folder to store the results and exported files
%
%   See the example files in the /exampledata/ folder for details.
%
%   USAGE:
%       [MRSCont] = LCGannetJob(jobFile);
%
%   INPUTS:
%       jobFile     = File containing a correct LCGannet job definition.
%                     Accepted file formats are .m.
%
%   OUTPUTS:
%       MRSCont     = LCGannet MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-07-15)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-07-15: First version of the code.

% Close any remaining open figures
close all;
warning('off','all');

%%% 1. INITIALISE DATA CONTAINER WITH DEFAULT SETTINGS
[MRSCont] = LCG_Settings;


%%% 2. CHECK JOB INPUT FILE FORMAT %%%
last3char   = jobFile((end-2):end);
last2char   = jobFile((end-1):end);
if strcmpi(last2char,'.m')
    jobFileFormat = 'm';
elseif strcmpi(last3char,'csv')
    jobFileFormat = 'csv';
else
    error('Unrecognized job file datatype. Job files need to end in .CSV or .M');
end


%%% 3a. LOAD AND RUN M FILE %%%
if strcmp(jobFileFormat, 'm')
    run(jobFile);
end

%%% 3b. LOAD CSV FILE %%%
if strcmp(jobFileFormat,'csv')
    jobCSV = readtable(jobFile, 'Delimiter', ',');
    % UTF-8 encoding may screw up the first column name; correct here
    if ~strcmp(jobCSV.Properties.VariableNames{1}, 'files')
        jobCSV.Properties.VariableNames{1} = 'files';
    end
    % Convert to a struct array
    jobStruct = table2struct(jobCSV);
    
    % Check whether the relevant fieldnames have been entered,
    % and save them as separate cells to be saved into the MRSCont
    % container.
    if isfield(jobStruct, 'files')
        files = {jobStruct.files};
    else
        error('Invalid job file! A job file needs to contain at least metabolite data in the field ''files''.');
    end
    if isfield(jobStruct, 'files_ref')
        files_ref = {jobStruct.files_ref};
    end
    if isfield(jobStruct, 'files_w')
        files_w = {jobStruct.files_w};
    end
    if isfield(jobStruct, 'files_nii')
        files_nii = {jobStruct.files_nii};
    end
    if isfield(jobStruct, 'outputFolder')
        outputFolder = {jobStruct.outputFolder};
    else
        error('Invalid job file! A job file needs to specify an output folder.');
    end
end


%%% 4. SAVE SETTINGS INTO MRSCONT %%%
switch seqType
    case 'unedited'
        MRSCont.flags.isUnEdited    = 1;
    case 'MEGA'
        MRSCont.flags.isMEGA        = 1;
    case 'HERMES'
        MRSCont.flags.isHERMES      = 1;
    case 'HERCULES'
        MRSCont.flags.isHERCULES    = 1;
    otherwise
        error('Invalid job file! seqType must be ''unedited'', ''MEGA'', ''HERMES'', or ''HERCULES''.');
end
MRSCont.opts = opts;


%%% 5. SAVE FILE/FOLDER NAMES INTO MRSCONT %%%
% Make sure that the mandatory fields (metabolite data; output folder) are
% included in the job file.
if exist('files','var')
    MRSCont.files = files;
else
    error('Invalid job file! A job file needs to contain at least metabolite data in the field ''files''.');
end
if exist('files_ref','var')
    MRSCont.files_ref = files_ref;
end
if exist('files_w','var')
    MRSCont.files_w = files_w;
end
if exist('files_nii','var')
    MRSCont.files_nii = files_nii;
end
if exist('outputFolder','var')
    MRSCont.outputFolder = outputFolder;
else
    error('Invalid job file! A job file needs to specify an output folder.');
end

% Check that each array has an identical number of entries
fieldNames = {'files', 'files_ref', 'files_w', 'files_nii'};
ctr = 0;
for kk = 1:length(fieldNames)
    if isfield(MRSCont, fieldNames{kk})
        if ~isempty(MRSCont.(fieldNames{kk}))
            ctr = ctr + 1;
            whichFieldNames{ctr} = fieldNames{kk};
            numDataSets(ctr)     = length(MRSCont.(fieldNames{kk}));
        end
    end
end

% Check whether the number of entries is identical
isUnique = unique(numDataSets);
if length(isUnique) ~= 1
    msg = sprintf('''%s'' has %i entries, but ', whichFieldNames{1}, numDataSets(1));
    for ll = 2:length(whichFieldNames)
        msg2 = sprintf(' ''%s'' has %i entries, ', whichFieldNames{ll}, numDataSets(ll));
        msg = strcat(msg,msg2);
    end
    
    error('MyComponent:invalidJob', ['Invalid job file! ' msg '\b.']);
end


%%% 6. SET FLAGS %%%
MRSCont.flags.didLoadJob    = 1;
MRSCont.loadedJob           = jobFile;


%%% 7. CHECK IF OUTPUT STRUCTURE ALREADY EXISTS IN OUTPUT FOLDER %%%
[~,jobfilename,jobfileext]  = fileparts(jobFile);
outputFile                  = strrep([jobfilename jobfileext], jobfileext, '.mat');
MRSCont.outputFile          = outputFile;
if isfile(fullfile(outputFolder, outputFile))
    fprintf('Your selected output folder ''%s'' already contains an LCGannet output structure: \t%s.\n', outputFolder, outputFile);
    fprintf('You are about to load the job: \t%s.\n', jobFile);
    askOverWriteJob = input('Do you want to overwrite the existing job (y/n)? (Warning: This will delete all data in the data container!)   [y]   ','s');
    if isempty(askOverWriteJob)
        askOverWriteJob = 'y';
    end
    if askOverWriteJob=='n' || askOverWriteJob=='N'
        disp('Aborted! No new job loaded.');
        return;
    elseif askOverWriteJob=='y' || askOverWriteJob=='y'
        disp('Continue with loading new job, overwriting existing job.');
    end
end

%%% 8. SAVE THE OUTPUT STRUCTURE TO THE OUTPUT FOLDER %%%
% Determine output folder
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
save(fullfile(outputFolder, outputFile), 'MRSCont');

% Close any remaining open figures
close all;

end