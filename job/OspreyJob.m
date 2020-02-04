function [MRSCont] = OspreyJob(jobFile,GUI)
%% [MRSCont] = OspreyJob(jobFile)
%   This function loads the Osprey job defined in jobFile.
%   The job can be submitted in the following formats:
%       - .m
%
%   A valid Osprey job contains four distinct classes of items:
%       1. basic information on the MRS sequence used
%       2. several settings for data handling and modeling
%       3. a list of MRS (and, optionally, structural imaging) data files
%          to be loaded
%       4. an output folder to store the results and exported files
%
%   See the example files in the /exampledata/ folder for details.
%
%   USAGE:
%       [MRSCont] = OspreyJob(jobFile);
%
%   INPUTS:
%       jobFile     = File containing a correct Osprey job definition.
%                     Accepted file formats are .m.
%       GUI         = flag to decide whether plot is used in GUI
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-07-15)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-07-15: First version of the code.

% Close any remaining open figures and parse input arguments
close all;
warning('off','all');
if nargin<2
    GUI = 0;
end

%%% 1. INITIALISE DATA CONTAINER WITH DEFAULT SETTINGS
[MRSCont] = OspreySettings;


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
        outputFolder = jobStruct(1).outputFolder;
    else
        error('Invalid job file! A job file needs to specify an output folder.');
    end
    if isfield(jobStruct,'seqType')
        seqType = jobStruct(1).seqType;
    else
        fprintf('Sequence type is set to unedited (default). Please indicate otherwise in the csv-file or the GUI \n');
        seqType = 'unedited';
    end
    if isfield(jobStruct,'editTarget')
        opts.editTarget = jobStruct(1).editTarget;
    else
        fprintf('Editing target is set to none (default). Please indicate otherwise in the csv-file or the GUI \n');
        opts.editTarget = 'unedited';
    end
    if isfield(jobStruct,'saveLCM')
        opts.saveLCM = jobStruct(1).saveLCM;
    else
        fprintf('LCModel-readable files will be saved (default). Please indicate otherwise in the csv-file or the GUI \n');
        opts.saveLCM = 1;
    end
    if isfield(jobStruct,'saveJMRUI')
        opts.saveJMRUI = jobStruct(1).saveJMRUI;
    else
        fprintf('jMRUI-readable files will be saved (default). Please indicate otherwise in the csv-file or the GUI \n');
        opts.saveJMRUI = 1;
    end
    if isfield(jobStruct,'saveVendor')
        opts.saveVendor = jobStruct(1).saveVendor;
    else
        fprintf('Vendor-specific files (SDAT/SPAR, RDA, P) will be saved (default). Please indicate otherwise in the csv-file or the GUI \n');
        opts.saveVendor = 1;
    end
    if isfield(jobStruct,'method')
        opts.fit.method = jobStruct(1).method;
    else
        fprintf('Fitting algorithm is set to Osprey (default). Please indicate otherwise in the csv-file or the GUI \n');
        opts.fit.method = 'Osprey';
    end
    if isfield(jobStruct,'method')
        opts.fit.style = jobStruct(1).style;
    else
        fprintf('Fitting style is set to Concatenated (default). Please indicate otherwise in the csv-file or the GUI \n');
        opts.fit.style = 'Concatenated';
    end
    if isfield(jobStruct,'lolim_range') && isfield(jobStruct,'uplim_range')
        opts.fit.range = [jobStruct(1).lolim_range jobStruct(1).uplim_range];
    else
        fprintf('Fitting range is set to [0.2 4.2] ppm (default). Please indicate otherwise in the csv-file or the GUI \n');
        opts.fit.range = [0.2 4.2] ;
    end
    if isfield(jobStruct,'lolim_range') && isfield(jobStruct,'uplim_range')
        opts.fit.range = [jobStruct(1).lolim_range jobStruct(1).uplim_range];
    else
        fprintf('Fitting range is set to [0.2 4.2] ppm (default). Please indicate otherwise in the csv-file or the GUI \n');
        opts.fit.range = [0.2 4.2] ;
    end
    if isfield(jobStruct,'KnotSpace')
        opts.fit.bLineKnotSpace = jobStruct(1).KnotSpace;
    else
        fprintf('Baseline knot spacing is set to 0.4 ppm (default). Please indicate otherwise in the csv-file or the GUI \n');
        opts.fit.bLineKnotSpace = 0.4;
    end
    if isfield(jobStruct,'fitMM')
        opts.fit.fitMM = jobStruct(1).fitMM;
    else
        fprintf('Adding macromolecule and lipid basis functions to the fit (default). Please indicate otherwise in the csv-file or the GUI \n');
        opts.fit.fitMM = 1;
    end
end


%%% 4. SAVE SETTINGS & STAT FILE INTO MRSCONT  %%%
switch seqType
    case 'unedited'
        MRSCont.flags.isUnEdited    = 1;
        opts.editTarget             = 'none';
        if strcmp(opts.fit.style, 'Concatenated')
            fprintf('Fitting style was changed to Separate, because this is unedited data. Please indicate otherwise in the csv-file or the GUI \n');
            opts.fit.style = 'Separate';
        end        
    case 'MEGA'
        MRSCont.flags.isMEGA        = 1;
        opts.editTarget             = editTarget;
    case 'HERMES'
        MRSCont.flags.isHERMES      = 1;
        opts.editTarget             = editTarget;
    case 'HERCULES'
        MRSCont.flags.isHERCULES    = 1;
        opts.editTarget             = editTarget;
    otherwise
        error('Invalid job file! seqType must be ''unedited'', ''MEGA'', ''HERMES'', or ''HERCULES''.');
end
MRSCont.opts = opts;
if exist('file_stat','var')
    MRSCont.file_stat = file_stat;
    MRSCont.flags.hasStatfile = 1;
else
    MRSCont.flags.hasStatfile = 0;
end


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


%%% 6. SET UP DEFAULT OSPREY COLORMAP %%%
% Default colormap
colormap.Background     = [255/255 254/255 254/255];
colormap.LightAccent    = [110/255 136/255 164/255];
colormap.Foreground     = [11/255 71/255 111/255];
colormap.Accent         = [11/255 71/255 111/255];
MRSCont.colormap        = colormap;


%%% 7. SET FLAGS %%%
MRSCont.flags.didLoadJob    = 1;
MRSCont.loadedJob           = jobFile;


%%% 8. CHECK IF OUTPUT STRUCTURE ALREADY EXISTS IN OUTPUT FOLDER %%%
[~,jobfilename,jobfileext]  = fileparts(jobFile);
outputFile                  = strrep([jobfilename jobfileext], jobfileext, '.mat');
MRSCont.outputFile          = outputFile;
if ~GUI
    if exist(fullfile(outputFolder, outputFile), 'file') == 2
        fprintf('Your selected output folder ''%s'' already contains an Osprey output structure: \t%s.\n', outputFolder, outputFile);
        fprintf('You are about to load the job: \t%s.\n', jobFile);
        askOverWriteJob = input('Do you want to overwrite the existing job (y/n)? (Warning: This will delete all data in the data container!)   [y]   ','s');
        if isempty(askOverWriteJob)
            askOverWriteJob = 'y';
        end
        if askOverWriteJob=='n' || askOverWriteJob=='N'
             askloadMRSCont = input('Do you want to load the corresponding MRS Container (y/n)? [y]   ','s');
             if isempty(askloadMRSCont)
                askloadMRSCont = 'y';
             end
             if askloadMRSCont=='n' || askloadMRSCont=='N'
                disp('Aborted! No new job loaded.');
                return;
                else if askloadMRSCont=='y' || askloadMRSCont=='y'
                        load(fullfile(outputFolder, outputFile));
                    end
             end
        elseif askOverWriteJob=='y' || askOverWriteJob=='y'
            disp('Continue with loading new job, overwriting existing job.');
        end
    end
else
    opts.Interpreter = 'tex';
    opts.Default = 'Yes';
    if exist(fullfile(outputFolder, outputFile), 'file') == 2
        askOverWriteJob = questdlg(['Your selected output folder already contains an Osprey output structure: ' outputFile '\newline' ...
                                      'You are about to load the job: ', strrep(outputFile,'.mat','.m') '\newline' ...
                                      'Do you want to overwrite the existing job? \newline (Warning: This will delete all data in the data container!)'], ...
        'Load jobFile', 'Yes ','No', opts);
        if strcmp(askOverWriteJob, 'No')
             askloadMRSCont = questdlg( 'Do you want to load the corresponding MRS Container?', ...
        'Load MRS Container', 'Yes ','No','Yes');

             if strcmp(askloadMRSCont, 'No')
                disp('Aborted! No new job loaded.');
                return;
                else if strcmp(askloadMRSCont, 'Yes')
                        load(fullfile(outputFolder, outputFile));
                    end
             end
        elseif strcmp(askOverWriteJob, 'Yes')
            disp('Continue with loading new job, overwriting existing job.');
        end
    
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
