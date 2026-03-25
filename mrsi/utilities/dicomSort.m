function dicomSort(studyPath,varargin)
% dicomSort is a recursive dicom sorting tool.
%   dicomsort(input) recursively sorts all dicom files in a study folder
%   regardless of directory structure or hierarchy.
%
%   Syntax:
%   dicomsort(input)
%   dicomsort(input,output)
%   dicomsort(__,Name,Value)
%
%   Description:
%   dicomSort(input) sorts all dicom files in a study regardless of folder
%   hierarchy
%
%   dicomsort(input,'output',dir) sorts all dicom files in a subject folder
%   and outputs them to a specified folder
%
%   dicomsort(__,Name,Value) uses additional name-value pairs to customize
%   sorting
%
%   Name-Value Pair Arguments:
%   'output' --> output path
%       Define an output directory for sorting
%       Type: char | string
%
%   'preserve' --> true (default) | false
%       Specifies whether to preserve original files or directories.
%       Setting to true deletes all original files
%       Type: logical
%
%   'compression' --> type of compression
%       Specify if compression applied on old files. Recommended if
%       preserve is set to 'false' to prevent loss of data
%       Options: 'none' (default), 'zip', 'tar', 'gzip'
%       Type: char | string
%
%   'prefix' --> prefix before subject
%       String to attach prior to subject name during creation of subject
%       folder
%       Type: char | string
%
%   'suffix' --> suffix after subject
%       String to attach after subject name during creation of subject
%       folder
%
%   Example:
%   dicomsort('~/example_study/data');
%       Sort all dicom files in ~/example_study/data
%
%   dicomsort('~/example_study/data','output','~/example_study/sorted');
%       Sort all dicom files in ~/example_study/data and output sorted data
%       into ~/example_study/sorted
%
%   dicomsort(~/example_study/data,'output','~/example_study/sorted',...
%       'preserve',false','compression','tar');
%       Sort all dicom files in ~/example_study/data and output sorted data
%       into ~/example_study/sorted, whilst tar compressing original data
%       directory hierarchy and removing original files
%
%   Author: Siddhartha Dhiman
%   Email: dhiman@musc.edu
%   First created on 01/28/2019 using MATLAB 2018b
%   Last modified on 01/31/2019 using MATLAB 2018b
%
%   SEE ALSO DIR DICOMINFO COPYFILE MOVEFILE

warning off;
%% Parse Inputs
defaultComp = 'none';
defaultPreserve = true;
expectedLog = {'yes','no'};
expectedComp = {'none','zip','tar','gzip'};
p = inputParser;
addRequired(p,'studyPath',@isstr);
addOptional(p,'output',@isstr);
addParameter(p,'preserve',defaultPreserve,@islogical);
addParameter(p,'compression',defaultComp,...
    @(s) any(validatestring(s,expectedComp)));
addOptional(p,'prefix',@isstring);
addOptional(p,'suffix',@isstring);

parse(p,studyPath,varargin{:});

%% Perform Tests
%   Check whether input path exists
if ~isfolder(studyPath)
    error('Input path not found, ensure it exists');
elseif isstr(p.Results.output)
    outPath = p.Results.output;
    if isdir(p.Results.output)
        ;
    else
        mkdir(p.Results.output)
    end
elseif ~isstr(p.Results.output)
    outPath = studyPath;
end

%% Throw Errors during Impossible Cases and Summarize Results
if ~isstr(p.Results.output) && p.Results.preserve == 0
    error(sprintf('Cannot preserve original folder structure when "output" is undefined.\nPlease define an "output" directory or disable "preserve"'));
elseif ~isstr(p.Results.output) && ~strcmp(p.Results.compression,'none')
    error('Study directory cannot be preserved after altering it, define an output directory to enable compression.');
else
    ;
end

%% Check and Decompress Files if Present
studyDirFolders = dir(studyPath);
compressedFiles = vertcat(dir(fullfile(studyPath,'**/*.zip')),...
    dir(fullfile(studyPath,'**/*.tar')),...
    dir(fullfile(studyPath,'**/*.gz')));

if length(compressedFiles) >= 1
    parfor i = 1:length(compressedFiles)
        %   Check file extension
        [~,name,ext] = fileparts(fullfile(compressedFiles(i).folder,...
            compressedFiles(i).name));
        mkdir(fullfile(compressedFiles(i).folder,name));
        try
            if ext == '.zip';
                unzip(fullfile(compressedFiles(i).folder,...
                    compressedFiles(i).name),...
                    fullfile(compressedFiles(i).folder));
            elseif ext == '.tar'
                untar(fullfile(compressedFiles(i).folder,...
                    compressedFiles(i).name),...
                    fullfile(compressedFiles(i).folder));
            elseif ext == '.gz'
                gunzip(fullfile(compressedFiles(i).folder,...
                    compressedFiles(i).name),...
                    fullfile(compressedFiles(i).folder));
            else
            end
        catch
            fprintf('Corruption detected: skipping uncompression of %s',name);
            continue
        end
    end
end

%% Tunable Function Variables
studyDir = vertcat(dir(fullfile(studyPath,'**/*')),...
    dir(fullfile(studyPath,'**/*.dcm'))); %   Recursive directory listing

rmPattern = {'.','.DS_Store'};   %   Remove files beginning with

%% Clean up Main Study Dir Listing
rmIdx = zeros(1,length(studyDirFolders));
for i = 1:length(studyDirFolders)
    if any(startsWith(studyDirFolders(i).name,'.'));
        rmIdx(i) = 1;
        
    else
        %   If nothing found, don't mark for deletion
        rmIdx(i) = 0;
    end
end
studyDirFolders(rmIdx ~= 0) = [];
%   Create folder listing for compression
for i = 1:length(studyDirFolders)
    compPaths{i} = fullfile(studyDirFolders(i).folder,...
        studyDirFolders(i).name);
end
nComp = length(compPaths);

%% Clean-up Dicom Directory Listing
rmIdx = zeros(1,length(studyDir));
for i = 1:length(studyDir)
    %  Check for directories
    if studyDir(i).isdir == 1
        rmIdx(i) = 1;
        
        %   Check for files starting with'.'
    elseif any(startsWith(studyDir(i).name,'.'));
        rmIdx(i) = 1;
        
    else
        %   If nothing found, don't mark for deletion
        rmIdx(i) = 0;
    end
end
studyDir(rmIdx ~= 0) = [];   %   Apply deletion filter
nFiles = length(studyDir);
fprintf('Found %d dicom files\n',nFiles);

%% Read Parallel Pool Properties
curCluster = parcluster('local');
%   Initialize empty vector same size as number of worrkers for ETA
%   caclulation
etaVec = zeros(nFiles,1);
%   Create exponential moving average object
% movAvg = dsp.MovingAverage('Method','Exponential weighting',...
%     'WindowLength',curCluster.NumWorkers);

%%   Initialize Parallel Data Queue
parQ = parallel.pool.DataQueue;
%   Initialize progress waitbar
parWaitBar = waitbar(0,'Initializing sorting algorithm...',...
    'Name','Sorting dicoms');
%   After receiving new data, update_progress() will be called
fprintf('Sorting dicoms...');
afterEach(parQ,@parProgress);
n_completed = 0;

%% Sort Dicom Files
%   Run in parent parfor for speed
preSorted = cell(nFiles,1);
j = 1;
parfor i = 1:nFiles
    %   Test whether file is readable by dicom. If not, move on to next
    %   file instead of throwing an error
    tic;
    try
        tmp = dicominfo(fullfile(studyDir(i).folder,studyDir(i).name));
    catch
        continue
    end
    
    %   Check whether file contains a '.dcm' extension. If it does, make no
    %   changes. If not, append '.dcm' at the end
    if ~contains(studyDir(i).name,'.dcm')
        newName = [studyDir(i).name '.dcm'];
    else
        newName = studyDir(i).name;
    end
    
    %   Check for prefix flag and append if it does. Skip otherwise
    if isstr(p.Results.prefix)
        tmp.PatientID = [p.Results.prefix tmp.PatientID];
    else
        tmp.PatientID = tmp.PatientID;
    end
    
    %   Check for suffix flag and append if it does. Skip otherwise
    if isstr(p.Results.suffix)
        tmp.PatientID = [tmp.PatientID p.Results.suffix];
    else
        tmp.PatientID = tmp.PatientID;
    end
    
    %   Define copy/move folder and replace all '.' with -_-
    if ~isempty(tmp.ProtocolName)
    defineFolder = fullfile(outPath,tmp.PatientID,...
        sprintf('%0.2d_%s',tmp.SeriesNumber,tmp.ProtocolName));
    else
        defineFolder = fullfile(outPath,tmp.PatientID,...
        sprintf('%0.2d',tmp.SeriesNumber));
    end
    defineFolder = strrep(defineFolder,'.','_');
    
    %   Check if a directory 'PatientID/Protocol' exists. If it does, do
    %   nothing. Otherwise, make directry
    if ~exist(defineFolder,'dir')
        mkdir(defineFolder);
    else
        ;
    end
    
    %   Check whether current file is already sorted. If it is, skip that
    %   file and take note of patient ID
    if contains(tmp.Filename,defineFolder);
        preSorted{i} = tmp.PatientID;
        continue;
    else
        ;
    end
    
    %   Initiate file copy if another output path present
    if isstr(p.Results.output)
        copyfile(tmp.Filename,fullfile(defineFolder,newName));
        %   Move them otherwise
    else
        movefile(tmp.Filename,fullfile(defineFolder,newName));
    end
    
    %   Update parallel process progress
    etaVec(i) = toc;
    send(parQ,i);
end
fprintf('\n');

%   Check compression flag and act accordingly
if strcmp(p.Results.compression,'none');
    fprintf('Skipping compression\n');
elseif strcmp(p.Results.compression,'zip')
    fprintf('Zipping files...');
    zip(fullfile(studyPath,'study_original_files.zip'),compPaths);
    fprintf('saved as %s',fullfile(studyPath,'study_original_files.zip\n'));
elseif strcmp(p.Results.compression,'tar')
    fprintf('Tarring files...');
    tar(fullfile(studyPath,'study_original_files.tar'),compPaths);
    fprintf('saved as %s',fullfile(studyPath,'study_original_files.tar\n'));
elseif strcmp(p.Results.compression,'gzip')
    fprintf('Gunziping files...');
    gzip(compPaths,studyPath);
    fprintf('saved as %s',fullfile(studyPath,'gzip','study_original_files.gz\n'));
else
    fprintf('Not sure what your compression options are');
end

%   From the list of PatientIDs that were already sorted, create a unique
%   list of IDs to that acts as a do-not-delete list. Check whether output
%   diretory has been not been defined. If both checks succeed, delete
%   files not on the DND.
if exist('preSorted','var') && ~isstr(p.Results.output)
    preSorted = unique(preSorted(~cellfun('isempty',preSorted)));
    %   Remove old file and structure if the output directory was unspecified
    fprintf('Not preserving files, removing folders:\n');
    for i = 1:nComp
        %   Check if file to be deleted is a member of do-no-delete list.
        %   If it's not, delete.
        if ~ismember(studyDirFolders(i).name,preSorted)
            fprintf('%d/%d:    %s\n',i,nComp,...
                fullfile(studyDirFolders(i).folder,studyDirFolders(i).name));
            if studyDirFolders(i).isdir
                rmdir(fullfile(studyDirFolders(i).folder,studyDirFolders(i).name),'s');
            else
                delete(fullfile(studyDirFolders(i).folder,studyDirFolders(i).name));
            end
        else
            ;
        end
    end
    %   Then check if an output directory is undefined and no DND list
    %   exists. If checks pass, delete all old files
elseif ~exist('preSorted','var') && ~isstr(p.Results.output)
    for i = 1:nComp
        fprintf('%d/%d:    %s\n',i,nComp,...
            fullfile(studyDirFolders(i).folder,studyDirFolders(i).name));
        rmdir(fullfile(studyDirFolders(i).folder,studyDirFolders(i).name),'s');
    end
else
    fprintf('Preserving orignal files\n');
end

%% PARFOR Progress Calculation
    function parProgress(~)
        if ~exist('n_completed','var')
            n_completed = 0;
        else
            n_completed = n_completed + 1;
        end
        %   Calculate percentage progress
        parPercentage = n_completed/nFiles*100;
        %   Calculate ETA
        % meanCurTime = movAvg(nonzeros(etaVec));
        meanCurTime=1;
        eta = meanCurTime(end) * (nFiles - n_completed);
        %   Update waitbar
        waitbar(n_completed/nFiles,parWaitBar,...
            sprintf('%0.1f%% completed\nETA:%0.0f sec (%0.1f min)',...
            parPercentage,eta,eta/60));
    end
end