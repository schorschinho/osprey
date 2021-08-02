function [MRSCont, retMsg] = osp_detDataType(MRSCont)
%% [MRSCont, retMsg] = osp_detDataType(MRSCont)
%   This function determines the MRI vendor and datatype of the filenames
%   provided in the MRSCont.files cells. A warning is flagged if the
%   vendor (and filetype) is not identical between all files.
%
%   USAGE:
%       [MRSCont, retMsg] = osp_detDataType(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%       retMsg      = Return message.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-19)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-02-19: First version of the code.

% Concatenate all data including MM, water and reference scans
files = MRSCont.files;
if isfield(MRSCont, 'files_mm')
    files = [files MRSCont.files_mm];
end
if isfield(MRSCont, 'files_ref')
    files = [files MRSCont.files_ref];
end
if isfield(MRSCont, 'files_w')
    files = [files MRSCont.files_w];
end

% Determine data and vendor for each file
for kk = 1:length(files)
    % Determine whether the provided array entries are files or folders
    whatIs = exist(files{kk});
    switch whatIs
        case 7
            % If the input data is provided as folders, there are two 
            % possible formats: 1) Bruker, 2) DICOM.
            % Let us first check whether the provided folder has a
            % subfolder 'pdata', which would indicate a Bruker dataset:
            dirFolder = dir([files{kk}]);
            ispdata = isempty(find(strcmpi({dirFolder.name}, 'pdata')));
            if ~ispdata
                % 1) BRUKER
                buffer.vendor{kk}       = 'Bruker';
                buffer.datatype{kk}     = 'fid';
            else
                % 2) DICOM
                % check that all files in this folder have the
                % same type
                filesInFolder = dirFolder(~[dirFolder.isdir]);
                filesInFolder = filesInFolder(~ismember({filesInFolder.name}, {'.','..','.DS_Store'}));
                hidden = logical(ones(1,length(filesInFolder)));
                for jj = 1:length(filesInFolder)
                    if strcmp(filesInFolder(jj).name(1),'.')
                        hidden(jj) = 0;
                    end
                end
                filesInFolder = filesInFolder(hidden);%delete hidden files
                for rr = 1:length(filesInFolder)
                    [~,~,ext_fold{rr}] = fileparts(filesInFolder(rr).name);
                end
                % If not, throw an error
                if length(unique(ext_fold)) > 1
                    retMsg = sprintf('Error during loading. Folder %s contains data in more than one file format.\n', files{kk});
                    sprintf('Error during loading. Folder %s contains data in more than one file format.\n', files{kk});
                end
                % If all files have the same extension, pick the first one to
                % determine the format
                fileToDet = fullfile(files{kk}, filesInFolder(1).name);
            end
            
        case 2
            % If files, just take the filename.
            fileToDet = files{kk};
        case 0
            % If it doesn't exist, throw an error
            sprintf('Error during loading. File or folder %s does not exist. Check the job file!\n', files{kk});
            error('Error during loading. File or folder %s does not exist. Check the job file!\n', files{kk});
    end
    
    % Parse file extension
    if exist('fileToDet')
        [~,~,ext] = fileparts(fileToDet);
        if strcmpi(ext,'.7')
            buffer.vendor{kk}       = 'GE';
            buffer.datatype{kk}     = 'P';
            MRSCont.flags.hasRef = 1;
        elseif strcmpi(ext,'.RDA')
            buffer.vendor{kk}       = 'Siemens';
            buffer.datatype{kk}     = 'RDA';
        elseif strcmpi(ext,'.IMA') || strcmpi(ext,'.DCM') || strcmpi(ext,'')
            buffer.vendor{kk}       = 'Siemens';
            buffer.datatype{kk}     = 'DICOM';
        elseif strcmpi(ext,'.DAT')
            buffer.vendor{kk}       = 'Siemens';
            buffer.datatype{kk}     = 'TWIX';
        elseif strcmpi(ext,'.RAW')
            buffer.vendor{kk}       = 'LCModel';
            buffer.datatype{kk}     = 'RAW';
        elseif strcmpi(ext,'.LAB')
            buffer.vendor{kk}       = 'Philips';
            buffer.datatype{kk}     = 'RAW';
        elseif strcmpi(ext,'.SDAT')
            buffer.vendor{kk}       = 'Philips';
            buffer.datatype{kk}     = 'SDAT';
        elseif strcmpi(ext,'.DATA')
            buffer.vendor{kk}       = 'Philips';
            buffer.datatype{kk}     = 'DATA';
        elseif strcmpi(ext,'.gz') || strcmpi(ext,'.nii')
            % For now, leave the vendor field empty; we'll fill it while
            % loading the actual data.
            buffer.vendor{kk}       = '';
            buffer.datatype{kk}     = 'NIfTI-MRS';
        else
            retMsg = 'Unrecognized datatype. Filenames need to end in .7 .SDAT .DATA .RAW .RDA .IMA .DCM .DAT .GZ or .NII.GZ!';
            fprintf(retMsg);
        end
    end
end

% Make sure that all provided files share the same filetype
seq_params = {'vendor','datatype'};
for pp = 1:length(seq_params)
    unique_params = unique(buffer.(seq_params{pp}));
    if length(unique_params) > 1
        retMsg = 'WARNING! One or more datatypes or vendors are not the same across all input files.';
        sprintf(retMsg);
    else
        MRSCont.(seq_params{pp}) = unique_params{1};
    end
end

% Print final success message
if ~exist('retMsg','var')
    retMsg = sprintf('All provided datafiles are of the %s %s format.', MRSCont.vendor, MRSCont.datatype);
    sprintf('All provided datafiles are of the %s %s format.\n', MRSCont.vendor, MRSCont.datatype);
end

end