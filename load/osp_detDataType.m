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
            % If folders, then check that all files in this folder have the
            % same type
            dirFolder = dir([files{1}]);
            filesInFolder = dirFolder(~[dirFolder.isdir]);
            for rr = 1:length(filesInFolder)
                [~,~,ext{rr}] = fileparts(filesInFolder(rr).name);
            end
            % If not, throw an error
            if length(unique(ext)) > 1
                retMsg = sprintf('Error during loading. Folder %s contains data in more than one file format.\n', files{kk});
                sprintf('Error during loading. Folder %s contains data in more than one file format.\n', files{kk});
            end
            % If all files have the same extension, pick the first one to
            % determine the format
            fileToDet = [files{kk} filesInFolder(1).name];
        case 2
            % If files, just take the filename.
            fileToDet = files{kk};
        case 0
            % If it doesn't exist, throw an error
            sprintf('Error during loading. File or folder %s does not exist. Check the job file!\n', files{kk});
            error('Error during loading. File or folder %s does not exist. Check the job file!\n', files{kk});
    end
    last2char   = fileToDet((end-1):end);
    last3char   = fileToDet((end-2):end);
    last4char   = fileToDet((end-3):end);
    if strcmpi(last2char,'.7')
        buffer.vendor{kk}       = 'GE';
        buffer.datatype{kk}     = 'P';
        MRSCont.flags.hasRef = 1;
    elseif strcmpi(last3char,'RDA')
        buffer.vendor{kk}       = 'Siemens';
        buffer.datatype{kk}     = 'RDA';
    elseif strcmpi(last3char,'IMA') || strcmpi(last3char,'DCM')
        buffer.vendor{kk}       = 'Siemens';
        buffer.datatype{kk}     = 'DICOM';
    elseif strcmpi(last4char,'.DAT')
        buffer.vendor{kk}       = 'Siemens';
        buffer.datatype{kk}     = 'TWIX';
    elseif strcmpi(last3char,'RAW')
        buffer.vendor{kk}       = 'Philips';
        buffer.datatype{kk}     = 'RAW';
    elseif strcmpi(last4char,'SDAT')
        buffer.vendor{kk}       = 'Philips';
        buffer.datatype{kk}     = 'SDAT';
    elseif strcmpi(last4char,'DATA')
        buffer.vendor{kk}       = 'Philips';
        buffer.datatype{kk}     = 'DATA';
    else
        retMsg = 'Unrecognized datatype. All filenames need to end .7 .SDAT .DATA .RAW .RDA .IMA .DCM or .DAT';
        sprintf(fileID,retMsg);
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