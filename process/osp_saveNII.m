function [MRSCont] = osp_saveNII(MRSCont)
%% [MRSCont] = osp_saveNII(MRSCont)
%   This function writes all processed MRS data to NIfTI-MRS files.
%   
%   One file is produced for unedited MRS data (PRESS, sLASER, etc.), four 
%   files are produced for MEGA-edited data (A, B, sum, difference), and 
%   seven files are produced for HERMES/HERCULES-edited data
%   (A, B, C, D, sum, diff1, diff2).
%
%   If reference scans and short-TE water scans are provided, one file
%   is produced, independent of the type of sequence. In all cases, the
%   water scan will be the sum of all water-unsuppressed scans.
%
%   USAGE:
%       [MRSCont] = osp_saveNII(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-08-06)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)

% Close any remaining open figures
close all;

% Set up saving location
saveDestination = fullfile(MRSCont.outputFolder, 'NIfTIMRS');
outputFunction = @io_writeniimrs;

%% Export files

% Loop over all datasets
for kk = 1:MRSCont.nDatasets
    
    % Set up saving location
    if ~exist(saveDestination,'dir')
        mkdir(saveDestination);
    end
        
    % Write NIfTI-MRS files depending on sequence type
    % Get the input file name
    [path,filename,ext]       = fileparts(MRSCont.files{kk});
    if strcmp(ext, '.gz')
        % If compressed, remove the .nii part as well
        filename = strrep(filename, '.nii', '');
    end
    
    % For batch analysis, get the last three sub-folders (e.g. subject,
    % session, and scan)
    path_split          = regexp(path,filesep,'split');
    if length(path_split) >= 3
        name = [path_split{end-2} '_' path_split{end-1} '_' path_split{end} '_' filename];
    elseif length(path_split) == 2
        name = [path_split{end-1} '_' path_split{end} '_' filename];
    elseif length(path_split) == 1
        name = [path_split{1} '_' filename];
    end
    
    if MRSCont.flags.isUnEdited
        outfile         = fullfile(saveDestination, [name '_A.nii.gz']);
        RF              = outputFunction(MRSCont.processed.A{kk},outfile);
    elseif MRSCont.flags.isMEGA
        outfileA        = fullfile(saveDestination, [name '_A.nii.gz']);
        RF              = outputFunction(MRSCont.processed.A{kk},outfileA);
        outfileB        = fullfile(saveDestination, [name '_B.nii.gz']);
        RF              = outputFunction(MRSCont.processed.B{kk},outfileB);
        outfileDiff1    = fullfile(saveDestination, [name '_DIFF1.nii.gz']);
        RF              = outputFunction(MRSCont.processed.diff1{kk},outfileDiff1);
        outfileSum      = fullfile(saveDestination, [name '_SUM.nii.gz']);
        RF              = outputFunction(MRSCont.processed.sum{kk},outfileSum);
    elseif MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
        outfileA        = fullfile(saveDestination, [name '_A.nii.gz']);
        RF              = outputFunction(MRSCont.processed.A{kk},outfileA);
        outfileB        = fullfile(saveDestination, [name '_B.nii.gz']);
        RF              = outputFunction(MRSCont.processed.B{kk},outfileB);
        outfileC        = fullfile(saveDestination, [name '_C.nii.gz']);
        RF              = outputFunction(MRSCont.processed.C{kk},outfileC);
        outfileD        = fullfile(saveDestination, [name '_D.nii.gz']);
        RF              = outputFunction(MRSCont.processed.D{kk},outfileD);
        outfileDiff1    = fullfile(saveDestination, [name '_DIFF1.nii.gz']);
        RF              = outputFunction(MRSCont.processed.diff1{kk},outfileDiff1);
        outfileDiff2    = fullfile(saveDestination, [name '_DIFF2.nii.gz']);
        RF              = outputFunction(MRSCont.processed.diff2{kk},outfileDiff2);
        outfileSum      = fullfile(saveDestination, [name '_SUM.nii.gz']);
        RF              = outputFunction(MRSCont.processed.sum{kk},outfileSum);
    else
        error('No flag set for sequence type!');
    end
    
    % Check if reference scans exist, if so, write LCM .RAW file
    if MRSCont.flags.hasRef
        % Get the input file name. For GE, the water reference is
        % already contained in the P file.
        if strcmpi(MRSCont.vendor, 'GE')
            [path_ref,filename_ref,ext_ref]   = fileparts(MRSCont.files{kk});
        else
            [path_ref,filename_ref,ext_ref]   = fileparts(MRSCont.files_ref{kk});
        end
        if strcmp(ext_ref, '.gz')
            % If compressed, remove the .nii part as well
            filename_ref = strrep(filename_ref, '.nii', '');
        end
        % For batch analysis, get the last two sub-folders (e.g. site and
        % subject)
        path_ref_split          = regexp(path_ref,filesep,'split');
        if length(path_ref_split) >= 3
            name_ref = [path_ref_split{end-2} '_' path_ref_split{end-1} '_' path_ref_split{end} '_' filename_ref];
        elseif length(path_ref_split) == 2
            name_ref = [path_ref_split{end-1} '_' path_ref_split{end} '_' filename_ref];
        elseif length(path_ref_split) == 1
            name_ref = [path_ref_split{1} '_' filename_ref];
        end
        outfileRef      = fullfile(saveDestination, [name_ref '_REF.nii.gz']);
        RF              = outputFunction(MRSCont.processed.ref{kk},outfileRef);
    end
    
    % Now do the same for the (short-TE) water signal
    if MRSCont.flags.hasWater
        % Get TE and the input file name
        [path_w,filename_w,ext_w]   = fileparts(MRSCont.files_w{kk});
        % For batch analysis, get the last two sub-folders (e.g. site and
        % subject)
        if strcmp(ext_w, '.gz')
            % If compressed, remove the .nii part as well
            filename_w = strrep(filename_w, '.nii', '');
        end
        path_w_split          = regexp(path_w,filesep,'split');
        if length(path_w_split) >= 3
            name_w = [path_w_split{end-2} '_' path_w_split{end-1} '_' path_w_split{end} '_' filename_w];
        elseif length(path_w_split) == 2
            name_w = [path_w_split{end-1} '_' path_w_split{end} '_' filename_w];
        elseif length(path_w_split) == 1
            name_w = [path_w_split{1} '_' filename_w];
        end
        outfileW        = fullfile(saveDestination, [name_w '_W.nii.gz']);
        RF              = outputFunction(MRSCont.processed.w{kk},outfileW);
    end
end

% Set exit flags
MRSCont.flags.didNIIWrite           = 1;

end