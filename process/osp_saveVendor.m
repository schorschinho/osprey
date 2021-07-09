function [MRSCont] = osp_saveVendor(MRSCont)
%% [MRSCont] = osp_saveVendor(MRSCont)
%   This function writes all MRS data loaded by OspreyLoad to separate
%   vendor-specific raw data file formats:
%
%   -   If the input format is by Philips (SDAT/SPAR, DATA/LIST, SIN/LAB/RAW),
%       then this function will create SDAT/SPAR files.
%   -   If the input format is by Siemens (RDA, TWIX, DCM, IMA), then this 
%       function will create RDA files.
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
%       [MRSCont] = osp_saveVendor(MRSCont);
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
%
%   HISTORY:
%       2019-08-06: First version of the code.

% Close any remaining open figures
close all;

% Set up saving location
switch MRSCont.vendor
    case 'Philips'
        saveDestination = fullfile(MRSCont.outputFolder, 'SDATFiles');
        outputFunction = @(in,outfile) io_writesdatspar(in,outfile);
    case 'Siemens'
        saveDestination = fullfile(MRSCont.outputFolder, 'RDAFiles');
        outputFunction = @(in,outfile) io_writerda(in,outfile);
    case 'GE'
        % coming soon
        % outputFolder = 'PFiles';
        return
    otherwise
        warning('Saving processed data in %s format is currently not supported. Please contact the Gannet team (gabamrs@gmail.com).', MRSCont.vendor);
        return
end


%% Export files

% Loop over all datasets
for kk = 1:MRSCont.nDatasets
    
    % Set up saving location
    if ~exist(saveDestination,'dir')
        mkdir(saveDestination);
    end
        
    % Write vendor-specific files depending on sequence type
    % Get the input file name
    [path,filename,~]       = fileparts(MRSCont.files{kk});
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
        outfile         = fullfile(saveDestination, [name '_A']);
        RF              = outputFunction(MRSCont.processed.A{kk},outfile);
    elseif MRSCont.flags.isMEGA
        outfileA        = fullfile(saveDestination, [name '_A']);
        RF              = outputFunction(MRSCont.processed.A{kk},outfileA);
        outfileB        = fullfile(saveDestination, [name '_B']);
        RF              = outputFunction(MRSCont.processed.B{kk},outfileB);
        outfileDiff1    = fullfile(saveDestination, [name '_DIFF1']);
        RF              = outputFunction(MRSCont.processed.diff1{kk},outfileDiff1);
        outfileSum      = fullfile(saveDestination, [name '_SUM']);
        RF              = outputFunction(MRSCont.processed.sum{kk},outfileSum);
    elseif MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
        outfileA        = fullfile(saveDestination, [name '_A']);
        RF              = outputFunction(MRSCont.processed.A{kk},outfileA);
        outfileB        = fullfile(saveDestination, [name '_B']);
        RF              = outputFunction(MRSCont.processed.B{kk},outfileB);
        outfileC        = fullfile(saveDestination, [name '_C']);
        RF              = outputFunction(MRSCont.processed.C{kk},outfileC);
        outfileD        = fullfile(saveDestination, [name '_D']);
        RF              = outputFunction(MRSCont.processed.D{kk},outfileD);
        outfileDiff1    = fullfile(saveDestination, [name '_DIFF1']);
        RF              = outputFunction(MRSCont.processed.diff1{kk},outfileDiff1);
        outfileDiff2    = fullfile(saveDestination, [name '_DIFF2']);
        RF              = outputFunction(MRSCont.processed.diff2{kk},outfileDiff2);
        outfileSum      = fullfile(saveDestination, [name '_SUM']);
        RF              = outputFunction(MRSCont.processed.sum{kk},outfileSum);
    else
        error('No flag set for sequence type!');
    end
    
    % Check if reference scans exist, if so, write LCM .RAW file
    if MRSCont.flags.hasRef
        % Get the input file name. For GE, the water reference is
        % already contained in the P file.
        if strcmpi(MRSCont.vendor, 'GE')
            [path_ref,filename_ref,~]   = fileparts(MRSCont.files{kk});
        else
            [path_ref,filename_ref,~]   = fileparts(MRSCont.files_ref{kk});
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
        outfileRef      = fullfile(saveDestination, [name_ref '_REF']);
        RF              = outputFunction(MRSCont.processed.ref{kk},outfileRef);
    end
    
    % Now do the same for the (short-TE) water signal
    if MRSCont.flags.hasWater
        % Get TE and the input file name
        [path_w,filename_w,~]   = fileparts(MRSCont.files_w{kk});
        % For batch analysis, get the last two sub-folders (e.g. site and
        % subject)
        path_w_split          = regexp(path_w,filesep,'split');
        if length(path_w_split) >= 3
            name_w = [path_w_split{end-2} '_' path_w_split{end-1} '_' path_w_split{end} '_' filename_w];
        elseif length(path_w_split) == 2
            name_w = [path_w_split{end-1} '_' path_w_split{end} '_' filename_w];
        elseif length(path_w_split) == 1
            name_w = [path_w_split{1} '_' filename_w];
        end
        outfileW        = fullfile(saveDestination, [name_w '_W']);
        RF              = outputFunction(MRSCont.processed.w{kk},outfileW);
    end
end

% Set exit flags
MRSCont.flags.didVendorWrite           = 1;

end