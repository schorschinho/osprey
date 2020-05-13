function [MRSCont] = osp_saveJMRUI(MRSCont)
%% [MRSCont] = osp_saveJMRUI(MRSCont)
%   This function writes all MRS data loaded by OspreyLoad to separate
%   jMRUI-readable .TXT files.
%   
%   One .TXT file is produced for unedited MRS data
%   (PRESS, sLASER, etc.), four .TXT files are produced for MEGA-edited
%   data (A, B, sum, difference), and seven .TXT files are produced for
%   HERMES/HERCULES-edited data (A, B, C, D, sum, diff1, diff2).
%
%   If reference scans and short-TE water scans are provided, one .TXT file
%   is produced, independent of the type of sequence. In all cases, the
%   water scan will be the sum of all water-unsuppressed scans.
%
%   USAGE:
%       [MRSCont] = osp_saveJMRUI(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-07-23)
%       goeltzs1@jhmi.edu
%   
%   CREDITS:    
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-07-23: First version of the code.

% Close any remaining open figures
close all;

% Set up saving location
saveDestination = fullfile(MRSCont.outputFolder, 'jMRUIFiles');
if ~exist(saveDestination,'dir')
    mkdir(saveDestination);
end

%% Export files

% Loop over all datasets
for kk = 1:MRSCont.nDatasets
    
    % Write jMRUI .TXT files depending on sequence type
    % Get TE and the input file name
    te                  = MRSCont.processed.A{kk}.te;
    [path,filename,~]   = fileparts(MRSCont.files{kk});
    % For batch analysis, get the last two sub-folders (e.g. site and
    % subject)
    path_split          = regexp(path,filesep,'split');
    if length(path_split) > 2
        name = [path_split{end-1} '_' path_split{end} '_' filename];
    end
    if MRSCont.flags.isUnEdited
        outfile         = fullfile(saveDestination, [name '_jMRUI_A.TXT']);
        RF              = io_writejmrui(MRSCont.processed.A{kk},outfile);
    elseif MRSCont.flags.isMEGA
        outfileA        = fullfile(saveDestination, [name '_jMRUI_A.TXT']);
        RF              = io_writejmrui(MRSCont.processed.A{kk},outfileA);
        outfileB        = fullfile(saveDestination, [name '_jMRUI_B.TXT']);
        RF              = io_writejmrui(MRSCont.processed.B{kk},outfileB);
        outfileDiff1    = fullfile(saveDestination, [name '_jMRUI_DIFF1.TXT']);
        RF              = io_writejmrui(MRSCont.processed.diff1{kk},outfileDiff1);
        outfileSum      = fullfile(saveDestination, [name '_jMRUI_SUM.TXT']);
        RF              = io_writejmrui(MRSCont.processed.sum{kk},outfileSum);
    elseif MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
        outfileA        = fullfile(saveDestination, [name '_jMRUI_A.TXT']);
        RF              = io_writejmrui(MRSCont.processed.A{kk},outfileA);
        outfileB        = fullfile(saveDestination, [name '_jMRUI_B.TXT']);
        RF              = io_writejmrui(MRSCont.processed.B{kk},outfileB);
        outfileC        = fullfile(saveDestination, [name '_jMRUI_C.TXT']);
        RF              = io_writejmrui(MRSCont.processed.C{kk},outfileC);
        outfileD        = fullfile(saveDestination, [name '_jMRUI_D.TXT']);
        RF              = io_writejmrui(MRSCont.processed.D{kk},outfileD);
        outfileDiff1    = fullfile(saveDestination, [name '_jMRUI_DIFF1.TXT']);
        RF              = io_writejmrui(MRSCont.processed.diff1{kk},outfileDiff1);
        outfileDiff2    = fullfile(saveDestination, [name '_jMRUI_DIFF2.TXT']);
        RF              = io_writejmrui(MRSCont.processed.diff2{kk},outfileDiff2);
        outfileSum      = fullfile(saveDestination, [name '_jMRUI_SUM.TXT']);
        RF              = io_writejmrui(MRSCont.processed.sum{kk},outfileSum);
    else
        error('No flag set for sequence type!');
    end
    
    % Check if reference scans exist, if so, write jMRUI .TXT file
    if MRSCont.flags.hasRef
        % Get TE and the input file name. For GE, the water reference is
        % already contained in the P file.
        if strcmpi(MRSCont.vendor, 'GE')
            te_ref                      = MRSCont.processed.A{kk}.te;
            [path_ref,filename_ref,~]   = fileparts(MRSCont.files{kk});
        else
            te_ref                      = MRSCont.processed.ref{kk}.te;
            [path_ref,filename_ref,~]   = fileparts(MRSCont.files_ref{kk});
        end
        % For batch analysis, get the last two sub-folders (e.g. site and
        % subject)
        path_ref_split          = regexp(path_ref,filesep,'split');
        if length(path_ref_split) > 2
            name_ref = [path_ref_split{end-1} '_' path_ref_split{end} '_' filename_ref];
        end
        outfileRef      = fullfile(saveDestination, [name_ref '_jMRUI_REF.TXT']);
        RF              = io_writejmrui(MRSCont.processed.ref{kk},outfileRef);
    end
    
    % Now do the same for the (short-TE) water signal
    if MRSCont.flags.hasWater
        % Get TE and the input file name
        te_w                = MRSCont.processed.w{kk}.te;
        [path_w,filename_w,~]   = fileparts(MRSCont.files_w{kk});
        % For batch analysis, get the last two sub-folders (e.g. site and
        % subject)
        path_w_split          = regexp(path_w,filesep,'split');
        if length(path_w_split) > 2
            name_w = [path_w_split{end-1} '_' path_w_split{end} '_' filename_w];
        end
        outfileW        = fullfile(saveDestination, [name_w '_jMRUI_W.TXT']);
        RF              = io_writejmrui(MRSCont.processed.w{kk},outfileW);
    end
end

% Set exit flags
MRSCont.flags.didjMRUIWrite           = 1;

end