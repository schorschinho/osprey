function [MRSCont] = LCGannetProcess(MRSCont)
%% [MRSCont] = LCGannetProcess(MRSCont)
%   This function pre-processes MRS data from all major vendors.
%   Data is read from the provided input filenames. It is shaped,
%   preprocessed, aligned, etc. according to the type of sequence
%   (un-edited data, MEGA-edited (ON/OFF), HERMES/HERCULES (A/B/C/D),
%   etc.).
%
%   USAGE:
%       MRSCont = LCGannetProcess(MRSCont);
%
%   INPUTS:
%       MRSCont     = LCGannet MRS data container.
%
%   OUTPUTS:
%       MRSCont     = LCGannet MRS data container.
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

% Check that LCGannetLoad has been run before
if ~MRSCont.flags.didLoadData
    error('Trying to process data, but raw data has not been loaded yet. Run LCGannetLoad first.')
end

% Post-process raw data depending on sequence type
if MRSCont.flags.isUnEdited
    [MRSCont] = LCG_processUnEdited(MRSCont);
elseif MRSCont.flags.isMEGA
    error('Coming soon!');
    % [MRSCont] = LCG_processMEGA(MRSCont);
elseif MRSCont.flags.isHERMES
    error('Coming soon!');
    % [MRSCont] = LCG_processHERMES(MRSCont);
elseif MRSCont.flags.isHERCULES
    error('Coming soon!');
    %[MRSCont] = LCG_processHERCULES(MRSCont);
else
    error('No flag set for sequence type!');
end

% Start visualization output of load here

%% Clean up and save
% Set exit flags
MRSCont.flags.didProcess           = 1;

% Optional: write edited files to LCModel .RAW files
if MRSCont.opts.saveLCM
    [MRSCont] = LCG_saveLCM(MRSCont);
end

% Optional: write edited files to jMRUI .txt files
if MRSCont.opts.saveJMRUI
    [MRSCont] = LCG_saveJMRUI(MRSCont);
end

end