function [MRSCont] = OspreyProcess(MRSCont)
%% [MRSCont] = OspreyProcess(MRSCont)
%   This function pre-processes MRS data from all major vendors.
%   Data is read from the provided input filenames. It is shaped,
%   preprocessed, aligned, etc. according to the type of sequence
%   (un-edited data, MEGA-edited (ON/OFF), HERMES/HERCULES (A/B/C/D),
%   etc.).
%
%   USAGE:
%       MRSCont = OspreyProcess(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
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

% Check that OspreyLoad has been run before
if ~MRSCont.flags.didLoadData
    error('Trying to process data, but raw data has not been loaded yet. Run OspreyLoad first.')
end

% Post-process raw data depending on sequence type
if MRSCont.flags.isUnEdited
    [MRSCont] = osp_processUnEdited(MRSCont);
elseif MRSCont.flags.isMEGA
    [MRSCont] = osp_processMEGA(MRSCont);
elseif MRSCont.flags.isHERMES
    [MRSCont] = osp_processHERMES(MRSCont);
elseif MRSCont.flags.isHERCULES
    % For now, process HERCULES like HERMES data
    [MRSCont] = osp_processHERMES(MRSCont);
else
    error('No flag set for sequence type!');
end

%% Clean up and save
% Set exit flags
MRSCont.flags.didProcess           = 1;

% Save the output structure to the output folder
% Determine output folder
outputFolder    = MRSCont.outputFolder;
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
save(fullfile(outputFolder, outputFile), 'MRSCont');


% Optional: write edited files to LCModel .RAW files
if MRSCont.opts.saveLCM
    [MRSCont] = osp_saveLCM(MRSCont);
end

% Optional: write edited files to jMRUI .txt files
if MRSCont.opts.saveJMRUI
    [MRSCont] = osp_saveJMRUI(MRSCont);
end

% Optional: write edited files to vendor specific format files readable to
% LCModel and jMRUI
% SPAR/SDAT if Philips
% RDA if Siemens
if MRSCont.opts.saveVendor
    [MRSCont] = osp_saveVendor(MRSCont);
end


end