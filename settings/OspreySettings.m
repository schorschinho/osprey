function [MRSCont] = OspreySettings
%% [MRSCont] = OspreySettings
%   This function initialises default settings and variables required to
%   run Osprey.
%
%   USAGE:
%       [MRSCont] = OspreySettings;
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-02-18)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-02-19: First version of the code.

close all;

%%% 1. REQUIRED INPUT FLAGS %%%
% Initialise MRSCont struct with default settings
MRSCont = struct;
MRSCont.flags.isUnEdited    = 0;
MRSCont.flags.isMEGA        = 0;
MRSCont.flags.isHERMES      = 0;
MRSCont.flags.isHERCULES    = 0;
MRSCont.flags.isPRIAM       = 0;
MRSCont.flags.isMRSI        = 0;
MRSCont.flags.addImages        = 0;
MRSCont.opts.savePDF                = 0;
MRSCont.opts.saveLCM                = 0;
MRSCont.opts.savejMRUI              = 0;
MRSCont.opts.saveNII              = 0;
MRSCont.opts.saveVendor             = 0;
MRSCont.opts.fit.includeMetabs      = {'default'};      % Options: 'default', 'full', custom set       
MRSCont.opts.fit.method             = 'Osprey';         % Options: 'Osprey' (default), 'LCModel'
MRSCont.opts.fit.range              = [0.2 4.2];        % Default: [0.2 4.2]
MRSCont.opts.fit.rangeWater         = [2.0 7.4];        % Default: [2.0 7.4]
MRSCont.opts.fit.style              = 'Separate';       % Options: 'Separate' (Default - will fit DIFF and OFF separately), 'Concatenated' (will fit DIFF and SUM simultaneously), 
MRSCont.opts.fit.bLineKnotSpace     = 0.4;              % Baseline spline knot spacing [ppm]. Default: 0.4.
MRSCont.opts.fit.fitMM              = 1;                % Add MM and lipid basis functions to basis set? Default: 1.
MRSCont.opts.fit.coMM3              = 'none';           % Add co-edited MM3 peak model for GABA editing? Default: none.
MRSCont.opts.fit.FWHMcoMM3          = 14;               % FWHM [Hz] of the co-edited peak Default: 14 Hz.

%%% 2. FIND AND SET PATHS %%%
% Osprey
[settingsFolder,~,~] = fileparts(which('OspreySettings.m'));
allFolders      = strsplit(settingsFolder, filesep);
ospFolder       = strjoin(allFolders(1:end-1), filesep); % parent folder (= Osprey folder)
matlabFolder    = strjoin(allFolders(1:end-2), filesep); % parent-parent folder (usually MATLAB folder)
addpath(genpath(ospFolder));

% SPM
addpath(genpath([matlabFolder filesep 'spm12' filesep]));    % SPM path
% Check if SPM12 is installed
spmversion = fileparts(which('spm'));
if isempty(spmversion)
    warning('SPM not found! Please install SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12) and set the path in OspreySettings.');
elseif strcmpi(spmversion(end-3:end),'spm8')
    warning(['SPM8 detected, but only SPM12 is supported. ' ...
           'Please install SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12) and set the path in OspreySettings.']);
end

%%% 3. INITIALISE MRSCONT FIELDS AND FLAGS %%%
% Set default fields
MRSCont.ospFolder           = ospFolder;
MRSCont.files               = {};
MRSCont.files_mm            = {};
MRSCont.files_ref           = {};
MRSCont.files_w             = {};
MRSCont.flags.hasMM         = 0;
MRSCont.flags.hasRef        = 0;
MRSCont.flags.hasWater      = 0;
% Set default flags
MRSCont.flags.didLCMWrite   = 0;
MRSCont.flags.didjMRUIWrite = 0;
MRSCont.flags.didVendorWrite= 0;
MRSCont.flags.didJob        = 0;
MRSCont.flags.didLoadData   = 0;
MRSCont.flags.didProcess    = 0;
MRSCont.flags.didCoreg      = 0;
MRSCont.flags.didSeg        = 0;
MRSCont.flags.didFit        = 0;
MRSCont.flags.didQuantify   = 0;
MRSCont.flags.didOverview   = 0;


%%% DO NOT EDIT BELOW - DO NOT EDIT BELOW - DO NOT EDIT BELOW - DO NOT EDIT BELOW %%%

%%% 4. INITIALISE ALLOWED FILE TYPES (FOR GUI USE) %%% 
global globalDefaults
globalDefaults.supportedFileTypes.DeIdentify = {'*.sdat','Philips SDAT files (*.sdat)'; ...
    '*.rda','Siemens RDA files (*.rda)'; ...
    '*.dat','Siemens TWIX files (*.dat)'; ...
    '*.7','GE P-files (*.7)'};
globalDefaults.supportedFileTypes.Load = {'*.sdat','Philips SDAT/SPAR files (*.sdat)'; ...
    '*.data','Philips DATA/LIST files (*.data)'; ...
    '*.raw','Philips RAW/SIN/LAB files (*.raw)'; ...
    '*.rda','Siemens RDA files (*.rda)'; ...
    '*.dat','Siemens TWIX files (*.dat)'; ...
    '*.dcm;*.ima','Siemens DICOM files (*.dcm,*.ima)'; ...
    '*.7','GE P-files (*.7)'};
globalDefaults.supportedFileTypes.Nifti = {'*.nii','NIfTI files (*.nii)'};    
globalDefaults.supportedFileTypes.Mat = {'*.mat','MATLAB files (*.mat)'};
%%% DO NOT EDIT - DO NOT EDIT - DO NOT EDIT - DO NOT EDIT %%%

end