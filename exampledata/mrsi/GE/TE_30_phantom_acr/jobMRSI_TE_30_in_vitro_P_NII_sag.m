%% jobSDAT.m
%   This function describes an Osprey job defined in a MATLAB script.
%
%   A valid Osprey job contains four distinct classes of items:
%       1. basic information on the MRS sequence used
%       2. several settings for data handling and modeling
%       3. a list of MRS (and, optionally, structural imaging) data files
%          to be loaded
%       4. an output folder to store the results and exported files
%
%   The list of MRS and structural imaging files is provided in the form of
%   cell arrays. They can simply be provided explicitly, or from a more
%   complex script that automatically determines file names from a given
%   folder structure.
%
%   Osprey distinguishes between four sets of data:
%       - metabolite (water-suppressed) data
%           (MANDATORY)
%           Defined in cell array "files"
%       - water reference data acquired with the SAME sequence as the
%           metabolite data, just without water suppression RF pulses. This
%           data is used to determine complex coil combination
%           coefficients, and perform eddy current correction.
%           (OPTIONAL)
%           Defined in cell array "files_ref"
%       - additional water data used for water-scaled quantification,
%           usually from short-TE acquisitions due to reduced T2-weighting
%           (OPTIONAL)
%           Defined in cell array "files_w"
%       - Structural image data used for co-registration and tissue class
%           segmentation (usually a T1 MPRAGE). These files need to be
%           provided in the NIfTI format (*.nii) or, for GE data, as a
%           folder containing DICOM Files (*.dcm).
%           (OPTIONAL)
%           Defined in cell array "files_nii"
%
%   Files in the formats
%       - .7 (GE)
%       - .SDAT, .DATA/.LIST, .RAW/.SIN/.LAB (Philips)
%       - .DAT (Siemens)
%   usually contain all of the acquired data in a single file per scan. GE
%   systems store water reference data in the same .7 file, so there is no
%   need to specify it separately under files_ref.
%
%   Files in the formats
%       - .DCM (any)
%       - .IMA, .RDA (Siemens)
%   may contain separate files for each average. Instead of providing
%   individual file names, please specify folders. Metabolite data, water
%   reference data, and water data need to be located in separate folders.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-07-15)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-07-15: First version of the code.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 1. SPECIFY SEQUENCE INFORMATION %%%

% Specify sequence type
seqType = 'unedited';           % OPTIONS:    - 'unedited' (default)
                                %             - 'MEGA'
                                %             - 'HERMES'
                                %             - 'HERCULES'
                                % Specify Multi voxel type (optional)
MultiVoxel = 'MRSI';           % OPTIONS:    - 'PRIAM' (default)

dataScenario = 'phantom';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2. SPECIFY DATA HANDLING AND MODELING OPTIONS %%%
% Save LCModel-exportable files for each spectrum?
opts.saveLCM                = 0;                % OPTIONS:    - 0 (no, default)
                                                %             - 1 (yes)
% Save jMRUI-exportable files for each spectrum?
opts.saveJMRUI              = 0;                % OPTIONS:    - 0 (no, default)
                                                %             - 1 (yes)

% Save processed spectra in vendor-specific format (SDAT/SPAR, RDA, P)?
opts.saveVendor             = 0;                % OPTIONS:    - 0 (no, default)
                                                %             - 1 (yes)

% Choose the fitting algorithm
opts.fit.method             = 'Osprey';       % OPTIONS:    - 'Osprey' (default)
                                                %           - 'AQSES' (planned)
                                                %           - 'LCModel' (planned)
                                                %           - 'TARQUIN' (planned)

% Choose the fitting style for difference-edited datasets (MEGA, HERMES, HERCULES)
% (only available for the Osprey fitting method)
opts.fit.style              = 'Separate';   % OPTIONS:  - 'Concatenated' (default) - will fit DIFF and SUM simultaneously)
                                                %           - 'Separate' - will fit DIFF and OFF separately

% Determine fitting range (in ppm) for the metabolite and water spectra
opts.fit.range              = [0.5 4.0];        % [ppm] Default: [0.2 4.2]
opts.fit.rangeWater         = [2.0 7.4];        % [ppm] Default: [2.0 7.4]

% Determine the baseline knot spacing (in ppm) for the metabolite spectra
opts.fit.bLineKnotSpace     = 0.4;              % [ppm] Default: 0.4.

% Add macromolecule and lipid basis functions to the fit?
opts.fit.fitMM              = 1;                % OPTIONS:    - 0 (no)
                                                %             - 1 (yes, default)


% Now we define MRSI specific processing options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do you want interactive feedback during the processing?
opts.MRSI.interactive = 0;

% You can add voxel shifts and flips to the nifti output if needed
opts.MRSI.nii_shifts = [0.5 0.5 0.5];
opts.MRSI.nii_flip.cc = 1;

% Export into single nii-file
opts.MRSI.pseudo3D = 1;

% This will generate amplitude integral maps for quick inspection. You can
% define different regions and spectra to be used.

opts.MRSI.quickMapsList{1}.specs = {'raw'};
opts.MRSI.quickMapsList{1}.target = 'raw';
opts.MRSI.quickMapsList{1}.names.raw = {'tNAA','tCr','tCho','tLip','H2O'};
opts.MRSI.quickMapsList{1}.limits.raw = [1.95, 2.1; 2.95, 3.12; 3.12, 3.25;0, 1.95;4.1, 6;];
opts.MRSI.quickMapsList{1}.abs = 1;             % Use magnitude spec
opts.MRSI.quickMapsList{1}.interpolation = 1;   % Spatial interpolation

% You can also add a square mask to accelerate processing
% Top right (confirm) corner has index [min(x) max(y)]
% for example:
% opts.MRSI.outerMask.x = [2 33];
% opts.MRSI.outerMask.y = [4 42];
% opts.MRSI.outerMask.z = [1 5];

% the default is no outer mask
opts.MRSI.outerMask.x = [];
opts.MRSI.outerMask.y = [];
opts.MRSI.outerMask.z = [];

% Or you can interactivley define the mask
opts.MRSI.outerMask.Interactive = 0;


% The data is preprocessed prior to modeling. This includes eddy current
% correction (optional), removal of residual water (HSVD or L2-basis),
% removal of lipid (L2-basis or L2-mask), and frequency correction.

polResidNAA = 1; % Polarity NAA

% Options for auto phasing
% opts.MRSI.phase.type = 'none';
% opts.MRSI.phase.type = 'Cr-Cho';
opts.MRSI.phase.type = 'auto_phase';
opts.MRSI.phase.limits = [1.7,2.2];

% Options for MRSI nuisance signal removal are:
% No water removal
% opts.MRSI.NuisanceRemoval.water.type = 'none';

% HSVD water removal
% opts.MRSI.NuisanceRemoval.water.type = 'HSVD';
% opts.MRSI.NuisanceRemoval.water.comp = 32; % Max number of components removed in HSVD

% L2 removal with simulated water basis functions
opts.MRSI.NuisanceRemoval.water.type = 'L2-basis';
opts.MRSI.NuisanceRemoval.water.basisArguments.Components = 2000;
opts.MRSI.NuisanceRemoval.water.basisArguments.lineWidthRange= [10 80];
opts.MRSI.NuisanceRemoval.water.basisArguments.PPMRange = [4.2 5.1];
opts.MRSI.NuisanceRemoval.water.basisArguments.beta = 10;
opts.MRSI.NuisanceRemoval.water.basisArguments.plotBasis = false;

% No lipid removal
opts.MRSI.NuisanceRemoval.lipid.type = 'none';

% L2 removal with simulated lipid basis functions
% opts.MRSI.NuisanceRemoval.lipid.type = 'L2-basis';
% opts.MRSI.NuisanceRemoval.lipid.basisArguments.Components = 2000;
% opts.MRSI.NuisanceRemoval.lipid.basisArguments.lineWidthRange= [10 80];
% opts.MRSI.NuisanceRemoval.lipid.basisArguments.PPMRange = [0.3 1.900];
% opts.MRSI.NuisanceRemoval.lipid.basisArguments.beta = 100;
% opts.MRSI.NuisanceRemoval.lipid.basisArguments.plotBasis = false;

% L2 removal with lipids defined by lip mask, requires segmentation!
% opts.MRSI.NuisanceRemoval.lipid.type = 'L2-mask';
% opts.MRSI.NuisanceRemoval.lipid.basisArguments.Components = 32;
% opts.MRSI.NuisanceRemoval.lipid.basisArguments.lineWidthRange= [10 80];
% opts.MRSI.NuisanceRemoval.lipid.basisArguments.PPMRange = [0.3 1.85];
% opts.MRSI.NuisanceRemoval.lipid.basisArguments.beta = 0.01;
% opts.MRSI.NuisanceRemoval.lipid.basisArguments.plotBasis = false;

% Options for MRSI frequency alignment are:
% No frequency alignment
% opts.MRSI.FreqAlign.type = 'none';

% Cross-correlation alignment of frequencies defined below
opts.MRSI.FreqAlign.type = 'CC';
opts.MRSI.FreqAlign.frequencies = [2.01,3.03,3.22];
opts.MRSI.FreqAlign.polarity = [1,1,1];
opts.MRSI.FreqAlign.lim = [1.85,4];
opts.MRSI.FreqAlign.realpart = 0;
opts.MRSI.FreqAlign.zerofill = 0;

% Cross-correlation alignment of frequencies defined below after Wavelet filter of baseline if lipid/noise > thresh
% opts.MRSI.FreqAlign.type = 'CCwithLipRemoval';
% opts.MRSI.FreqAlign.thresh = 10;
% opts.MRSI.FreqAlign.frequencies = [2.01,3.03,3.22];
% opts.MRSI.FreqAlign.polarity = [1,1,1];
% opts.MRSI.FreqAlign.lim = [1.85,4];
% opts.MRSI.FreqAlign.realpart = 0;
% opts.MRSI.FreqAlign.zerofill = 0;

% Flags for maximum echo processing
opts.MRSI.MaxEcho.separate = 0;

% This will generate amplitude integral maps for quick inspection. You can
% define different regions and spectra to be used.
opts.MRSI.quickMapsList{2}.specs = {'A'};
opts.MRSI.quickMapsList{2}.target = 'processed';
opts.MRSI.quickMapsList{2}.names.A = {'tNAA','tCr','tCho','tLip','H2O'};
opts.MRSI.quickMapsList{2}.limits.A = [1.95, 2.1; 2.95, 3.12; 3.12, 3.25;0, 1.95;4.1, 6;];
opts.MRSI.quickMapsList{2}.abs = 1;             % Use magnitude spec
opts.MRSI.quickMapsList{2}.interpolation = 1;   % Spatial interpolation


% These are the settings for cosmetically enhanced outputs of the processed
% spectra including line broadening and zero-filling. This will only be
% applied to the exported spectra and not to the spectra that are modelled
% in the next step
opts.MRSI.cosemtics.GaussianLB = 3;
opts.MRSI.cosemtics.ZeroFillFactor = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 3. SPECIFY MRS DATA AND STRUCTURAL IMAGING FILES %%
% When using single-average Siemens RDA or DICOM files, specify their
% folders instead of single files!

% Specify metabolite data
% (MANDATORY)
files       = {which(fullfile('exampledata','mrsi','GE','TE_30_phantom_acr','mrs_nii','mrsi_sag','1sl_metabolite_te_30_mrsi_raw_act_noID.nii.gz'))};

% Specify water reference data for eddy-current correction (same sequence as metabolite data!)
% (OPTIONAL)
% Leave empty for GE P-files (.7) - these include water reference data by
% default.
files_ref   =  {};

% Specify water data for quantification (e.g. short-TE water scan)
% (OPTIONAL)
files_w     = {};

% Specify T1-weighted structural imaging data
% (OPTIONAL)
% Link to single NIfTI (*.nii) files for Siemens and Philips data
% Link to DICOM (*.dcm) folders for GE data
files_nii   = {which(fullfile('exampledata','mrsi','GE','TE_30_phantom_acr','anat_nii','T1w_anat.nii.gz'))};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 4. SPECIFY STAT FILE %%%
% Supply location of a csv file, which contains possible correlation
% measures and group variables. Each column must start with the name of the
% measure. For the grouping variable use 'group' and numbers between 1 and
% the number of included groups. If no group is supplied the data will be
% treated as one group.

file_stat = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5. SPECIFY OUTPUT FOLDER %%
% The Osprey data container will be saved as a *.mat file in the output
% folder that you specify below. In addition, any exported files (for use
% with jMRUI, TARQUIN, or LCModel) will be saved in sub-folders.

% Specify output folder
% (MANDATORY)
data_folder = fileparts(which(fullfile('exampledata','mrsi','GE','TE_30_phantom_acr','jobMRSI_TE_30_in_vitro_P_NII_sag.m')));
outputFolder = fullfile(data_folder, 'derivatives_sag');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
