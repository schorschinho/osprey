function [MRSCont] = OspreyLoad(MRSCont)
%% [MRSCont] = OspreyLoad(MRSCont)
%   This function loads the raw MRS data from all major vendors.
%   Data is read from the provided input filenames. It is shaped according 
%   to the type of sequence (un-edited data, MEGA-edited (ON/OFF), 
%   HERMES/HERCULES (A/B/C/D), etc.).
%
%   USAGE:
%       MRSCont = OspreyLoad(MRSCont);
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

% Set flags
if ~isempty(MRSCont.files)
    MRSCont.flags.hasFiles = 1;
end
if ~isempty(MRSCont.files_mm)       %re_mm adding functionality to load MM data
    MRSCont.flags.hasMM = 1;        %re_mm 
end                                 %re_mm
if ~isempty(MRSCont.files_ref)
    MRSCont.flags.hasRef = 1;
end
if ~isempty(MRSCont.files_w)
    MRSCont.flags.hasWater = 1;
end

% Version check and updating log file
outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));
[~,MRSCont.ver.CheckOsp ] = osp_Toolbox_Check ('OspreyLoad',MRSCont.flags.isGUI);

% Determine data types
[MRSCont, retMsg] = osp_detDataType(MRSCont);

% Determine number of datasets
MRSCont.nDatasets = length(MRSCont.files);

% Load raw data (call loaders depending on file type)
switch MRSCont.vendor
    case 'Siemens'
        switch MRSCont.datatype
            case 'TWIX'
                [MRSCont] = osp_LoadTwix(MRSCont);
            case 'RDA'
                [MRSCont] = osp_LoadRDA(MRSCont);
            case 'DICOM'
                [MRSCont] = osp_LoadDICOM(MRSCont);
            otherwise
                msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                fprintf(msg);
                error(msg);
        end
    case 'Philips'
        switch MRSCont.datatype
            case 'SDAT'
                [MRSCont] = osp_LoadSDAT(MRSCont);
            case 'DATA'
                if ~MRSCont.flags.isMRSI
                    [MRSCont] = osp_LoadDATA(MRSCont);
                else
                    [MRSCont] = load_mrsi_data(MRSCont);
                end
            case 'LAB'
                error('Support for Philips RAW/LAB/SIN files coming soon!');
                %[MRSCont] = osp_LoadRAW(MRSCont);
            otherwise
                msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                fprintf(msg);
                error(msg);
        end
    case 'GE'
        switch MRSCont.datatype
            case 'P'
                [MRSCont] = osp_LoadP(MRSCont);
            otherwise
                msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                fprintf(msg);
                error(msg);
        end
    case 'Bruker'
        switch MRSCont.datatype
            case 'fid'
                [MRSCont] = osp_LoadBrukerFid(MRSCont);
            otherwise
                msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                fprintf(msg);
                error(msg);
        end
        
    case 'LCModel'
        switch MRSCont.datatype
            case 'RAW'
                [MRSCont] = osp_LoadRAW(MRSCont);
            otherwise
                msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                fprintf(msg);
                error(msg);
        end
    
    case ''
        % We left the vendor field empty for NIfTI-MRS data
        switch MRSCont.datatype
            case 'NIfTI-MRS'
                [MRSCont] = osp_LoadNII(MRSCont);
            otherwise
                msg = 'Data type not supported. Please contact the Osprey team (gabamrs@gmail.com).';
                fprintf(msg);
                error(msg);
        end
        
    otherwise
        msg = 'Vendor not supported. Please contact the Osprey team (gabamrs@gmail.com).';
        fprintf(msg);
        error(msg);
end

% Perform coil combination (SENSE-based reconstruction if PRIAM flag set)
if ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    if sum(strcmp(MRSCont.datatype, {'DATA', 'LAB', 'P'})) == 1 || ~MRSCont.flags.coilsCombined
        [MRSCont] = osp_combineCoils(MRSCont);
    else
        if ~strcmp(MRSCont.datatype, 'TWIX')
            fprintf('Data type %s %s is already coil-combined.\n', MRSCont.vendor, MRSCont.datatype);
        end
    end
elseif MRSCont.flags.isPRIAM
    [MRSCont] = osp_senseRecon(MRSCont);
elseif MRSCont.flags.isMRSI && ~strcmp(MRSCont.datatype,'DATA')
    [MRSCont] = osp_MRSIRecon(MRSCont);
end

%% If DualVoxel or MRSI we want to extract y-axis scaling
% Creates y-axis range to align the process plots between datasets
if MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI
    MRSCont.plot.load.match = 1; % Scaling between datasets is turned off by default
else
    MRSCont.plot.load.match = 0; % Scaling between datasets is turned off by default
end
MRSCont = osp_scale_yaxis(MRSCont,'OspreyLoad');
%% Clean up and save
% Set exit flags and version
MRSCont.flags.didLoadData           = 1;
diary off

% Save the output structure to the output folder
% Determine output folder
outputFile      = MRSCont.outputFile;
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end


% Optional: Create all pdf figures
if MRSCont.opts.savePDF
    osp_plotAllPDF(MRSCont, 'OspreyLoad')
end

% Gather some more information from the processed data;
if  MRSCont.flags.isGUI
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
else
   save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
end

end