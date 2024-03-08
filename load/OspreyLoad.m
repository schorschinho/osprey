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
if ~isempty(MRSCont.files_mm_ref)
    MRSCont.flags.hasMMRef = 1;
end
if ~isempty(MRSCont.files_ref)
    MRSCont.flags.hasRef = 1;
end
if ~isempty(MRSCont.files_w)
    MRSCont.flags.hasWater = 1;
end

% Version check and updating log file
outputFolder = MRSCont.outputFolder;
diary(fullfile(outputFolder, 'LogFile.txt'));
[~,MRSCont.ver.CheckOsp ] = osp_CheckRunPreviousModule(MRSCont, 'OspreyLoad');

% Determine data types
[MRSCont, retMsg,reordered] = osp_detDataType(MRSCont);
MRSCont.flags.reordered = reordered;
% Parse ECC flag entry
if ~MRSCont.flags.isSERIES
    if length(MRSCont.opts.ECC.raw) == 1
        MRSCont.opts.ECC.raw = ones(size(MRSCont.files)) * MRSCont.opts.ECC.raw;
    end
    if length(MRSCont.opts.ECC.mm) == 1 && ~isempty(MRSCont.files_mm)
        MRSCont.opts.ECC.mm = ones(size(MRSCont.files_mm)) * MRSCont.opts.ECC.mm;
    end
end



tempDatasets = [];
MRSCont.opts.MultipleSpectra.metab = fliplr(size(MRSCont.files));
if length(MRSCont.opts.MultipleSpectra.metab) > 1
    tempDatasets(end+1) = MRSCont.opts.MultipleSpectra.metab(2);
else
    MRSCont.opts.MultipleSpectra.metab(2) = 1;
end
MRSCont.opts.MultipleSpectra.mm = [];
MRSCont.opts.MultipleSpectra.mm_ref = [];
MRSCont.opts.MultipleSpectra.ref = [];
MRSCont.opts.MultipleSpectra.w = [];
if MRSCont.flags.hasMM
    MRSCont.opts.MultipleSpectra.mm = fliplr(size(MRSCont.files_mm));
    if length(MRSCont.opts.MultipleSpectra.metab) > 1
        tempDatasets(end+1) = MRSCont.opts.MultipleSpectra.mm(2);
    else
        MRSCont.opts.MultipleSpectra.mm(2) = 1;
    end
end
if MRSCont.flags.hasMMRef
    MRSCont.opts.MultipleSpectra.mm_ref = fliplr(size(MRSCont.files_mm_ref));
    if length(MRSCont.opts.MultipleSpectra.mm_ref) > 1
        tempDatasets(end+1) = MRSCont.opts.MultipleSpectra.mm_ref(2);
    else
        MRSCont.opts.MultipleSpectra.mm_ref(2) = 1;
    end
end
if MRSCont.flags.hasRef
    MRSCont.opts.MultipleSpectra.ref = fliplr(size(MRSCont.files_ref));
    if length(MRSCont.opts.MultipleSpectra.ref) > 1
        tempDatasets(end+1) = MRSCont.opts.MultipleSpectra.ref(2);
    else
        MRSCont.opts.MultipleSpectra.ref(2) = 1;
    end
end
if MRSCont.flags.hasWater
    MRSCont.opts.MultipleSpectra.w = fliplr(size(MRSCont.files_w));
    if length(MRSCont.opts.MultipleSpectra.w) > 1
        tempDatasets(end+1) = MRSCont.opts.MultipleSpectra.w(2);
    else
        MRSCont.opts.MultipleSpectra.w(2) = 1;
    end
end

[maxDatasets,~] = max(tempDatasets);

% Determine number of datasets
MRSCont.nDatasets = size(MRSCont.files,2);
if maxDatasets > 1
    MRSCont.nDatasets(2) = maxDatasets;
    if ~isfield(MRSCont.opts, 'extras')
        MRSCont.opts.extras.names = {};
        MRSCont.opts.extras.exp_var = [];
        for ex = 1 : MRSCont.nDatasets(2)
            MRSCont.opts.extras.names{end+1} = ['Exp_' num2str(ex)];
            MRSCont.opts.extras.exp_var(end+1) = 1;
        end
    end
else
    MRSCont.nDatasets(2) = 1;
end
if MRSCont.opts.MultipleSpectra.metab(2) < maxDatasets
    MRSCont.opts.MultipleSpectra.metab = ones(1,maxDatasets);
else if maxDatasets > 1
    MRSCont.opts.MultipleSpectra.metab = [1:MRSCont.opts.MultipleSpectra.metab(2)];
    else
    MRSCont.opts.MultipleSpectra.metab = [1:MRSCont.nDatasets(1)];
end
    if ~MRSCont.flags.isSERIES
        if length(MRSCont.opts.ECC.raw) == 1
            MRSCont.opts.ECC.raw = ones(size(MRSCont.files)) * MRSCont.opts.ECC.raw;
        else
            MRSCont.opts.ECC.raw = repmat(MRSCont.opts.ECC.raw, [ MRSCont.nDatasets(2) 1]);
        end  
    end
end
if MRSCont.flags.hasMM
    if MRSCont.opts.MultipleSpectra.mm(2) < maxDatasets
        MRSCont.opts.MultipleSpectra.mm = ones(1,maxDatasets);
    else if maxDatasets > 1
        MRSCont.opts.MultipleSpectra.mm = [1:MRSCont.opts.MultipleSpectra.mm(2)];
        else
        MRSCont.opts.MultipleSpectra.mm = [1:MRSCont.nDatasets(1)];
        end
    end
end
if MRSCont.flags.hasMMRef
    if MRSCont.opts.MultipleSpectra.mm_ref(2) < maxDatasets
        MRSCont.opts.MultipleSpectra.mm_ref = ones(1,maxDatasets);
    else if maxDatasets > 1
        MRSCont.opts.MultipleSpectra.mm_ref = [1:MRSCont.opts.MultipleSpectra.mm_ref(2)];
        else
        MRSCont.opts.MultipleSpectra.mm_ref = [1:MRSCont.nDatasets(1)];
    end
    if ~MRSCont.flags.isSERIES
        if length(MRSCont.opts.ECC.mm) == 1 && ~isempty(MRSCont.files_mm)
            MRSCont.opts.ECC.mm = ones(size(MRSCont.mm)) * MRSCont.opts.ECC.mm;
        else
            MRSCont.opts.ECC.mm = repmat(MRSCont.opts.ECC.mm, [ MRSCont.nDatasets(2) 1]);
        end  
    end
    end
end
if MRSCont.flags.hasRef
    if MRSCont.opts.MultipleSpectra.ref(2) < maxDatasets
        MRSCont.opts.MultipleSpectra.ref = ones(1,maxDatasets);
    else if maxDatasets > 1
        MRSCont.opts.MultipleSpectra.ref = [1:MRSCont.opts.MultipleSpectra.ref(2)];
        else
        MRSCont.opts.MultipleSpectra.ref = [1:MRSCont.nDatasets(1)];
        end
    end
end
if MRSCont.flags.hasWater
    if MRSCont.opts.MultipleSpectra.w(2) < maxDatasets
        MRSCont.opts.MultipleSpectra.w = ones(1,maxDatasets);
    else if maxDatasets > 1
        MRSCont.opts.MultipleSpectra.w = [1:MRSCont.opts.MultipleSpectra.w(2)];
        else
        MRSCont.opts.MultipleSpectra.w = [1:MRSCont.nDatasets(1)];
        end
    end
end


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
                if MRSCont.flags.hasWater
                    MRSCont.opts.MultipleSpectra.w = fliplr(size(MRSCont.files_w));
                    if length(MRSCont.opts.MultipleSpectra.w) > 1
                        tempDatasets(end+1) = MRSCont.opts.MultipleSpectra.w(2);
                    else
                        MRSCont.opts.MultipleSpectra.w(2) = 1;
                    end
                    if MRSCont.opts.MultipleSpectra.w(2) < maxDatasets
                        MRSCont.opts.MultipleSpectra.w = ones(1,maxDatasets);
                    else if maxDatasets > 1
                        MRSCont.opts.MultipleSpectra.w = [1:MRSCont.opts.MultipleSpectra.w(2)];
                        else
                        MRSCont.opts.MultipleSpectra.w = [1:MRSCont.nDatasets(1)];
                        end
                    end
                end
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

if MRSCont.flags.isUnEdited
    for kk = 1:MRSCont.nDatasets(1)
        raw                         = MRSCont.raw{kk};                                          % Get the kk-th dataset

        %%% MERGE MULTIPLE DIMENSIONS %%%
        % If the dimensionality of the dataset isn't just along the
        % 'averages' dimension, we re-sort the other dimensions into the
        % 'averages' dimension here. This can be the case when some
        % sequences divide the entire acquisition into dynamics, which are
        % then again divided into transients.
        % This does NOT apply to SPECIAL-localized data - here, the
        % different dimensions have an explicit meaning and need to be
        % preserved.
        % This is also NOT the case for data acquired with the series
        % function!
        if ~MRSCont.flags.isSPECIAL && ~MRSCont.flags.isSERIES
            if raw.dims.extras ~= 0
                % Generate empty struct
                temp = struct;
                % Extract extras and add to the temporary struct
                for pp = 1:raw.sz(raw.dims.extras)
                    extrasToAdd = op_takeextras(raw, pp);
                    temp = op_concatAverages(temp, extrasToAdd);
                end
                % Save back to MRSCont
                raw = temp;
                MRSCont.raw{kk} = raw;
            elseif raw.dims.subSpecs ~= 0
                % Generate empty struct
                temp = struct;
                % Extract subspecs and add to the temporary struct
                for pp = 1:raw.sz(raw.dims.subSpecs)
                    subspecsToAdd = op_takesubspec(raw, pp);
                    temp = op_concatAverages(temp, subspecsToAdd);
                end
                % Save back to MRSCont
                raw = temp;
                MRSCont.raw{kk} = raw;
            end
        end
    end
end

%% Seperate scans for Series scans
if MRSCont.flags.isSERIES
    MRSCont.nDatasets(2) = MRSCont.raw{kk}.sz(MRSCont.raw{kk}.dims.extras);
    tempRaw = MRSCont.raw;
    if MRSCont.flags.hasMM
        tempRawMM = MRSCont.raw_mm;
    end
    if MRSCont.flags.hasMMRef
        tempRawMMref = MRSCont.raw_mm_ref;
    end
    if MRSCont.flags.hasRef
        tempRawRef = MRSCont.raw_ref;
    end
    if MRSCont.flags.hasWater
        tempRawW = MRSCont.raw_w;
    end
    for kk = 1:MRSCont.nDatasets(1)
        for ll = 1 : MRSCont.nDatasets(2)
            MRSCont.opts.MultipleSpectra.metab(ll) = ll;
            extras = op_takeextra(tempRaw{kk}, ll);
             MRSCont.raw{ll,kk} = extras;
        end    
        if MRSCont.flags.hasMM
            for ll = 1 : MRSCont.nDatasets(2)
                MRSCont.opts.MultipleSpectra.mm(ll) = ll;
                extras = op_takeextra(tempRawMM{kk}, ll);
                 MRSCont.raw_mm{ll,kk} = extras;
            end
        end
        if MRSCont.flags.hasMMRef
            for ll = 1 : MRSCont.nDatasets(2)
                MRSCont.opts.MultipleSpectra.mm_ref(ll) = ll;
                extras = op_takeextra(tempRawMMref{kk}, ll);
                 MRSCont.raw_mm_ref{ll,kk} = extras;
            end
        end
        if MRSCont.flags.hasRef
            for ll = 1 : MRSCont.nDatasets(2)
                MRSCont.opts.MultipleSpectra.ref(ll) = ll;
                extras = op_takeextra(tempRawRef{kk}, ll);
                 MRSCont.raw_ref{ll,kk} = extras;
            end
        end
        if MRSCont.flags.hasWater
            for ll = 1 : MRSCont.nDatasets(2)
                MRSCont.opts.MultipleSpectra.w(ll) = ll;
                extras = op_takeextra(tempRawW{kk}, ll);
                 MRSCont.raw_w{ll,kk} = extras;
            end
        end
    end 

    % Parse ECC flag entry
    if length(MRSCont.opts.ECC.raw) == 1
        MRSCont.opts.ECC.raw = ones(size(MRSCont.raw)) * MRSCont.opts.ECC.raw;
    end
    if length(MRSCont.opts.ECC.mm) == 1 && ~isempty(MRSCont.files_mm)
        MRSCont.opts.ECC.mm = ones(size(MRSCont.raw_mm)) * MRSCont.opts.ECC.mm;
    end
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
MRSCont.flags.didLoad           = 1;
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
