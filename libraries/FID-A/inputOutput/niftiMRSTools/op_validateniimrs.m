% op_validateniimrs.m
% Georg Oeltzschner, Johns Hopkins University 2024.
%
% USAGE:
% op_validateniimrs(nii_file);
% 
% DESCRIPTION:
% Validates a FID-A struct similar to the 'nifti_mrs/validator' class by
% William T. Clarke (https://github.com/wtclarke/nifti_mrs_tools/blob/master/src/nifti_mrs/validator.py)
%
% Requires a valid 'definitions.json' file from https://github.com/wtclarke/mrs_nifti_standard
%
% See the NIfTI-MRS specification under
% https://docs.google.com/document/d/1tC4ugzGUPLoqHRGrWvOcGCuCh_Dogx_uu0cxKub0EsM/edit
%
% DEPENDENCIES:
% This function requires the dcm2nii toolbox (Xiangrui Li) to be on the
% MATLAB path
% https://github.com/xiangruili/dicm2nii
%
% INPUTS:
% nii_file_or_struct   = filename of NIfTI MRS file (*.nii or *.nii.gz) to load
%                        OR NIfTI struct (loaded with nii_tool from dicm2nii toolbox).
%
% OUTPUTS:
% nii         = Same as input. Not used. 

function op_validateniimrs(nii_file_or_struct)
%function op_validateniimrs(nii_file_or_struct)

%%% Read in the data using the dicm2nii toolbox %%%
% (https://github.com/xiangruili/dicm2nii)

if isstruct(nii_file_or_struct)
    nii = nii_file_or_struct;
else
    try
        nii = nii_tool('load', nii_file_or_struct);
    catch ME
        switch ME.identifier
            case 'MATLAB:UndefinedFunction'
                error(['Cannot find the function ''nii_tool.m''.' ...
                    ' Please ensure that you have downloaded the required', ...
                    ' dcm2nii toolbox (https://github.com/xiangruili/dicm2nii)', ...
                    ' and added it to your MATLAB path.']);
            otherwise
                rethrow(ME);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Validate data inside the NIfTI image %%%

% 1. Check data is complex
if isreal(nii.img)
    error('Data must be complex.');
end

% 2. Check number of dimensions is between 4 and 7
if ndims(nii.img) < 4 || ndims(nii.img) > 7
    error('Data must have between 4 and 7 dimensions. It has %i.\n', ndims(nii.img));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Validate NIfTI header %%%

% 1. Check data type is complex (32 or 1792)
if ~ismember(nii.hdr.datatype, [32 1792])
    error('Data type is not complex.');
end

% 2. Check dwell time
if nii.hdr.pixdim(5) < 0 || nii.hdr.pixdim(5) > 1
    error('Dwell time (%0.5g s) is unrealistic. \n', nii.hdr.pixdim(5));
end

% 3. Check intent name is of format mrs_vMajor_minor
pat = regexpPattern('mrs_v\d+_\d+');
if length(extract(nii.hdr.intent_name, pat)) < 1
    error('Intent string %s does not match ''mrs_vMajor_minor''. \n', nii.hdr.intent_name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Validate NIfTI header extension %%%

% 1. Check that it is json formatted string
try
    hdr_ext = jsondecode(nii.ext.edata_decoded);
catch
    error('Header extension is not json deserialisable.');
end

% 2. Check that it contains the two required bits of metadata
% SpectrometerFrequency
if isfield(hdr_ext, 'SpectrometerFrequency')
    if ~isa(hdr_ext.SpectrometerFrequency, 'double')
        error('SpectrometerFrequency must be a float of double precision (or an array of floats).');
    end
else
    error('Header extension must contain SpectrometerFrequency.');
end

% ResonantNucleus
if isfield(hdr_ext, 'ResonantNucleus')
    if ~isa(hdr_ext.ResonantNucleus, 'cell')
        error('ResonantNucleus must be a cell array of strings.');
    else
        if ~isa(hdr_ext.ResonantNucleus{1}, 'char')
            error('ResonantNucleus must be a cell array of strings.');
        end
    end
else
    error('Header extension must contain ResonantNucleus.');
end

% 3. Check that it contains any required dimension information.
% Calculate implied size
dimension_sizes     = size(nii.img);
data_dimensions     = length(dimension_sizes);
hdr_ext_fields      = fieldnames(hdr_ext);

% Load JSON definitions
if exist('definitions2.json', 'file')
    % Load the proper NIfTI-MRS definition file
    fid         = fopen('definitions.json');
    rawtxt      = fread(fid,inf);
    fclose(fid);
    jsonString  = char(rawtxt');
    jsonStruct  = jsondecode(jsonString);
    % Extract the necessary keys and values
    nifti_mrs_version   = jsonStruct.nifti_mrs_version;
    required            = jsonStruct.required;
    standard_defined    = jsonStruct.standard_defined;
    dimension_tags      = jsonStruct.dimension_tags;

    % Convert dimension tags to key/value map (consistent with the original
    % (pre-JSON) definition in the nifti_mrs_definitions sub-function
    % below:
    dimension_tags      = containers.Map(fieldnames(dimension_tags), struct2cell(dimension_tags));

    % Now we need to map JSON data types to MATLAB data types
    standard_defined    = remapJSONDataTypes(standard_defined);
    % Need to add 'struct' for 'ProcessingApplied' as the MATLAB JSON
    % reader loads an array of JSON object as a struct, not a cell (as we
    % stored it)
    standard_defined.ProcessingApplied;
else
    % Fall back on the (potentially non-compliant) version in this function
    [nifti_mrs_version, dimension_tags, required, standard_defined] = nifti_mrs_definitions;
end

% Check that N-dimensional data has all dim_(up-to-n) tags 
for ddx = 5:8
    
    if data_dimensions > ddx-1

        if ismember(sprintf('dim_%i', ddx), hdr_ext_fields)
            % Check if dimension tags have compliant values
            if ~isKey(dimension_tags, hdr_ext.(sprintf('dim_%i', ddx)))
                error('dim_%i must be a defined tag.\n', ddx);
            end
            % Check if (optional) dimension info tags are strings
            if ismember(sprintf('dim_%i_info', ddx), hdr_ext_fields)
                if ~isa(hdr_ext.(sprintf('dim_%i_info', ddx)), 'char')
                    error('dim_%i_info must be a string.\n', ddx);
                end
            end
            % Check if (optional) dimension info tags are structs (might already be
            % caught by JSON conversion)
            if ismember(sprintf('dim_%i_header', ddx), hdr_ext_fields)
                if ~isa(hdr_ext.(sprintf('dim_%i_header', ddx)), 'struct')
                    error('dim_%i_info must be a struct.\n', ddx);
                end
            end

        else
            % Data with N dimensions needs dim_(5 to n)
            error('With %i dimensions the header extension must contain dim_%i.\n', data_dimensions, ddx)
        end
    else
        % This information shouldn't exist as it refers to data in a dimension higher than that specified
        if ismember(sprintf('dim_%i', ddx), hdr_ext_fields) || ...
           ismember(sprintf('dim_%i_info', ddx), hdr_ext_fields) || ...
           ismember(sprintf('dim_%i_header', ddx), hdr_ext_fields)
            error('dim_%i tags exceed specified dimensions (%i).\n', ddx, data_dimensions);
        end
    end

end

% Additional check that dim{0-4} tags don't exist
for ddx = 0:4
    if ismember(sprintf('dim_%i', ddx), hdr_ext_fields) || ...
       ismember(sprintf('dim_%i_info', ddx), hdr_ext_fields) || ...
       ismember(sprintf('dim_%i_header', ddx), hdr_ext_fields)
        error('dim_%i{_info,_header} tags are forbidden; dim_N... can only take values 5-7.\n', ddx);
    end
end

% 4. Check standard defined data types
standard_fields = fieldnames(standard_defined);
% Check every header extension field against the list of standard defined
% data types
for kk = 1:length(hdr_ext_fields)
    idx = find(ismember(standard_fields, hdr_ext_fields{kk}));
    if ~isempty(idx)
        % If the header extension field is a standard defined one, check if
        % the data types are compliant 
        fprintf('Checking %s ... \n', hdr_ext_fields{kk});
        if ~check_type(hdr_ext.(hdr_ext_fields{kk}), standard_defined.(standard_fields{idx}).type)
            error('%s must be a %s, but is a %s.\n', hdr_ext_fields{kk}, standard_defined.(standard_fields{idx}).type, class(hdr_ext.(hdr_ext_fields{kk})));
        end
    end
end

% 5. Check user-defined format
pat = regexpPattern('^dim_[567](_((info)|(header)))?$');
for kk = 1:length(hdr_ext_fields)
    idx = find(ismember(standard_fields, hdr_ext_fields{kk}));
    % If it's not a standard key, a dim_ key or one of the required ones...
    if isempty(idx) && ...
            ~strcmpi(hdr_ext_fields{kk}, 'SpectrometerFrequency') && ...
            ~strcmpi(hdr_ext_fields{kk}, 'ResonantNucleus') && ...
            length(extract(hdr_ext_fields{kk}, pat)) < 1
        % ... then it must be a user-defined key
        % Check if it is a struct and has a field 'Description' 
        if ~isa(hdr_ext.(hdr_ext_fields{kk}), 'struct')
            error('User-defined key %s must be a JSON object. \n', hdr_ext_fields{kk});
        end
        if ~ismember('Description', fieldnames(hdr_ext.(hdr_ext_fields{kk})))
            error('User-defined key %s must be a JSON object and include a ''Description''. \n', hdr_ext_fields{kk});
        end
    end
end

% 6. Check dynamic header validity
for ddx = 5:8

    if ismember(sprintf('dim_%i_header', ddx), hdr_ext_fields)
        % Allowed formats:
        % - array of the same length as the dimension
        % - json/struct with 'start' and 'increment' fields
        % - if non-standard header, require a nested 'Value' + 'Description'
        % json/struct

        dim_header_fields = fieldnames(hdr_ext.(sprintf('dim_%i_header', ddx)));
        for kk = 1:length(dim_header_fields)
            if ismember(dim_header_fields{kk}, standard_fields)
                % Test standard formats
                test_dyn_header_format(hdr_ext.(sprintf('dim_%i_header', ddx)).(dim_header_fields{kk}), ddx, dimension_sizes);
            else
                % For non-stndard dim_header tags, check they have 'Value'
                % and 'Description' fields
                if ismember('Value', fieldnames(hdr_ext.(sprintf('dim_%i_header', ddx)).(dim_header_fields{kk}))) && ...
                   ismember('Description', fieldnames(hdr_ext.(sprintf('dim_%i_header', ddx)).(dim_header_fields{kk})))
                    test_dyn_header_format(hdr_ext.(sprintf('dim_%i_header', ddx)).(dim_header_fields{kk}).Value, ddx, dimension_sizes);
                else
                    error('dim_%i_header with non-standard tag must contain a ''Value'' and ''Description'' key.\n', ddx);
                end
            end
        end

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Validate spectral width %%%

% If a SpectralWidth field is present, check that it matches the dwell
% time (stored in pixdim(5)).
if ismember('SpectralWidth', hdr_ext_fields)

    if abs(hdr_ext.SpectralWidth - 1/nii.hdr.pixdim(5)) > 1e-2
        error('SpectralWidth (%0.2f Hz) does not match 1/dwelltime (%0.2f Hz). \n', hdr_ext.SpectralWidth, 1/nii.hdr.pixdim(5));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Print success message
fprintf('NIfTI-MRS container successfully validated! \n');

end


%%%% INLINE FUNCTIONS BELOW %%%%

function [nifti_mrs_version, dimension_tags, required, standard_defined] = nifti_mrs_definitions

% Definitions of NIfTI-MRS standard meta data and dimension tags.

% Define nifti-mrs version number here.
% First element is major version, second is minor
nifti_mrs_version = [0, 9];

% Possible dimension tags and descriptions
dimension_tags = containers.Map({'DIM_COIL', 'DIM_DYN', 'DIM_INDIRECT_0', 'DIM_INDIRECT_1', 'DIM_INDIRECT_2', ...
    'DIM_PHASE_CYCLE', 'DIM_EDIT', 'DIM_MEAS', 'DIM_USER_0', 'DIM_USER_1', 'DIM_USER_2', 'DIM_ISIS', 'DIM_METCYCLE'}, ...
    {'For storage of data from each individual receiver coil element.', ...
    'For storage of each individual acquisition transient. E.g. for post-acquisition B0 drift correction.', ...
    'The indirect detection dimension - necessary for 2D (and greater) MRS acquisitions.', ...
    'The indirect detection dimension - necessary for 2D (and greater) MRS acquisitions.', ...
    'The indirect detection dimension - necessary for 2D (and greater) MRS acquisitions.', ...
    'Used for increments of phase-cycling, for example in dephasing unwanted coherence order pathways, or TPPI for 2D spectra.', ...
    'Used for edited MRS techniques such as MEGA or HERMES.', ...
    'Used to indicate multiple repeats of the full sequence contained within the same original data file.', ...
    'User defined dimension.', ...
    'User defined dimension.', ...
    'User defined dimension.', ...
    'Dimension for storing image-selected in vivo spectroscopy (ISIS) acquisitions.', ...
    'For storing metabolite cycled data.'});

% Required metadata fields
required.SpectrometerFrequency = struct('doc', 'Precession frequency in MHz of the nucleus being addressed for each spectral axis.', 'type', {[], 'double'});
required.ResonantNucleus = struct('doc', 'Must be one of the DICOM recognised nuclei “1H”, “3HE”, “7LI”, “13C”, “19F”, “23NA”, “31P”, “129XE” or one named in the specified format. I.e. Mass number followed by the chemical symbol in uppercase.', 'type', {{}, 'char'});

% Defined metadata fields
% 5.1 MRS specific Tags
% 'SpectralWidth'
% 'EchoTime'
% 'RepetitionTime'
% 'InversionTime'
% 'MixingTime'
% 'AcquisitionStartTime'
% 'ExcitationFlipAngle'
% 'TxOffset'
% 'VOI'
% 'WaterSuppressed'
% 'WaterSuppressionType'
% 'SequenceTriggered'
% 5.2 Scanner information
% 'Manufacturer'
% 'ManufacturersModelName'
% 'DeviceSerialNumber'
% 'SoftwareVersions'
% 'InstitutionName'
% 'InstitutionAddress'
% 'TxCoil'
% 'RxCoil'
% 5.3 Sequence information
% 'SequenceName'
% 'ProtocolName'
% 5.4 Sequence information
% 'PatientPosition'
% 'PatientName'
% 'PatientID'
% 'PatientWeight'
% 'PatientDoB'
% 'PatientSex'
% 5.5 Provenance and conversion metadata
% 'ConversionMethod'
% 'ConversionTime'
% 'OriginalFile'
% 5.6 Spatial information
% 'kSpace'
% 5.7 Editing Pulse information structure
% 'EditCondition'
% 'EditPulse'
% 5.8 Processing Provenance
% 'ProcessingApplied'

% These fields are optional but must not be redefined.
% Format is a struct of structs containing type, unit string, doc string, anonymisation state

standard_defined = struct( ...
    'SpectralWidth', struct('type', 'double', 'unit', 'Hz', 'doc', 'The spectral bandwidth of the MR signal that is sampled. Inverse of the dwell time. NIfTI-MRS standard compliant software will always preferentially infer the spectral width from the dwell time stored in the NIfTI pixdim field. Units: hertz', 'anonymisation_state', false), ...
    'EchoTime', struct('type', 'double', 'unit', 's', 'doc', 'Time from centroid of excitation to start of FID or centre of echo. Units: Seconds', 'anonymisation_state', false), ...
    'RepetitionTime', struct('type', 'double', 'unit', 's', 'doc', 'Sequence repetition time. Units: Seconds', 'anonymisation_state', false), ...
    'InversionTime', struct('type', 'double', 'unit', 's', 'doc', 'Inversion time. Units: Seconds', 'anonymisation_state', false), ...
    'MixingTime', struct('type', 'double', 'unit', 's', 'doc', 'Mixing time in e.g. STEAM sequence. Units: Seconds', 'anonymisation_state', false), ...
    'AcquisitionStartTime', struct('type', 'double', 'unit', 's', 'doc', 'Time, relative to EchoTime, that the acquisition starts. Positive values indicate a time after the EchoTime, negative indicate before the EchoTime, a value of zero indicates no offset. Units: Seconds', 'anonymisation_state', false), ...
    'ExcitationFlipAngle', struct('type', 'double', 'unit', 'degrees', 'doc', 'Nominal excitation pulse flip-angle', 'anonymisation_state', false), ...
    'TxOffset', struct('type', 'double', 'unit', 'ppm', 'doc', 'Transmit chemical shift offset from SpectrometerFrequency', 'anonymisation_state', false), ...
    'VOI', struct('type', {{'double', 'double', 'double'}, []}, 'unit', [], 'doc', 'VoI localisation volume for MRSI sequences. Stored as a 4 x 4 affine using identical conventions to the xform NIfTI affine matrix. Not defined for data stored with a single spatial voxel', 'anonymisation_state', false), ...
    'WaterSuppressed', struct('type', 'logical', 'unit', [], 'doc', 'Boolean value indicating whether data was collected with (True) or without (False) water suppression.', 'anonymisation_state', false), ...
    'WaterSuppressionType', struct('type', 'char', 'unit', [], 'doc', 'Type of water suppression used.', 'anonymisation_state', false), ...
    'SequenceTriggered', struct('type', 'logical', 'unit', [], 'doc', 'Boolean value indicating whether the sequence is triggered. If triggered the repetition time might not be constant.', 'anonymisation_state', false), ...
    'Manufacturer', struct('type', 'char', 'unit', [], 'doc', 'Manufacturer of the device. DICOM tag (0008,0070).', 'anonymisation_state', false), ...
    'ManufacturersModelName', struct('type', 'char', 'unit', [], 'doc', "Manufacturer's model name of the device. DICOM tag (0008,1090).", 'anonymisation_state', true), ...
    'DeviceSerialNumber', struct('type', 'char', 'unit', [], 'doc', "Manufacturer's serial number of the device. DICOM tag (0018,1000).", 'anonymisation_state', true), ...
    'SoftwareVersions', struct('type', 'char', 'unit', [], 'doc', "Manufacturer's designation of the software version. DICOM tag (0018,1020).", 'anonymisation_state', false), ...
    'InstitutionName', struct('type', 'char', 'unit', [], 'doc', "Institution's Name. DICOM tag (0008,0080).", 'anonymisation_state', false), ...
    'InstitutionAddress', struct('type', 'char', 'unit', [], 'doc', "Institution's address. DICOM tag (0008,0081).", 'anonymisation_state', false), ...
    'TxCoil', struct('type', 'char', 'unit', [], 'doc', "Name or description of transmit RF coil.", 'anonymisation_state', false), ...
    'RxCoil', struct('type', 'char', 'unit', [], 'doc', "Name or description of receive RF coil.", 'anonymisation_state', false), ...
    'SequenceName', struct('type', 'char', 'unit', [], 'doc', "User defined name. DICOM tag (0018,0024).", 'anonymisation_state', false), ...
    'ProtocolName', struct('type', 'char', 'unit', [], 'doc', "User-defined description of the conditions under which the Series was performed. DICOM tag (0018,1030).", 'anonymisation_state', false), ...
    'PatientPosition', struct('type', 'char', 'unit', [], 'doc', "Patient position descriptor relative to the equipment. DICOM tag (0018,5100). Must be one of the DICOM defined code strings e.g. HFS, HFP.", 'anonymisation_state', false), ...
    'PatientName', struct('type', 'char', 'unit', [], 'doc', "Patient's full name. DICOM tag (0010,0010).", 'anonymisation_state', true), ...
    'PatientID', struct('type', 'char', 'unit', [], 'doc', "Patient identifier. DICOM tag (0010,0020).", 'anonymisation_state', true), ...
    'PatientWeight', struct('type', 'double', 'unit', 'kg', 'doc', "Weight of the Patient in kilograms. DICOM tag (0010,1030).", 'anonymisation_state', false), ...
    'PatientDoB', struct('type', 'char', 'unit', [], 'doc', "Date of birth of the named Patient. YYYYMMDD. DICOM tag (0010,0030).", 'anonymisation_state', true), ...
    'PatientSex', struct('type', 'char', 'unit', [], 'doc', "Sex of the named Patient. 'M', 'F', 'O'. DICOM tag (0010,0040).", 'anonymisation_state', false), ...
    'ConversionMethod', struct('type', 'char', 'unit', [], 'doc', "Description of the process or program used for conversion. May include additional information like software version.", 'anonymisation_state', false), ...
    'ConversionTime', struct('type', 'char', 'unit', [], 'doc', "Time and date of conversion. ISO 8601 compliant format.", 'anonymisation_state', false), ...
    'OriginalFile', struct('type', {{'cell','char'}}, 'unit', [], 'doc', "Name and extension of the original file(s).", 'anonymisation_state', true), ...
    'kSpace', struct('type', 'logical', 'unit', [], 'doc', "Three element list, corresponding to the first three spatial dimensions. If True the data is stored as a dense k-space representation.", 'anonymisation_state', false), ...
    'EditCondition', struct('type', 'char', 'unit', [], 'doc', "List of strings that index the entries of the EditPulse structure that are used in this data acquisition. Typically used in dynamic headers (dim_N_header).", 'anonymisation_state', false), ...
    'EditPulse', struct('type', 'char', 'unit', [], 'doc', "Structure defining editing pulse parameters for each condition. Each condition must be assigned a key.", 'anonymisation_state', false), ...
    'ProcessingApplied', struct('type', 'struct', 'unit', [], 'doc', "Describes and records the processing steps applied to the data.", 'anonymisation_state', false));
    
end

function trueOrFalse = check_type(value, json_type)
% Checks that `values` is of type json_type

if ischar(json_type)
    if isa(value, json_type)
        trueOrFalse = true;
    else
        trueOrFalse = false;
    end
else
    for kk = 1:length(json_type)
        if isa(value, json_type{kk})
            trueOrFalse = true;
            return;
        end
    end
    % If this is reached, 'value' did not match any of the 'json_type{kk}'
    error('%s is not a %s, but is a %s. \n', value, json_type{kk}, class(value));
end

end

function test_dyn_header_format(x, ddx, dimension_sizes)
% Checks that the dynamic header have correct content and formatting

% Make sure it's either an array of values or a json/struct
if ~isa(x, 'struct') && ~isvector(x)
    error('dim_%i_header not an array or struct/json. \n', ddx);
end

% Make sure that it has 'start' and 'increment' fields if json/struct
if isa(x, 'struct')
    if ~ismember('start', fieldnames(x)) || ~ismember('increment', fieldnames(x))
        error('dim_%i_header is a struct/json but does not contain ''start'' or ''increment''. \n', ddx);
    end
end

% Make sure that it has same length of its assigned dimension if array
if ~isa(x, 'struct') && isvector(x)
    dim_size = dimension_sizes(ddx);
    if length(x) ~= dim_size
        error('dim_%i_header is an array but its size %i does not match the dimension size %i. \n', ddx, length(x), dim_size);
    end
end

end

function standard_defined    = remapJSONDataTypes(standard_defined)
% Changes strings describing JSON datatype to matching strings describing
% MATLAB datatypes

standard_fields = fieldnames(standard_defined);
% Loop through all fields
for rr = 1:length(standard_fields)
    
    % How many types are allowed?
    standard_defined.(standard_fields{rr}).type = standard_defined.(standard_fields{rr}).type'; % for some weird reason, need to transpose here (??)
    nTypes = length(standard_defined.(standard_fields{rr}).type);
    
    for nn = 1:nTypes
        nthType = standard_defined.(standard_fields{rr}).type{nn};

        switch nthType
            case 'number'
                standard_defined.(standard_fields{rr}).type{nn} = 'double';
            case 'string'
                standard_defined.(standard_fields{rr}).type{nn} = 'char';
            case 'array'
                standard_defined.(standard_fields{rr}).type{nn} = 'cell';
            case 'bool'
                standard_defined.(standard_fields{rr}).type{nn} = 'logical';
            case 'object'
                standard_defined.(standard_fields{rr}).type{nn} = 'struct';
        end

    end

end

end

