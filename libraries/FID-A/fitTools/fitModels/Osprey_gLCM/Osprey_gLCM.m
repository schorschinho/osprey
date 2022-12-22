function [ModelParameter] = Osprey_gLCM(DataToModel,JsonModelFile)
%% Global function for new Osprey LCM
% Inputs:   DataToModel - FID-A/Osprey struct with data or cell of structs
%           JsonModelFile - Master model file for all steps
% Outputs:  explicit - Struct with model parameters
%           implicit - NII results file 


%% What happens here:
%   1. Decode model json file 
%   2. Prepare basisset matrix (and export as NII)
%       2a. Do on-the-fly generation of MMs
%   3. Prepare data according to model json & Run steps defined in model json
%   4. Save results (and export as NII)

%% 1. Decode model json file
% Read the json file and generate a ModelProcedure struct from it which will guide the
% rest of the analysis. Catch missing parameters here?

fid = fopen(JsonModelFile); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
if strcmp('win',osp_platform('filesys'))
    str = strrep(str,'\','\\');
end
str = replace(str, whitespacePattern + '"', '"');
ModelProcedure  = jsondecode(str);

%% 2. Prepare basisset matrix (and export as NII)
% Load basisset files, add MMs if needed, resample basis sets according to
% the DataToModel. Generate a basisset matrix for each step? including the
% indirect dimensions for MSM.

basisSet = load(ModelProcedure.Steps{1}.basisset.file);
basisSet = basisSet.BASIS;
basisSet = osp_recalculate_basis_specs(basisSet);          % HZ re-calculate specs
basisSet = fit_sortBasisSet(basisSet);
%% 3. Prepare data according to model json and model data
% Prepare data for each fit step, again, including the indirect dimensions
% for MSM
% Create the spline basis functions for the given resolution, fit range,
% and knot spacing parameter.
for ss = 1 : length(ModelProcedure.Steps)
    opts.dkntmn             = ModelProcedure.Steps{ss}.fit_opts.baseline.dkntmn; % minimum spacing between two neighboring spline knots
    opts.optimDomain        = ModelProcedure.Steps{ss}.fit_opts.optimDomain; % do the least-squares optimization in the frequency domain
    opts.optimSignalPart    = ModelProcedure.Steps{ss}.fit_opts.optimSignalPart; % do the least-squares optimization over the real part of the spectrum
    opts.optimFreqFitRange  = ModelProcedure.Steps{ss}.fit_opts.ppm; % set the frequency-domain fit range to the fit range specified in the Osprey container
    if isfield(ModelProcedure.Steps{ss},'initials')
        opts.initials  = ModelProcedure.Steps{ss}.initials; % specify initials constructor
    end
    if ~iscell(DataToModel)
        % Create an instance of the class
        if ss == 1
            ModelParameter{1} = FitObject(DataToModel, basisSet, opts);
        else
            ModelParameter{1}.updateOptsAccordingToStep(opts);
        end
        ModelParameter{1}.excludeBasisFunctionFromFit('all');
        ModelParameter{1}.includeBasisFunctionInFit(ModelProcedure.Steps{ss}.basisset.include);
        % Run steps defined ModelProcedure struct
        % Loop across all steps with anonymous calls defined in the ModelProcedure struct
%         ModelParameter{1}.initFit;
        ModelParameter{1}.createModel;
    else
        for kk = 1 : length(DataToModel)
            % Create an instance of the class
            ModelParameter{kk} = FitObject(DataToModel{kk}, basisSet, opts);
            ModelParameter{kk}.excludeBasisFunctionFromFit('all');
            ModelParameter{kk}.includeBasisFunctionInFit(ModelProcedure.Steps{ss}.basisset.include);
            % Run steps defined ModelProcedure struct
            % Loop across all steps with anonymous calls defined in the ModelProcedure struct
%             ModelParameter{kk}.initFit;
            ModelParameter{kk}.createModel;
        end
    end
end
    %% 4. Save parameter results (and export as NII)
    % Save ModelParameter to be includeded in the container and final
    % results in NII LCM format for easy reading and visualization
end