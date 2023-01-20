function [ModelParameter] = Osprey_gLCM(DataToModel, JsonModelFile, CheckGradient)
%% Global function for new Osprey LCM
% Inputs:   DataToModel - FID-A/Osprey struct with data or cell of structs
%           JsonModelFile - Master model file for all steps
%           CheckGradient - Do a gradient check in the lsqnonlin solver    
% Outputs:  explicit - Struct with model parameters
%           implicit - NII results file 

%% 0. Check inputs
arguments
    % Argument validation introduced in relatively recent MATLAB versions
    % (2019f?)
    DataToModel {isStructOrCell}
    JsonModelFile string
    CheckGradient double {mustBeNumeric} = 0; % optional
end

%% What happens here:
%   1. Decode model json file 
%   2. Prepare basisset matrix (and export as NII)
%       2a. Do on-the-fly generation of MMs
%   3. Prepare data according to model json & Run steps defined in model json
%   4. Save results (and export as NII)

%% 1. Decode model json file
% Read the json file and generate a ModelProcedure struct from it which will guide the
% rest of the analysis. Catch missing parameters here?
ModelProcedure = jsonToStruct(JsonModelFile);

%% 2. Prepare basisset matrix (and export as NII)
% Load basisset files, add MMs if needed, resample basis sets according to
% the DataToModel. Generate a basisset matrix for each step? including the
% indirect dimensions for MSM.
if length(ModelProcedure.basisset.file) == 1  
    basisSet = load(ModelProcedure.basisset.file{1});
    basisSet = basisSet.BASIS;
    basisSet = recalculateBasisSpecs(basisSet);
    basisSet = fit_sortBasisSet(basisSet);
else
    for bb = 1 : length(ModelProcedure.basisset.file)
        if bb == 1
            basisSet = load(ModelProcedure.basisset.file{bb});
            basisSet = basisSet.BASIS;
            basisSet = recalculateBasisSpecs(basisSet);
            basisSet = fit_sortBasisSet(basisSet);
        else
            basisSetToAdd = load(ModelProcedure.basisset.file{bb});
            basisSetToAdd = basisSetToAdd.BASIS;
            basisSetToAdd = recalculateBasisSpecs(basisSetToAdd);
            basisSetToAdd = fit_sortBasisSet(basisSetToAdd);
            basisSet.fids = cat(3,basisSet.fids,basisSetToAdd.fids);
            basisSet.specs = cat(3,basisSet.specs,basisSetToAdd.specs);
        end
    end
    basisSet.sz = size(basisSet.fids);
    basisSet.nExtra = basisSet.sz(3);
end

%% 3. Prepare data according to model json and model data
% Prepare data for each fit step, again, including the indirect dimensions
% for MSM
% Create the spline basis functions for the given resolution, fit range,
% and knot spacing parameter.
for ss = 1 : length(ModelProcedure.Steps)
    opts.ModelFunction      = ModelProcedure.Steps(ss).ModelFunction; % get model function
    opts.baseline             = ModelProcedure.Steps(ss).fit_opts.baseline; % setup baseline options
    opts.optimDomain        = ModelProcedure.Steps(ss).fit_opts.optimDomain; % do the least-squares optimization in the frequency domain
    opts.optimSignalPart    = ModelProcedure.Steps(ss).fit_opts.optimSignalPart; % do the least-squares optimization over the real part of the spectrum
    opts.optimFreqFitRange  = ModelProcedure.Steps(ss).fit_opts.ppm; % set the frequency-domain fit range to the fit range specified in the Osprey container
    opts.solver  = ModelProcedure.Steps(ss).fit_opts.solver; % set the solver specified in the Osprey container
    opts.CheckGradient = CheckGradient; % Do a gradient check in the lsqnonlin solver
    if isfield(ModelProcedure.Steps(ss),'parameter')
        opts.parameter  = ModelProcedure.Steps(ss).parameter; % specify parmetrizations constructor
    end
    if isfield(ModelProcedure.Steps(ss),'parametrizations')
        opts.parametrizations  = ModelProcedure.Steps(ss).parametrizations; % specify parmetrizations constructor
    end
    if ModelProcedure.Steps(ss).extra.flag == 1 && isfield(ModelProcedure.Steps(ss).extra,'DynamicModelJson')
        opts.paraIndirect  = jsonToStruct(ModelProcedure.Steps(ss).extra.DynamicModelJson); % specify parmetrizations constructor for indirect dimension
    end
    if ~iscell(DataToModel)
        % Create an instance of the class
        if ss == 1
            ModelParameter{1} = FitObject(DataToModel, basisSet, opts);
        else
            ModelParameter{1}.updateOptsAccordingToStep(opts);
        end
        ModelParameter{1}.excludeBasisFunctionFromFit('all');
        ModelParameter{1}.includeBasisFunctionInFit(ModelProcedure.Steps(ss).basisset.include);
        % Run steps defined ModelProcedure struct
        % Loop across all steps with anonymous calls defined in the ModelProcedure struct
        ModelParameter{1}.createModel;
    else
        for kk = 1 : length(DataToModel)
            % Create an instance of the class
            if ss == 1
                ModelParameter{kk} = FitObject(DataToModel{kk}, basisSet, opts);
            else
                ModelParameter{kk}.updateOptsAccordingToStep(opts);
            end
            ModelParameter{kk}.excludeBasisFunctionFromFit('all');
            ModelParameter{kk}.includeBasisFunctionInFit(ModelProcedure.Steps(ss).basisset.include);
            % Run steps defined ModelProcedure struct
            % Loop across all steps with anonymous calls defined in the ModelProcedure struct
            ModelParameter{kk}.createModel;
        end
    end
end
    %% 4. Save parameter results (and export as NII)
    % Save ModelParameter to be includeded in the container and final
    % results in NII LCM format for easy reading and visualization
end
function isStructOrCell(f)
    assert(iscell(f) || isstruct(f));
end