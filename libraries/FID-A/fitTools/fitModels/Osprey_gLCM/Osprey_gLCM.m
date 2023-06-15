function [ModelParameter] = Osprey_gLCM(DataToModel, JsonModelFile,average, NumericJacobian, CheckGradient, BasisSetStruct)
%% Global function for new Osprey LCM
% Inputs:   DataToModel - FID-A/Osprey struct with data or cell of structs
%           JsonModelFile - Master model file for all steps
%           average       - average data prior to modeling (NIfTI-MRS only)
%           NumericJacobian - use numerical jacobian flag
%           CheckGradient - Do a gradient check in the lsqnonlin solver  
%           BasisSetStruct - Include predefined basisset struct
% Outputs:  explicit - Struct with model parameters
%           implicit - NII results file 

%% 0. Check inputs
arguments
    % Argument validation introduced in relatively recent MATLAB versions
    % (2019f?)
    DataToModel {isStructOrCell}
    JsonModelFile string
    average double {mustBeNumeric} = 0; % optional
    NumericJacobian double {mustBeNumeric} = 0; % optional
    CheckGradient double {mustBeNumeric} = 0; % optional
    BasisSetStruct struct = []; %optional
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
if isstruct(ModelProcedure.Steps)
    ModelProcedureCell = cell(size(ModelProcedure.Steps));
    for ss = 1 : size(ModelProcedure.Steps,1)
        ModelProcedureCell{ss} = ModelProcedure.Steps(ss,:);
    end
    ModelProcedure.Steps = ModelProcedureCell;
end
% Check if the R part is used for optimization and a factor 2 zerofilling
zf = 0;
for ss = 1 : length(ModelProcedure.Steps)
    if isfield(ModelProcedure.Steps{ss}.fit_opts,'optimSignalPart')
        if strcmp(ModelProcedure.Steps{ss}.fit_opts.optimSignalPart,'R')
            zf = 0;
        end
    end
end
%% 2. Prepare basisset matrix (and export as NII)
% Load basisset files, add MMs if needed, resample basis sets according to
% the DataToModel. Generate a basisset matrix for each step? including the
% indirect dimensions for MSM.
if isempty(BasisSetStruct)                                                  % User supplied a recalculated basis set
    if length(ModelProcedure.basisset.file) == 1  
        basisSet = load(ConvertRelativePath(ModelProcedure.basisset.file{1})); % Load basis set
        basisSet = basisSet.BASIS;
        basisSet = recalculateBasisSpecs(basisSet);                         % Add ppm axis and frequency domain data
        basisSet = fit_sortBasisSet(basisSet);                              % Sort according to Osprey standard
    else
        for bb = 1 : length(ModelProcedure.basisset.file)           
            if bb == 1
                basisSet = load(ModelProcedure.basisset.file{bb});          % Load basis set
                basisSet = basisSet.BASIS;
                basisSet = recalculateBasisSpecs(basisSet);                 % Add ppm axis and frequency domain data
                basisSet = fit_sortBasisSet(basisSet);                      % Sort according to Osprey standard
            else
                basisSetToAdd = load(ModelProcedure.basisset.file{bb});     % Load basis set
                basisSetToAdd = basisSetToAdd.BASIS;
                basisSetToAdd = recalculateBasisSpecs(basisSetToAdd);       % Add ppm axis and frequency domain data
                basisSetToAdd = fit_sortBasisSet(basisSetToAdd);            % Sort according to Osprey standard
                basisSet.fids = cat(3,basisSet.fids,basisSetToAdd.fids);    % Concatenate time domain basis functions for 2D fit
                basisSet.specs = cat(3,basisSet.specs,basisSetToAdd.specs); % Concatenate frequency domain basis functions for 2D fit
            end
        end
        basisSet.sz = size(basisSet.fids);                                  % Recalculate size entry
        basisSet.nExtra = basisSet.sz(3);                                   % Update extra dimension
    end
    if isfield(ModelProcedure.basisset,'opts')                              % Apply options to basis, e.g. repeat for averages
        if isfield(ModelProcedure.basisset.opts,'repeat')                   % Repeat for each average
            basisSet.fids = repmat(basisSet.fids,[1 1 ModelProcedure.basisset.opts.repeat]);    % Repeat time domain basis functions for 2D fit
            basisSet.specs = repmat(basisSet.specs,[1 1 ModelProcedure.basisset.opts.repeat]);  % Repeat frequency domain basis functions for 2D fit
            basisSet.sz = size(basisSet.fids);                              % Recalculate size entry
            basisSet.nExtra = basisSet.sz(3);                               % Update extra dimension
        end        
    end
    
else
    basisSet = BasisSetStruct;                                              % Take user supplied basis set
    basisSet.nExtra = 0;                                                    % Normally has no extra dimension
end

%% 3. Prepare data according to model json and model data
% Prepare data for each fit step, again, including the indirect dimensions
% for MSM
% Create the spline basis functions for the given resolution, fit range,
% and knot spacing parameter.
if ~iscell(DataToModel)                                                       % Just a single file
    DataToModel{kk} = DataToModel;                                            % Default input is a cell array so we need to change it 
end


                                                                                   % Cell array of data which we will loop over
for kk = 1 : length(DataToModel)  
    if ~isstruct(DataToModel{kk})                                                  % It is not a FID-A struct?
        if contains(DataToModel{kk},'.nii')                                        % We do only accept .nii files for direct parsing  
                temp                = io_loadspec_niimrs(DataToModel{kk});          % Temporarily load the NIfTI-MRS file
                if average                                                        % Average before modeling
                    temp            = op_averaging(temp);                         % Average before modeling
                else
                    temp.dims.extras = temp.dims.averages;
                    temp.dims.averages = 0;                                       % Multi-spectral modelign only works in extra dim                   
                end
                DataToModel{kk} 	= temp;
        end            
    end          
    for ss = 1 : length(ModelProcedure.Steps)
        % Setup options from model procedure json file
        opts.ModelFunction      = ModelProcedure.Steps{ss}.ModelFunction;               % get model function
        opts.baseline             = ModelProcedure.Steps{ss}.fit_opts.baseline;         % setup baseline options
        opts.optimDomain        = ModelProcedure.Steps{ss}.fit_opts.optimDomain;        % do the least-squares optimization in the frequency domain
        opts.optimSignalPart    = ModelProcedure.Steps{ss}.fit_opts.optimSignalPart;    % do the least-squares optimization over the real part of the spectrum
        opts.optimFreqFitRange  = ModelProcedure.Steps{ss}.fit_opts.ppm;                % set the frequency-domain fit range to the fit range specified in the Osprey container
        opts.solver  = ModelProcedure.Steps{ss}.fit_opts.solver;                        % set the solver specified in the Osprey container
        opts.NumericJacobian = NumericJacobian;                                         % Use numerical jacobian instead of the analytical jacobian
        opts.CheckGradient = CheckGradient;                                             % Do a gradient check in the lsqnonlin solver
    
        if isfield(ModelProcedure.Steps{ss}.fit_opts,'InitialPick')
            opts.InitialPick  = ModelProcedure.Steps{ss}.fit_opts.InitialPick;          % Do inital fit on specific spectrum
        end
        if isfield(ModelProcedure.Steps{ss},'parameter')
            opts.parameter  = ModelProcedure.Steps{ss}.parameter;                       % specify parmetrizations constructor
        end
        if isfield(ModelProcedure.Steps{ss},'parametrizations')
            opts.parametrizations  = ModelProcedure.Steps{ss}.parametrizations;         % specify parmetrizations constructor
            parameter = {'ph0','ph1','gaussLB','lorentzLB','freqShift','metAmpl','baseAmpl'};
            opts.parametrizations  = orderfields(opts.parametrizations,parameter(ismember(parameter, fieldnames(opts.parametrizations)))); % order struct names according to standard
        end
        if ModelProcedure.Steps{ss}.extra.flag == 1 && isfield(ModelProcedure.Steps{ss}.extra,'DynamicModelJson')
            opts.paraIndirect  = jsonToStruct(ModelProcedure.Steps{ss}.extra.DynamicModelJson); % specify parmetrizations constructor for indirect dimension
        end
        if strcmp(ModelProcedure.Steps{ss}.module,'OptimReg')                           % Change tolerance values (good for regularization)
            opts.FunctionTolerance = 1e-3;
            opts.StepTolerance = 1e-3;
            opts.OptimalityTolerance = 1e-3;
        else
            opts.FunctionTolerance = 1e-6;
            opts.StepTolerance = 1e-6;
            opts.OptimalityTolerance = 1e-6;
        end
        
    
            clc
            fprintf('Running model procedure step %i. \n', ss);
            % Apply zero-filling if needed
            if zf &&  ~DataToModel{kk}.flags.zeropadded
                 DataToModel{kk} = op_zeropad(DataToModel{kk},2);                    % Zero-fill data if needed (real part optimization)
            end
    
            % Create an instance of the class
            if ss == 1
                ModelParameter{kk,1} = FitObject(DataToModel{kk}, basisSet, opts);  % Create OspreyFitObj instance (ss = 1)
            else
                ModelParameter{kk,1}.updateOptsAccordingToStep(opts);               % Update model options according to step
            end
    
            % Set basisset according to step
            ModelParameter{kk,1}.excludeBasisFunctionFromFit('all');
            ModelParameter{kk,1}.includeBasisFunctionInFit(ModelProcedure.Steps{ss}.basisset.include);
    
            % Is this modeling or optimization of the regularization parameter
            if ~strcmp(ModelProcedure.Steps{ss}.module,'OptimReg')
                fprintf('Running model of dataset #%i of %i \n', kk, length(DataToModel));
                % Run steps defined ModelProcedure struct
                % Loop across all steps with anonymous calls defined in the ModelProcedure struct
                ModelParameter{kk}.createModel;
            else
                fprintf('Optimize regularizer of dataset #%i of %i \n', kk, length(DataToModel));
                % Run steps defined in ModelProcedure struct N times to
                % optimize the regularization parameter                
                opts.regularizer    = ModelProcedure.Steps{ss}.regularizer;
                if ss == 1
                    ModelParameter{kk,1}.optimizeRegularization(opts);          
                else
                    ModelParameter{kk,1}.updateBaselineAccordingToStep;
                    ModelParameter{kk,1}.updateOptsAccordingToStep(opts);
                    ModelParameter{kk,1}.optimizeRegularization(opts);
                end
            end
    end
    DataToModel{kk} 	= [];                                               % Memory efficiency
end
    %% 4. Save parameter results (and export as NII)
    % Save ModelParameter to be includeded in the container and final
    % results in NII LCM format for easy reading and visualization

    % Lets reduce the file size    
    for kk = 2 : size(ModelParameter,1)
        ModelParameter{kk,1}.economizeStorage(1,1);                         % Remove basis set and jacobians
    end
    ModelParameter{1,1}.economizeStorage(0,1);                              % Keep first basis set but no jacobians
end

function isStructOrCell(f)
    assert(iscell(f) || isstruct(f));
end

function out_str = ConvertRelativePath(in_str)
    if strcmp(in_str(1:5),'which')                                          % Full path is not given we have to eval
        strdef = append('out_str = ', in_str, ';');
        pattern = "''"; % Match find double quotes
        replacement = "'"; % Replace the matched string with just a quote
        strdef = regexprep(strdef, pattern, replacement);
        eval(strdef);
    else
        out_str = in_str;
    end
    if ~isfile(out_str) && ~isempty(out_str)                                % Intercept if file doesn't exist
        error('basis file %s does not exist.', out_str);
    end
end