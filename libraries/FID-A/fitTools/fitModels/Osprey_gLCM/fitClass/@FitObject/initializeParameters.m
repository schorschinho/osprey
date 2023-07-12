function [init, lb, ub, ex, sd, fun] = initializeParameters(obj, init, lb, ub, ex, sd, fun, parameter)
%%  [init, lb, ub, ex, sd, fun] = initializeParameters(obj, init, lb, ub, ex, sd, fun, parameter)
%   This method sets the parametrization of the OspreyFitObj according to
%   the model procedure step as well as the indirect dimension link
%
%   USAGE:
%       [init, lb, ub, sd, fun] = initializeParameters(obj, init, lb, ub, sd, fun, parameter)
%
%   INPUTS:
%       init      = struct with initial parameters that is build for all parameters
%       lb        = struct with lower bounds that is build for all parameters
%       ub        = struct with upper bounds that is build for all parameters
%       ex        = struct with expectation values that is build for all parameters
%       sd        = struct with standard deviations that is build for all parameters
%       fun       = struct with function type parameters that is build for all parameters
%       parameter = parameter to initialize (e.g. p0, ph1, gaussLB...)
%
%   OUTPUTS:
%       obj       = OspreyFitObj with updated parametrizations.
%       init      = struct with initial parameters that is build for all parameters
%       lb        = struct with lower bounds that is build for all parameters
%       ub        = struct with upper bounds that is build for all parameters
%       ex        = struct with expectation values that is build for all parameters
%       sd        = struct with standard deviations that is build for all parameters
%       fun       = struct with function type parameters that is build for all parameters
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Write parametrizations in obj

nBasisFcts = sum(obj.BasisSets.includeInFit(obj.step,:));               % Number of basis functions
nBaselineComps = size(obj.BaselineBasis, 2);                            % Number of baseline parameters
if ~isfield(obj.Options{obj.step},'parametrizations')                   % No parametrizations field use defaults
    obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter);  % Write parameter defaults in object
else if ~isfield(obj.Options{obj.step}.parametrizations,parameter)      % No parametrization for the parameter
        obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter); % Write parameter defaults in object
else
    obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter,obj.Options{obj.step}.parametrizations.(parameter)); % Write parameter according to model procedure step
end
end
%% Define the number of parameters per spectrum and per indirect dimension for the parameter
switch parameter                            % Switch with parameter names
    case {'ph0', 'ph1', 'gaussLB'}
        % one parameter per spectrum
        nParamsPerSpec = 1;
    case {'metAmpl', 'freqShift', 'lorentzLB'}
        % one parameter per metabolite per spectrum
        nParamsPerSpec = nBasisFcts;
    case {'baseAmpl'}
        % one parameter per spline function per spectrum
        nParamsPerSpec = nBaselineComps;
end
nParamsPerIndir = size(obj.Data.fids,2);   % Assume that the parameters have to repeated along the indirect dimension to start with
type = 'free';                             % Assume that the indirect dimension parametrization is free to strat with

%% Update the parameterization according to the indirect dimension json file

if isfield(obj.Options{obj.step}, 'paraIndirect')   % json file for 2D model supplied
    if isfield(obj.Options{obj.step}.paraIndirect.parameters,parameter) % No parametrizations for parameter field use defaults
        if ~isfield(obj.Options{obj.step}.paraIndirect.parameters.(parameter),'parametrizations')  % No parametrizations for parameter field use defaults
            if ~isfield(obj.Options{obj.step}.paraIndirect.parameters,parameter)                   % No parametrizations for parameter field use defaults
                obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter); % Write parameter defaults in object
            else
                obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter,obj.Options{obj.step}.paraIndirect.parameters.(parameter)); % Write parameter according to 2D model procedure step
            end
        else
            obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter,obj.Options{obj.step}.paraIndirect.parameters.(parameter).parametrizations); % Write parameter according to 2D model procedure step
        end
        switch obj.Options{obj.step}.paraIndirect.parameters.(parameter).type   % Switch with parametrization function
            case 'fixed'
                nParamsPerIndir = 1;                                            % Parameters are linked across indirect dim
            case 'free'
                nParamsPerIndir = size(obj.Data.fids,2);                        % Parameters are free across indirect dim
            case 'dynamic'
                nParamsPerIndir = 1;
                obj.Options{obj.step}.parametrizations.(parameter).modulator = obj.Options{obj.step}.paraIndirect.parameters.(parameter).modulator;             %Add modulator for parameterization e.g. time
                obj.Options{obj.step}.parametrizations.(parameter).parameterNames = obj.Options{obj.step}.paraIndirect.parameters.(parameter).parameterNames;   %Add new parameter names
        end
    end
    type = obj.Options{obj.step}.paraIndirect.parameters.(parameter).type;  % Overwrite inital assuption of free links
end

%% Update the structs with the init, lb, ub, ex, sd, fun entries
% This is used to generate the init, lb, ub, ex, sd, fun struct according to the
% number of parameters needed. This is based on the nParamsPerSpec and
% nParamsPerIndir values. We basically have to repeat a single value
% according to those values. This gets more complicated when we are using
% fit ressults from a prior step. For example we can have lorentzianLB
% initial values from a single spectrum fit that we now have to repeat
% along the indirect dimension
if ~iscell(obj.Options{obj.step}.parametrizations.(parameter).init)  % This is numeric if defaults or not pulled from a prior step
    if ~strcmp(type,'dynamic')      % For free or fixed cases we just have to repeat the parameter value according to nParamsPerSpec and nParamsPerIndir
        structs = {'init','lb','ub','ex','sd'};  % Structs to write for output
        for st = 1 : length(structs)    % Loop over structs
            switch num2str([size(obj.Options{obj.step}.parametrizations.(parameter).(structs{st}),1) == nParamsPerIndir ...     % Check dimensions and repeat parameter accordingly
                    size(obj.Options{obj.step}.parametrizations.(parameter).(structs{st}),2) == nParamsPerSpec])
                case '1  1' % Correct dimensions do nothing
                    eval([ structs{st} '.' parameter ' = obj.Options{obj.step}.parametrizations.' parameter '.' structs{st} ';'])
                case '0  0' % Both dimensions are wrong repeat according to nParamsPerSpec and nParamsPerIndir
                    eval([ structs{st} '.' parameter ' = squeeze(repmat(obj.Options{obj.step}.parametrizations.' parameter '.' structs{st} ', [nParamsPerIndir, nParamsPerSpec]));'])
                case '0  1' % Indirect dimensions not correct repeat nParamsPerIndir times
                    eval([ structs{st} '.' parameter ' = squeeze(repmat(obj.Options{obj.step}.parametrizations.' parameter '.' structs{st} ', [nParamsPerIndir, 1]));'])
                case '1  0' % Params per spectrum dimensions not correct repeat nParamsPerSpec times
                    eval([ structs{st} '.' parameter ' = squeeze(repmat(obj.Options{obj.step}.parametrizations.' parameter '.' structs{st} ', [1, nParamsPerSpec]));'])
            end
        end
    else % For dynamic parametrization lb and ub are already defined from the set_parameter function
        if size(obj.Options{obj.step}.parametrizations.(parameter).init,2) == nParamsPerSpec    % Dimensions are correct
            init.(parameter) = obj.Options{obj.step}.parametrizations.(parameter).init;         % Add initials
            lb.(parameter) = obj.Options{obj.step}.parametrizations.(parameter).lb;             % Add upper bounds
            ub.(parameter) = obj.Options{obj.step}.parametrizations.(parameter).ub;             % Add lower bounds
            ex.(parameter) = obj.Options{obj.step}.parametrizations.(parameter).ex;             % Add expectation values
            sd.(parameter) = obj.Options{obj.step}.parametrizations.(parameter).sd;             % Add standard deviations
        else %  Both dimensions are wrong repeat according to nParamsPerSpec and nParamsPerIndir
            init.(parameter) = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).init, [nParamsPerIndir, nParamsPerSpec])); % Add initials
            lb.(parameter) = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).lb, [nParamsPerIndir, nParamsPerSpec])); % Add upper bounds
            ub.(parameter) = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).ub, [nParamsPerIndir, nParamsPerSpec])); % Add lower bounds
            ex.(parameter) = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).ex, [nParamsPerIndir, nParamsPerSpec])); % Add expectation values
            sd.(parameter) = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).sd, [nParamsPerIndir, nParamsPerSpec])); % Add standard deviations
        end
    end
else % When this is a cell array it contains the 'Step' reference (e.g. Step 2) to use results from a prior step
    step_to_get_ini = split(obj.Options{obj.step}.parametrizations.(parameter).init{1},' ');        % Get entry from model procedure e.g. Step 1
    step_to_get_ini = str2num(step_to_get_ini{2});                                                  % Get numeric value of step
    init.(parameter) = squeeze(repmat(obj.Model{step_to_get_ini}.parsOut.(parameter), [nParamsPerIndir, nParamsPerSpec]));  % Update init value accordingly
    lb.(parameter)   = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).lb,   [nParamsPerIndir , nParamsPerSpec])); % Repeat lb value according to nParamsPerIndir nParamsPerSpec
    ub.(parameter)   = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).ub,   [nParamsPerIndir , nParamsPerSpec])); % Repeat ub value according to nParamsPerIndir nParamsPerSpec
    ex.(parameter)   = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).ex,   [nParamsPerIndir , nParamsPerSpec])); % Repeat ex value according to nParamsPerIndir nParamsPerSpec
    sd.(parameter)   = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).sd,   [nParamsPerIndir , nParamsPerSpec])); % Repeat sd value according to nParamsPerIndir nParamsPerSpec
end
if isfield(obj.Options{obj.step}.parametrizations.(parameter),'type')              % If 2D we will have a type entry for the relation
    fun.(parameter) = obj.Options{obj.step}.parametrizations.(parameter).type;     % Add parametrization function type
else
    obj.Options{obj.step}.parametrizations.(parameter).type = type;               % If 2D we will have a type entry for the relation
    fun.(parameter) = obj.Options{obj.step}.parametrizations.(parameter).type;      % Add parametrization function type
end
end

% Description of default parameters are described below
function parametrizations = set_parameter(obj,parameter,predefined)
% This functions defines the parametrizations according to default values
% or predefined valeus parsed in the predefined variable
switch parameter                                            % Parameter switch with default parameters
    case 'ph0'
        % Initialize phi0 as constant with value 0
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = -pi;
        parametrizations.ub      = pi;
        parametrizations.init    = 0;
        parametrizations.ex      = 0;
        parametrizations.sd      = Inf;
        parametrizations.RegFun  = '';

    case 'ph1'
        % Initialize phi1 as constant with value 0
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = -pi/30;
        parametrizations.ub      = pi/30;
        parametrizations.init    = 0;
        parametrizations.ex      = 0;
        parametrizations.sd      = Inf;
        parametrizations.RegFun  = '';

    case 'gaussLB'
        % Initialize Gaussian LB as constant with value [0.04 *
        % hz/ppm]
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = 0;
        parametrizations.ub      = Inf;
        parametrizations.init    = 0.04 * obj.Data.txfrq*1e-6;
        parametrizations.ex      = 0.04 * obj.Data.txfrq*1e-6;
        parametrizations.sd      = Inf;
        parametrizations.RegFun  = '';

    case 'lorentzLB'
        % Initialize Lorentzian LB as constant with value
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = 0;
        parametrizations.ub      = Inf;
        parametrizations.init    = 0;
        parametrizations.ex      = 0;
        parametrizations.sd      = Inf;
        parametrizations.RegFun  = '';

    case 'freqShift'
        % Initialize frequency shifts as constant with value 0 Hz
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = -5;
        parametrizations.ub      = 5;
        parametrizations.init    = 0;
        parametrizations.ex      = 0;
        parametrizations.sd      = Inf;
        parametrizations.RegFun  = '';

    case 'metAmpl'
        % Initialize metabolite amplitudes as free with value 0
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = 0;
        parametrizations.ub      = Inf;
        parametrizations.init    = 0;
        parametrizations.ex      = 0;
        parametrizations.sd      = Inf;
        parametrizations.RegFun  = '';

    case 'baseAmpl'
        % Initialize baseline amplitudes as free with value 0
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = -Inf;
        parametrizations.ub      = Inf;
        parametrizations.init    = 0;
        parametrizations.ex      = 0;
        parametrizations.sd      = Inf;
        parametrizations.RegFun  = '';

    case 'x'
        % Initialize x (the external dependency vector) as natural numbers
        parametrizations.modulator = [1:1:size(obj.Data.fids,2)];
        parametrizations.name   = 'independentVariable';
        parametrizations.RegFun    = '';

    otherwise
        % Initialize as free parameter
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = -Inf;
        parametrizations.ub      = Inf;
        parametrizations.init    = 0;
        parametrizations.ex      = 0;
        parametrizations.sd      = Inf;
        parametrizations.RegFun  = '';

end

default_parametrization = parametrizations;                             % Store default parameterization as backup

% Update the parametrization according to the predefiend variable if it
% is defined
if nargin == 3
    pars = fields(predefined);                                          % Fields in predefined struct e.g. ub, lb, init...
    for ff = 1 : length(pars)                                           % Loop over fields in predefined struct
        if ~isstr(predefined.(pars{ff})) && ~iscell(predefined.(pars{ff}))  % Is a numeric value so we can just take it
            parametrizations.(pars{ff}) = predefined.(pars{ff});            % Add numeric value to parametrization
        else if ~iscell(predefined.(pars{ff}))                          % Is not a cell but a string this is need for Inf and -Inf parsing from json
                if ~isempty(str2num(predefined.(pars{ff})))         % Is Inf or -Inf string
                    parametrizations.(pars{ff}) = str2num(predefined.(pars{ff})); % Inf or -Inf to value
                else
                    parametrizations.(pars{ff}) = predefined.(pars{ff}); % Is a numeric value
                end
        else                                                        % If this is a cell it could be a 'Step' or 'Regularizer'
            for pp = 1 : length(predefined.(pars{ff}))          % Loop over entries in the predefined struct this is needed if we add different values per basis function
                if ischar(predefined.(pars{ff}){pp})            % If string it is either a Step reference or a Regularizer reference
                    if ~strcmp(predefined.(pars{ff}){pp}(1:2),'St')
                        parametrizations.(pars{ff})(pp,1) = str2num(predefined.(pars{ff}){pp});
                    else
                        step_to_get_ini = split(predefined.(pars{ff}){pp},' ');         % Get entry from model procedure e.g. Step 1
                        step_to_get_ini = str2num(step_to_get_ini{2});                  % Get numeric value of step
                        if ~strcmp(pars{ff}(1:2),'Re')                                  % Is Step
                            parametrizations.(pars{ff}) = squeeze(obj.Model{step_to_get_ini}.parsOut.(parameter));  % Update numerical value accordingly
                        else                                                            % Is Regularizer
                            parametrizations.(pars{ff}) = obj.Model{step_to_get_ini}.Regularization.OptimalRegPar; % Get regularization parameter
                        end
                    end
                else
                    if size(parametrizations.(pars{ff}),2) == 1                         % If single entry we will just add it
                        parametrizations.(pars{ff})(pp,1) = predefined.(pars{ff}){pp};  % Update numerical value accordingly
                    else                                                                % One value per basis function
                        parametrizations.(pars{ff})(pp,:) = repmat(predefined.(pars{ff}){pp},[1, size(parametrizations.(pars{ff}),2)]); % Repeat according to numbers of basis functions
                    end
                end
            end
        end
        end
    end                                                                 % End loop over fields in predefined struct
end

% This is needed if the number of parameters (basis functions or
% baseline parameters) changes between fit steps
if obj.step > 1                                                         % Only needed for steps > 1
    if strcmp(parameter, 'lorentzLB') || strcmp(parameter, 'freqShift') % Only do this for non amplitude parameters
        if sum(obj.BasisSets.includeInFit(obj.step,:)) ~= sum(obj.BasisSets.includeInFit(obj.step-1,:))     % Did the number of parameters change
            if sum(obj.BasisSets.includeInFit(obj.step,:)) > sum(obj.BasisSets.includeInFit(obj.step-1,:))  % More parameters then before
                for ff = 1 : length(pars)
                    if strcmp(pars{ff},'ub') || strcmp(pars{ff},'lb') || strcmp(pars{ff},'init')            % Only apply to init, ub, lb
                        parametrizations.(pars{ff}) = cat(2,parametrizations.(pars{ff}),...
                            default_parametrization.(pars{ff}) * ...
                            ones(1,sum(obj.BasisSets.includeInFit(obj.step,:))- sum(obj.BasisSets.includeInFit(obj.step-1,:))));
                    end
                end
            end
        end
    end
end
end