function [init, lb, ub, fun] = initializeParameters(obj, init, lb, ub, fun, parameter)
    
    % Collect number of basis functions and baseline spline elements
    nBasisFcts = sum(obj.BasisSets.includeInFit(obj.step,:));
    nBaselineComps = size(obj.BaselineBasis, 2);
    if ~isfield(obj.Options{obj.step},'parametrizations')
        obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter);
    else if ~isfield(obj.Options{obj.step}.parametrizations,parameter)
            obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter);
        else
            obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter,obj.Options{obj.step}.parametrizations.(parameter));
        end
    end
    if strcmp(obj.Options{obj.step}.parametrizations.(parameter).fun, 'free')
        switch parameter
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
        nParamsPerIndir = size(obj.Data.fids,2);  
        type = 'free';
        if isfield(obj.Options{obj.step}, 'paraIndirect')   % json file for 2D model supplied
            if isfield(obj.Options{obj.step}.paraIndirect.parameters,parameter)
                if ~isfield(obj.Options{obj.step}.paraIndirect.parameters.(parameter),'parametrizations')   % Update parametrization
                    obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter);
                else
                    obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter,obj.Options{obj.step}.paraIndirect.parameters.(parameter).parametrizations);
                end
                type = obj.Options{obj.step}.paraIndirect.parameters.(parameter).type; % Save parametrization 
                switch type
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
        end
        if ~iscell(obj.Options{obj.step}.parametrizations.(parameter).init)
            init.(parameter) = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).init, [nParamsPerIndir, nParamsPerSpec]));
        else
            step_to_get_ini = split(obj.Options{obj.step}.parametrizations.(parameter).init{1},' ');
            step_to_get_ini = str2num(step_to_get_ini{2});
            init.(parameter) = squeeze(repmat(obj.Model{step_to_get_ini}.parsOut.(parameter), [nParamsPerIndir, nParamsPerSpec]));
        end
        lb.(parameter)   = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).lb,   [nParamsPerIndir , nParamsPerSpec]));
        ub.(parameter)   = squeeze(repmat(obj.Options{obj.step}.parametrizations.(parameter).ub,   [nParamsPerIndir , nParamsPerSpec]));
        fun.(parameter) = type;
    else
        error('Continue coding here when you start to use a parametrization');
    end
    
end

function parametrizations = set_parameter(obj,parameter,predefined)
    switch parameter
        case 'ph0'
        % Initialize phi0 as constant with value 0
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = -pi;
        parametrizations.ub      = pi;
        parametrizations.init    = 0;
        
        case 'ph1'
        % Initialize phi1 as constant with value 0
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = -pi/30;
        parametrizations.ub      = pi/30;
        parametrizations.init    = 0;
        
        case 'gaussLB'
        % Initialize Gaussian LB as constant with value [0.04 *
        % hz/ppm]
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = 0;
        parametrizations.ub      = log(0.04 * obj.Data.txfrq*1e-6);
        parametrizations.init    = 0.04 * obj.Data.txfrq*1e-6;
        
        case 'lorentzLB'
        % Initialize Lorentzian LB as constant with value 2 Hz
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = 0;
        parametrizations.ub      = Inf;
        parametrizations.init    = log(2);
        
        case 'freqShift'
        % Initialize frequency shifts as constant with value 0 Hz
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = -5;
        parametrizations.ub      = 5;
        parametrizations.init    = 0;
        
        case 'metAmpl'
        % Initialize metabolite amplitudes as free with value 0
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = 0;
        parametrizations.ub      = Inf;
        parametrizations.init    = 0;
        
        case 'baseAmpl'
        % Initialize baseline amplitudes as free with value 0
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = -Inf;
        parametrizations.ub      = Inf;
        parametrizations.init    = 0;
        
        case 'x'
        % Initialize x (the external dependency vector) as natural numbers
        parametrizations.modulator = [1:1:size(obj.Data.fids,2)];
        parametrizations.name   = 'independentVariable';

        otherwise
        % Initialize as free parameter
        parametrizations.fun     = 'free';
        parametrizations.gradfun = 'free';
        parametrizations.lb      = -Inf;
        parametrizations.ub      = Inf;
        parametrizations.init    = 0;

    end
    if nargin == 3
        pars = fields(predefined);
        for ff = 1 : length(pars)
            if ~isstr(predefined.(pars{ff})) && ~iscell(predefined.(pars{ff}))
                parametrizations.(pars{ff}) = predefined.(pars{ff});
            else if ~iscell(predefined.(pars{ff}))
                parametrizations.(pars{ff}) = str2num(predefined.(pars{ff}));
            else
                for pp = 1 : length(predefined.(pars{ff}))
                    parametrizations.(pars{ff})(pp,1) = str2num(predefined.(pars{ff}){pp});
                end
            end
            end
        end
    end
    if ~strcmp(parametrizations.fun, 'free') && ~strcmp(parametrizations.fun, 'fixed')              % Dynamic parameter across indirect dim
        pars = {'lb', 'ub', 'init'};
        newPars = max([length(parametrizations.lb) length(parametrizations.ub) length(parametrizations.init)]); %Add new parameter names
        for ff = 1 : length(pars)
            if length(parametrizations.(pars{ff})) < newPars
                parametrizations.(pars{ff}) = repmat(parametrizations.(pars{ff}), [newPars, 1]);
            end
        end
    end
end