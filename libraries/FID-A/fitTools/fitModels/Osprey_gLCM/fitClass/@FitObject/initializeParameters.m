function [init, lb, ub] = initializeParameters(obj, init, lb, ub, parameter)
    
    % Collect number of basis functions and baseline spline elements
    nBasisFcts = sum(obj.BasisSets.includeInFit);
    nBaselineComps = size(obj.BaselineBasis, 2);
    if ~isfield(obj.Options{obj.step}.parametrizations,parameter)
        obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter);
    else
        obj.Options{obj.step}.parametrizations.(parameter) = set_parameter(obj,parameter,obj.Options{obj.step}.parametrizations.(parameter));
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
        if ~iscell(obj.Options{obj.step}.parametrizations.(parameter).init)
            init.(parameter) = repmat(obj.Options{obj.step}.parametrizations.(parameter).init, [size(obj.Data.fids,2), nParamsPerSpec]);
        else
            step_to_get_ini = split(obj.Options{obj.step}.parametrizations.(parameter).init{1},' ');
            step_to_get_ini = str2num(step_to_get_ini{2});
            init.(parameter) = repmat(obj.Model{step_to_get_ini}.parsOut.(parameter), [size(obj.Data.fids,2), nParamsPerSpec]);
        end
        lb.(parameter)   = repmat(obj.Options{obj.step}.parametrizations.(parameter).lb,   [size(obj.Data.fids,2), nParamsPerSpec]);
        ub.(parameter)   = repmat(obj.Options{obj.step}.parametrizations.(parameter).ub,   [size(obj.Data.fids,2), nParamsPerSpec]);
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
        parametrizations.ub      = 10;
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
        parametrizations.values = [1:1:size(obj.Data.fids,2)];
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
            if ~isstr(predefined.(pars{ff}))
                parametrizations.(pars{ff}) = predefined.(pars{ff});
            else
                parametrizations.(pars{ff}) = str2num(predefined.(pars{ff}));
            end
        end
    end
end