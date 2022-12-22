function [init, lb, ub] = initializeParameters(obj, init, lb, ub, parameter)
    
    % Collect number of basis functions and baseline spline elements
    nBasisFcts = sum(obj.BasisSets.includeInFit);
    nBaselineComps = size(obj.BaselineBasis, 2);
    
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
        if obj.step == 1
            init.(parameter) = repmat(obj.Options{obj.step}.parametrizations.(parameter).init, [size(obj.Data.fids,2), nParamsPerSpec]);
            lb.(parameter)   = repmat(obj.Options{obj.step}.parametrizations.(parameter).lb,   [size(obj.Data.fids,2), nParamsPerSpec]);
            ub.(parameter)   = repmat(obj.Options{obj.step}.parametrizations.(parameter).ub,   [size(obj.Data.fids,2), nParamsPerSpec]);
        else
            if ~iscell(obj.Options{obj.step}.initials.(parameter))
                init.(parameter) = repmat(obj.Options{obj.step}.initials.(parameter), [size(obj.Data.fids,2), nParamsPerSpec]);
            else
                step_to_get_ini = split(obj.Options{obj.step}.initials.(parameter){1},' ');
                step_to_get_ini = str2num(step_to_get_ini{2});
                init.(parameter) = repmat(obj.Model{step_to_get_ini}.parsOut.(parameter), [size(obj.Data.fids,2), nParamsPerSpec]);
            end
            lb.(parameter)   = repmat(obj.Options{obj.step}.parametrizations.(parameter).lb,   [size(obj.Data.fids,2), nParamsPerSpec]);
            ub.(parameter)   = repmat(obj.Options{obj.step}.parametrizations.(parameter).ub,   [size(obj.Data.fids,2), nParamsPerSpec]);
        end
    else
        error('Continue coding here when you start to use a parametrization');
    end
    
end