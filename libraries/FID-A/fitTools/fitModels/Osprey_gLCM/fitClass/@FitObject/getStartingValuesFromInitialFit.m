function [init] = getStartingValuesFromInitialFit(obj, init, parameter)
    
    % Collect the number of basis functions included in the fit
    nBasisFcts = sum(obj.BasisSets.includeInFit);
    
    if strcmp(obj.Options.parametrizations.(parameter).fun, 'free')
        switch parameter
            case {'ph0', 'ph1', 'gaussLB', 'metAmpl', 'baseAmpl'}
                % these were all determined during the initial fit
                init.(parameter) = obj.Model{1}.parsOut.(parameter)';
            case {'freqShift', 'lorentzLB'}
                % duplicate the initial value by the number of
                % metabolite basis functions
                init.(parameter) = repmat(obj.Model{1}.parsOut.(parameter), [1, nBasisFcts]);
        end
        
    else
        error('Continue coding here when you start to use a parametrization');
    end
    
end  