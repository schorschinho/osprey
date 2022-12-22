function x = pars2x(paramStruct)
            
    % Converts a parameter struct into a 1-D x vector that can be
    % passed on to solvers
    x = [paramStruct.ph0, paramStruct.ph1, paramStruct.gaussLB, ...
         paramStruct.lorentzLB, paramStruct.freqShift, ... 
         paramStruct.metAmpl, paramStruct.baseAmpl]';
    
end