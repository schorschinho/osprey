function paramStruct = x2pars(x, nBasisFcts)
            
    % Converts a 1-D x vector into a parameter struct
    paramStruct.ph0 = x(1);
    paramStruct.ph1 = x(2);
    paramStruct.gaussLB = x(3);
    paramStruct.lorentzLB = x(4:3+nBasisFcts);
    paramStruct.freqShift = x(4+nBasisFcts:3+2*nBasisFcts);
    paramStruct.metAmpl = x(4+2*nBasisFcts:3+3*nBasisFcts);
    paramStruct.baseAmpl = x(4+3*nBasisFcts:end);

    
end