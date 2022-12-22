function CRLB = getDiagonalElementsOfCRLB(x, nBasisFcts)
    
    % Converts a 1-D x vector into a parameter struct
    CRLB.ph0 = diag(x(1,1));
    CRLB.ph1 = diag(x(2,2));
    CRLB.gaussLB = diag(x(3,3));
    CRLB.lorentzLB = diag(x(4:3+nBasisFcts, 4:3+nBasisFcts));
    CRLB.freqShift = diag(x(4+nBasisFcts:3+2*nBasisFcts, 4+nBasisFcts:3+2*nBasisFcts));
    CRLB.metAmpl = diag(x(4+2*nBasisFcts:3+3*nBasisFcts, 4+2*nBasisFcts:3+3*nBasisFcts));
    CRLB.baseAmpl = diag(x(4+3*nBasisFcts:end, 4+3*nBasisFcts:end));
    
    
end