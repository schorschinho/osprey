function sse = initialFitLossFunction(x, tdData, basisSet, baselineBasis, ppm, t, fitRange)
            
    % get time-domain data, fft, define fit range
    fdData = fftshift(fft(tdData,[],1),1);
    [indMin, indMax] = ppmToIndex(ppm, fitRange);
    
    % get basis sets, find number of basis functions
    fids = basisSet.fids;
    nBasisFcts = size(fids,2);
    
    % apply non-linear parameters to basis set and fft
    ph0         = x(1);
    gaussLB     = exp(x(2));
    lorentzLB   = exp(x(3));
    freqShift   = x(4);
    specs       = FitObject.transformBasis(fids, gaussLB, lorentzLB, freqShift, t);
    
    % append basis set with the spline baseline basis functions
    fullBasis = cat(2, specs, baselineBasis);
    % apply zero-order phase
    fullBasis = exp(-1j*ph0) .* fullBasis;
    
    % crop to fit range
    fdData = fdData(indMin:indMax,:);
    fullBasis = fullBasis(indMin:indMax,:);
    
    % concatenate real and imaginary
    fdData = cat(1, real(fdData), imag(fdData));
    fullBasis = cat(1, real(fullBasis), imag(fullBasis));
    
    % solve the linear equation system
    beta = real(pinv(fullBasis) * fdData);
    
    % project metabolite coefficients to >0
    beta(beta(1:nBasisFcts) < 0) = 0;
    
    % make spectrum
    prediction = fullBasis * beta;
    
    % apply the forward model
    diffVec     = fdData - prediction;
    sse         = mean(abs(diffVec).*2);
    
end