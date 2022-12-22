% Forward models, loss functions, Jacobian and gradient functions
        % below
function sse = lossFunction(x, data, basisSet, baselineBasis, ppm, t, fitRange)
    
    [indMin, indMax] = FitObject.ppmToIndex(ppm, fitRange);
    
    % apply the forward model
    data        = fftshift(fft(data, [], 1),1);
    prediction  = FitObject.forwardModel(x, basisSet, baselineBasis, ppm, t);
    diffVec     = data(indMin:indMax) - prediction(indMin:indMax);
    sse         = real(sum(diffVec .* conj(diffVec)));
    
end