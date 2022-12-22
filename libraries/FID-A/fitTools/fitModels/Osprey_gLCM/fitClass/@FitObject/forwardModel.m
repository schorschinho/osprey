function [Y, baseline, metabs] = forwardModel(x, basisSet, baselineBasis, ppm, t)
            
    % Define the default 1-D forward model first
    fids = basisSet.fids;
    nBasisFcts = sum(basisSet.includeInFit);
    
    inputParams     = FitObject.x2pars(x, nBasisFcts);
    gaussLB     = inputParams.gaussLB;
    lorentzLB   = inputParams.lorentzLB;
    freqShift   = inputParams.freqShift;
    metAmpl     = inputParams.metAmpl;
    baseAmpl    = inputParams.baseAmpl;
    ph0         = inputParams.ph0;
    ph1         = inputParams.ph1;
    
    timeDomainMultiplier = zeros(size(fids));
    for ll = 1:nBasisFcts
        timeDomainMultiplier(:,ll) = exp(-(1i*freqShift(ll) + lorentzLB(ll) + gaussLB.^2.*t).*t)';    
    end
    
    Fl = timeDomainMultiplier .* fids;
    specs = fftshift(fft(Fl, [], 1),1);
    mets = metAmpl' * specs';
    baseline = baseAmpl' * baselineBasis';
    
    T_ph = exp(-1j .* (ph0 + ph1.*ppm));
    Y = T_ph .* (mets + baseline);  
    Y = Y';
    
    metabs = repmat(T_ph', [1 nBasisFcts]) .* specs .* metAmpl';
end