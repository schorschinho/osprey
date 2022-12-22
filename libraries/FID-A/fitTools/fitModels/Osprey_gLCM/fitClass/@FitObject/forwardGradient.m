function [grad, jac] = forwardGradient(x, data, basisSet, baselineBasis, ppm, t, fitRange)
            
    [indMin, indMax] = FitObject.ppmToIndex(ppm, fitRange);
    
    prediction  = FitObject.forwardModel(x, basisSet, baselineBasis, ppm, t);
    
    % Construct the Jacobian matrix of partial derivatives
    fids = basisSet.fids;
    nBasisFcts = size(fids,2);
    nBaselineComps = size(baselineBasis, 2);
    
    inputParams = FitObject.x2pars(x, nBasisFcts);
    gaussLB = inputParams.gaussLB;
    lorentzLB = inputParams.lorentzLB;
    freqShift = inputParams.freqShift;
    metAmpl = inputParams.metAmpl;
    baseAmpl = inputParams.baseAmpl;
    ph0 = inputParams.ph0;
    ph1 = inputParams.ph1;
    
    timeDomainMultiplier = zeros(size(fids));
    for ll = 1:nBasisFcts
        timeDomainMultiplier(:,ll) = exp(-(1i*freqShift(ll) + lorentzLB(ll) + gaussLB.^2.*t).*t)';    
    end
    
    T_ph = exp(-1j .* (ph0 + ph1.*ppm)');
    T_ph_basis = repmat(T_ph, [1, nBasisFcts]);
    T_ph_baseline = repmat(T_ph, [1, nBaselineComps]);
    T_t = repmat(t', [1, nBasisFcts]);
    T_tt = repmat(t'.*t', [1, nBasisFcts]);
    
    Fmet = timeDomainMultiplier .* fids;  

    
    Fmett = fftshift(fft(T_t .* Fmet, [], 1),1);
    Fmett2gauss = fftshift(fft(-2 .* gaussLB .* T_tt .* Fmet, [], 1),1);
    
    % partial derivatives
    dYdmetAmpl      = T_ph_basis .* fftshift(fft(Fmet, [], 1),1);
    dYdfreqShift    = T_ph_basis .* (-1j) .* Fmett .* metAmpl';
    dYdlorentzLB    = T_ph_basis .* (-1)  .* Fmett .* metAmpl';
    dYdgaussLB      = T_ph .* Fmett2gauss * metAmpl;
    dYdph0          = (1j) .* prediction;
    dYdph1          = (1j) .* ppm' .* prediction;
    dYdBj           = T_ph_baseline .* baselineBasis;
                
    % reduce to fit range
    dYdph0          = dYdph0(indMin:indMax,:);
    dYdph1          = dYdph1(indMin:indMax,:);
    dYdgaussLB      = dYdgaussLB(indMin:indMax,:);
    dYdlorentzLB    = dYdlorentzLB(indMin:indMax,:);
    dYdfreqShift    = dYdfreqShift(indMin:indMax,:);
    dYdmetAmpl      = dYdmetAmpl(indMin:indMax,:);
    dYdBj           = dYdBj(indMin:indMax,:);
    
    % reduce prediction to fit range;
    prediction = prediction(indMin:indMax);
    data = fftshift(fft(data,[],1),1);
    data = data(indMin:indMax);
    
    % fft data and reduce to fit range
    
    jac = cat(2, dYdph0, dYdph1, dYdgaussLB, dYdlorentzLB, dYdfreqShift, dYdmetAmpl, dYdBj);
    
    % grad = real(sum(prediction .* conj(jac) + conj(prediction) .* jac - conj(data) .* jac - data .* conj(jac),1))';
    
    % apply the chain rule
    % grad = real(sum(-(data .* conj(jac) - prediction .* conj(jac) + conj(data) .* jac - (conj(prediction) .* jac))))';
    grad = real(sum((data - prediction).*(-conj(jac)) + (-jac .* conj(data-prediction))))';

end