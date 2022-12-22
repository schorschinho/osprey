function specs = transformBasis(fids, gaussLB, lorentzLB, freqShift, t)
    timeDomainMultiplier = zeros(size(fids));
    nBasisFcts = size(fids,2);
    for ll = 1:nBasisFcts
        timeDomainMultiplier(:,ll) = exp(-(1i*freqShift + lorentzLB + gaussLB.^2.*t).*t)';
    end
    
    Fl = timeDomainMultiplier .* fids;
    specs = fftshift(fft(Fl, [], 1), 1);
end