function [ph0, gaussLB, lorentzLB, freqShift, metAmpl, baseAmpl] = initFit(obj)
        
        % This method runs a quick and dirty fit to get decent starting values.
        
        % Collect the basis functions
        step            = obj.step + 1;
        basisSet        = obj.BasisSets;
        baselineBasis   = obj.BaselineBasis;
        data            = obj.Data.fids;
        ppm             = obj.Data.ppm;
        t               = obj.Data.t;
        fitRange        = obj.Options{1}.optimFreqFitRange;
        [indMin, indMax] = ppmToIndex(ppm, fitRange);
        
        % Only use basis functions that are included
        basisSet.fids   = basisSet.fids(:, logical(basisSet.includeInFit));
        
        % Create x0
        % ph0, gaussLB, lorentzLB, freqShift
        gaussInit = 0.04*obj.Data.txfrq * 1e-6;
        lorentzInit = 2;
        x0 = [0,   log(gaussInit), log(lorentzInit), 0];
        lb = [-pi, 0,              0,                -Inf];
        ub = [pi,  Inf,            Inf,              Inf];
        
        
        tstart=tic;
        % Prepare the function wrapper
        fcn  = @(x) FitObject.initialFitLossFunction(x, data, basisSet, baselineBasis, ppm, t, fitRange);
        xk = fmincon(fcn, x0, [], [], [], [], lb, ub);
        
        time = toc(tstart);
        
        % Translate output
        fids = basisSet.fids;
        nBasisFcts = size(fids, 2);
        nBaselineComps = size(baselineBasis, 2);
        
        ph0 = xk(1);
        gaussLB = exp(xk(2));
        lorentzLB = exp(xk(3));
        freqShift = xk(4);
        
        specs = FitObject.transformBasis(fids, gaussLB, lorentzLB, freqShift, t);
        
        % append basis set with the spline baseline basis functions
        fullBasis = cat(2, specs, baselineBasis);
        % phase
        fullBasis = exp(-1j*ph0) .* fullBasis;
        
        % crop to fit range
        data = fftshift(fft(data, [], 1), 1);
        dataCrop = data(indMin:indMax,:);
        fullBasisCrop = fullBasis(indMin:indMax,:);
        
        % concatenate real and imaginary
        dataCrop = cat(1, real(dataCrop), imag(dataCrop));
        fullBasisCrop = cat(1, real(fullBasisCrop), imag(fullBasisCrop));
        
        % solve the linear equation system
        beta = real(pinv(fullBasisCrop) * dataCrop);
        
        % project metabolite coefficients to >=0
        beta(beta(1:nBasisFcts) < 0) = 0;
        
        % make spectrum
        prediction = fullBasis * beta;
        
        % separate metAmpl and baseAmpl
        metAmpl = beta(1:nBasisFcts);
        baseAmpl = beta(nBasisFcts+1:nBasisFcts+nBaselineComps);
        
        % Save modeling results
        parsOut.ph0 = ph0;
        parsOut.ph1 = 0;
        parsOut.gaussLB = gaussLB;
        parsOut.lorentzLB = lorentzLB;
        parsOut.freqShift = freqShift;
        parsOut.metAmpl = metAmpl;
        parsOut.baseAmpl = baseAmpl;
        obj.Model{step}.fit.fit = prediction;
        obj.Model{step}.fit.baseline = baselineBasis*baseAmpl;
        obj.Model{step}.fit.residual = data-prediction;
        obj.Model{step}.fit.metabs   = fullBasis(:,1:nBasisFcts).*metAmpl';
        obj.Model{step}.time = time;
        obj.Model{step}.parsOut = parsOut;
        obj.step =step;
        
    end