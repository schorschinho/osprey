function fh = BasicPhysicsModel
    fh.lossFunction = @lossFunction;
    fh.forwardGradient = @forwardGradient;
    fh.forwardJacobian = @forwardJacobian;
    fh.forwardModel = @forwardModel;
    fh.x2pars = @x2pars;
    fh.pars2x = @pars2x;
    fh.ppmToIndex = @ppmToIndex;
    fh.fminunc_wrapper = @fminunc_wrapper;
end

% Forward models, loss functions, Jacobian and gradient functions
function sse = lossFunction(x, data, basisSet, baselineBasis, ppm, t, fitRange, SignalPart, Domain, SSE)
    
    if strcmp(Domain,'FD')   
        [indMin, indMax] = ppmToIndex(ppm, fitRange);
    end
    
    if strcmp(Domain,'FD') 
        data        = fftshift(fft(data, [], 1),1);
    end

    prediction  = forwardModel(x, basisSet, baselineBasis, ppm, t);
    diffVec     = data(indMin:indMax) - prediction(indMin:indMax);

    switch SignalPart
        case 'R'
            residual = real(diffVec);
        case 'I'
            residual = imag(diffVec);
        case {'RI', 'IR'}
            residual = cat(1, real(diffVec), imag(diffVec));
        case 'A'
            residual = abs(diffVec);
    end
    
    switch SSE
        case 'res'
            sse = residual;
        case 'sos'
            sse = sum(residual.^2);
    end
    
end

function [grad] = forwardGradient(x, data, basisSet, baselineBasis, ppm, t, fitRange, SignalPart)
            
    [indMin, indMax] = ppmToIndex(ppm, fitRange);
    
    prediction  = forwardModel(x, basisSet, baselineBasis, ppm, t);
    
    % Construct the Jacobian matrix of partial derivatives
    fids = basisSet.fids;
    nBasisFcts = size(fids,2);
    nBaselineComps = size(baselineBasis, 2);
    
    inputParams = x2pars(x, nBasisFcts);
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

    
    Fmett = fftshift(fft(-T_t .* Fmet, [], 1),1);
    Fmett2gauss = fftshift(fft(-2 .* gaussLB .* T_tt .* Fmet, [], 1),1);
    
    % partial derivatives
    dYdmetAmpl      = T_ph_basis .* fftshift(fft(Fmet, [], 1),1);
    dYdfreqShift    = T_ph_basis .* (-1j) .* Fmett .* metAmpl';
    dYdlorentzLB    = T_ph_basis  .* Fmett .* metAmpl';
    dYdgaussLB      = T_ph .* Fmett2gauss * metAmpl;
    dYdph0          = (-1j) .* prediction;
    dYdph1          = (-1j) .* ppm' .* prediction;
    if nBaselineComps ~= 0
        dYdBj           = T_ph_baseline .* baselineBasis;
    else
        dYdBj = [];
    end
                
    % reduce to fit range
    dYdph0          = dYdph0(indMin:indMax,:);
    dYdph1          = dYdph1(indMin:indMax,:);
    dYdgaussLB      = dYdgaussLB(indMin:indMax,:);
    dYdlorentzLB    = dYdlorentzLB(indMin:indMax,:);
    dYdfreqShift    = dYdfreqShift(indMin:indMax,:);
    dYdmetAmpl      = dYdmetAmpl(indMin:indMax,:);
    if nBaselineComps ~= 0
        dYdBj           = dYdBj(indMin:indMax,:);
    end
    
    % reduce prediction to fit range;
    prediction = prediction(indMin:indMax);
    data = fftshift(fft(data,[],1),1);
    data = data(indMin:indMax);
    
    % fft data and reduce to fit range
    
    jac = cat(2, dYdph0, dYdph1, dYdgaussLB, dYdlorentzLB, dYdfreqShift, dYdmetAmpl, dYdBj);
    
    if strcmp(SignalPart,'R') 
        data        = real(data);
        prediction  = real(prediction);
        jac         = real(jac);
    end
    if strcmp(SignalPart,'I') 
        data        = imag(data);
        prediction  = imag(prediction);
        jac         = imag(jac);
    end
    if strcmp(SignalPart,'RI') 
        data        = cat(1, real(data), imag(data));
        prediction  = cat(1, real(prediction), imag(prediction));
        jac         = cat(1, real(jac), imag(jac));
    end
    if strcmp(SignalPart,'A') 
        data        = abs(data);
        prediction  = abs(prediction);
        jac         = abs(jac);
    end

    grad = sum((data - prediction).*(-conj(jac)) + (-jac .* conj(data-prediction)));
    grad = grad';

end

function [jac] = forwardJacobian(x, data, basisSet, baselineBasis, ppm, t, fitRange,SignalPart)
            
    [indMin, indMax] = ppmToIndex(ppm, fitRange);
    
    prediction  = forwardModel(x, basisSet, baselineBasis, ppm, t);
    
    % Construct the Jacobian matrix of partial derivatives
    fids = basisSet.fids;
    nBasisFcts = size(fids,2);
    nBaselineComps = size(baselineBasis, 2);
    
    inputParams = x2pars(x, nBasisFcts);
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

    
    Fmett = fftshift(fft(-T_t .* Fmet, [], 1),1);
    Fmett2gauss = fftshift(fft(-2 .* gaussLB .* T_tt .* Fmet, [], 1),1);
    
    % partial derivatives
    dYdmetAmpl      = T_ph_basis .* fftshift(fft(Fmet, [], 1),1);
    dYdfreqShift    = T_ph_basis .* (-1j) .* Fmett .* metAmpl';
    dYdlorentzLB    = T_ph_basis  .* Fmett .* metAmpl';
    dYdgaussLB      = T_ph .* Fmett2gauss * metAmpl;
    dYdph0          = (-1j) .* prediction;
    dYdph1          = (-1j) .* ppm' .* prediction;
    if nBaselineComps ~= 0
        dYdBj           = T_ph_baseline .* baselineBasis;
    else
        dYdBj = [];
    end
                
    % reduce to fit range
    dYdph0          = dYdph0(indMin:indMax,:);
    dYdph1          = dYdph1(indMin:indMax,:);
    dYdgaussLB      = dYdgaussLB(indMin:indMax,:);
    dYdlorentzLB    = dYdlorentzLB(indMin:indMax,:);
    dYdfreqShift    = dYdfreqShift(indMin:indMax,:);
    dYdmetAmpl      = dYdmetAmpl(indMin:indMax,:);
    if nBaselineComps ~= 0
        dYdBj           = dYdBj(indMin:indMax,:);
    end
    
    % reduce prediction to fit range;
    prediction = prediction(indMin:indMax);
    data = fftshift(fft(data,[],1),1);
    data = data(indMin:indMax);
    
    % fft data and reduce to fit range
    
    jac = cat(2, dYdph0, dYdph1, dYdgaussLB, dYdlorentzLB, dYdfreqShift, dYdmetAmpl, dYdBj);

    if strcmp(SignalPart,'R') 
        jac         = real(jac);
    end
    if strcmp(SignalPart,'I') 
        jac         = imag(jac);
    end
    if strcmp(SignalPart,'RI') 
        jac = cat(1, real(jac), imag(jac));
    end
    if strcmp(SignalPart,'A') 
        jac         = abs(jac);
    end
    if strcmp(SignalPart,'C')
        % just return the complex
    end
    
    jac = (-1) * jac;
end

function [Y, baseline, metabs] = forwardModel(x, basisSet, baselineBasis, ppm, t)
            
    % Define the default 1-D forward model first
    fids = basisSet.fids;
    nBasisFcts = sum(basisSet.includeInFit(end,:));
    
    inputParams     = x2pars(x, nBasisFcts);
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
    
    T_ph = exp(-1j .* (ph0 + ph1.*ppm)');
    mets = specs * metAmpl;
    baseline = baselineBasis * baseAmpl;

    % Return full model, metabolites, and baseline for plotting
    if ~isempty(baseline)
        Y = T_ph .* (mets + baseline);
        baseline = T_ph .* baseline;
    else
        Y = T_ph .* mets;
        baseline = zeros(size(specs));
    end
    metabs = repmat(T_ph, [1 size(specs,2)]) .* specs .* repmat(metAmpl', [size(specs,1) 1]);
    
end

function paramStruct = x2pars(x, nBasisFcts)
            
    % Converts a 1-D x vector into a parameter struct
    paramStruct.ph0 = x(1);
    paramStruct.ph1 = x(2);
    paramStruct.gaussLB = x(3);
    paramStruct.lorentzLB = x(4:3+nBasisFcts);
    paramStruct.freqShift = x(4+nBasisFcts:3+2*nBasisFcts);
    paramStruct.metAmpl = x(4+2*nBasisFcts:3+3*nBasisFcts);
    if ~(4+3*nBasisFcts > length(x))
        paramStruct.baseAmpl = x(4+3*nBasisFcts:end);
    else
        paramStruct.baseAmpl = [];
    end

    
end

function x = pars2x(paramStruct)
            
    % Converts a parameter struct into a 1-D x vector that can be
    % passed on to solvers
    if ~isempty(paramStruct.baseAmpl)
        x = [paramStruct.ph0, paramStruct.ph1, paramStruct.gaussLB, ...
             paramStruct.lorentzLB, paramStruct.freqShift, ... 
             paramStruct.metAmpl, paramStruct.baseAmpl]';
    else
        x = [paramStruct.ph0, paramStruct.ph1, paramStruct.gaussLB, ...
             paramStruct.lorentzLB, paramStruct.freqShift, ... 
             paramStruct.metAmpl]';
    end
    
end

function [f,g,h] = fminunc_wrapper(x,F,GJ,H)
    % [f,g,h] = fminunc_wrapper( x, F, GJ, H )
    % for use with Matlab's "fminunc"
    f = F(x);
    if nargin > 2 && nargout > 1
        g = GJ(x);
    end
    if nargin > 3 && nargout > 2
        h = H(x);
    end
end