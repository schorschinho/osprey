function fh = GeneralizedBasicPhysicsModel
    fh.lossFunction = @lossFunction;
    fh.forwardGradient = @forwardGradient;
    fh.forwardJacobian = @forwardJacobian;
    fh.updateJacobianBlock = @updateJacobianBlock;
    fh.forwardModel = @forwardModel;
    fh.x2pars = @x2pars;
    fh.pars2x = @pars2x;
    fh.ppmToIndex = @ppmToIndex;
    fh.fminunc_wrapper = @fminunc_wrapper;
end

% Forward models, loss functions, Jacobian and gradient functions
function sse = lossFunction(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange, SignalPart, Domain, SSE, parametrizations)
    
    if strcmp(Domain,'FD')   
        [indMin, indMax] = ppmToIndex(ppm, fitRange);
    end
    
    if strcmp(Domain,'FD') 
        data        = fftshift(fft(data, [], 1),1);
    end

    prediction  = forwardModel(x, basisSet, baselineBasis, ppm, t, parametrizations);
    diffVec     = data(indMin:indMax,:) - prediction(indMin:indMax,:);

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
    
    % Reshape the multidimenisonal case
    residual = reshape(residual,[],1);

    switch SSE
        case 'res'
            sse = residual;
        case 'sos'
            sse = sum(residual.^2);
    end

    % Not sure if this is working needs furhter testing
%     SigmaSquared = (NormNoise * length(data(indMin:indMax)))^2;
%     sse = sse / SigmaSquared;
    
end

% Create the forward gradient (2D case is not implemented yet)
function [grad] = forwardGradient(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange, SignalPart, parametrizations)
            
    [indMin, indMax] = ppmToIndex(ppm, fitRange);
    
    prediction  = forwardModel(x, basisSet, baselineBasis, ppm, t, parametrizations);
    
    % Construct the Jacobian matrix of partial derivatives
    fids = basisSet.fids;
    nBasisFcts = size(fids,2);
    nBaselineComps = size(baselineBasis, 2);
    secDim = size(data,2);
    
    inputParams = x2pars(x, secDim, parametrizations);
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
        dYdbaseAmpl           = T_ph_baseline .* baselineBasis;
    else
        dYdbaseAmpl = [];
    end
                
    % reduce to fit range
    dYdph0          = dYdph0(indMin:indMax,:);
    dYdph1          = dYdph1(indMin:indMax,:);
    dYdgaussLB      = dYdgaussLB(indMin:indMax,:);
    dYdlorentzLB    = dYdlorentzLB(indMin:indMax,:);
    dYdfreqShift    = dYdfreqShift(indMin:indMax,:);
    dYdmetAmpl      = dYdmetAmpl(indMin:indMax,:);
    if nBaselineComps ~= 0
        dYdbaseAmpl           = dYdbaseAmpl(indMin:indMax,:);
    end
    
    % reduce prediction to fit range;
    prediction = prediction(indMin:indMax);
    data = fftshift(fft(data,[],1),1);
    data = data(indMin:indMax);
    
    % fft data and reduce to fit range
    
    jac = cat(2, dYdph0, dYdph1, dYdgaussLB, dYdlorentzLB, dYdfreqShift, dYdmetAmpl, dYdbaseAmpl);
    
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

% Not sure if this is working needs furhter testing
%     SigmaSquared = (NormNoise * length(data))^2;
%     grad = grad / SigmaSquared;

end

% Create forward jacobian (2D case seems to work for a single basis
% function)
function [jac] = forwardJacobian(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange,SignalPart, parametrizations)
       
    % Initialize output and partial derivatives to allow loop across
    % indirect dimension
    dYdph0          = [];
    dYdph1          = [];
    dYdgaussLB      = [];
    dYdlorentzLB    = [];
    dYdfreqShift    = [];
    dYdmetAmpl      = [];
    dYdbaseAmpl     = [];

    [indMin, indMax] = ppmToIndex(ppm, fitRange);
    prediction  = forwardModel(x, basisSet, baselineBasis, ppm, t, parametrizations);
    
    
    fidsBasis = basisSet.fids;
    nBasisFcts = size(fidsBasis,2);
    nBaselineComps = size(baselineBasis, 2);
    secDim = size(data,2);
    inputParams = x2pars(x, secDim, parametrizations);

    % Construct the Jacobian matrix of partial derivatives. This is done in
    % a block-wise fashion for 2-D data. 

    % Loop over indirect dimension
    for sD = 1 : secDim
        
        fids        = squeeze(fidsBasis(:,:,sD));
        gaussLB = squeeze(inputParams.gaussLB(sD,:));
        lorentzLB = squeeze(inputParams.lorentzLB(sD,:));
        freqShift = squeeze(inputParams.freqShift(sD,:));
        metAmpl = squeeze(inputParams.metAmpl(sD,:))';
        baseAmpl = squeeze(inputParams.baseAmpl(sD,:))';
        ph0 = squeeze(inputParams.ph0(sD));
        ph1 = squeeze(inputParams.ph1(sD));
        
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
        dYdmetAmpl      = cat(3,dYdmetAmpl,T_ph_basis .* fftshift(fft(Fmet, [], 1),1));
        dYdfreqShift    = cat(3,dYdfreqShift,T_ph_basis .* (-1j) .* Fmett .* metAmpl');
        dYdlorentzLB    = cat(3,dYdlorentzLB,T_ph_basis  .* Fmett .* metAmpl');
        dYdgaussLB      = cat(3,dYdgaussLB,T_ph .* Fmett2gauss * metAmpl);
        dYdph0          = cat(3,dYdph0,(-1j) .* prediction(:,sD));
        dYdph1          = cat(3,dYdph1,(-1j) .* ppm' .* prediction(:,sD));
        if nBaselineComps ~= 0
            dYdbaseAmpl           = cat(3,dYdbaseAmpl,T_ph_baseline .* baselineBasis);
        else
            dYdbaseAmpl = cat(3,dYdbaseAmpl,[]);
        end
                           
    end
    
    % reduce to fit range
    dYdph0          = dYdph0(indMin:indMax,:,:);
    dYdph1          = dYdph1(indMin:indMax,:,:);
    dYdgaussLB      = dYdgaussLB(indMin:indMax,:,:);
    dYdlorentzLB    = dYdlorentzLB(indMin:indMax,:,:,:);
    dYdfreqShift    = dYdfreqShift(indMin:indMax,:,:,:);
    dYdmetAmpl      = dYdmetAmpl(indMin:indMax,:,:,:);
    if nBaselineComps ~= 0
        dYdbaseAmpl           = dYdbaseAmpl(indMin:indMax,:,:,:);
    end

    % update each block in the jacobian according to the 2-D parametrization
    if secDim > 1
        [dYdph0]        = updateJacobianBlock(dYdph0,'ph0', parametrizations,inputParams);
        [dYdph1]        = updateJacobianBlock(dYdph1,'ph1', parametrizations,inputParams);
        [dYdgaussLB]    = updateJacobianBlock(dYdgaussLB,'gaussLB', parametrizations,inputParams);
        [dYdlorentzLB]  = updateJacobianBlock(dYdlorentzLB,'lorentzLB', parametrizations,inputParams);
        [dYdfreqShift]  = updateJacobianBlock(dYdfreqShift,'freqShift', parametrizations,inputParams);
        [dYdmetAmpl]    = updateJacobianBlock(dYdmetAmpl,'metAmpl', parametrizations,inputParams);
        [dYdbaseAmpl]   = updateJacobianBlock(dYdbaseAmpl,'baseAmpl', parametrizations,inputParams);        
    end

    jac = cat(2, dYdph0, dYdph1, dYdgaussLB, dYdlorentzLB, dYdfreqShift, dYdmetAmpl, dYdbaseAmpl);

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

% Not sure if this is working needs furhter testing    
%     SigmaSquared = (NormNoise * size(jac,1))^2;
%     jac = jac / SigmaSquared;
end

% Update jacobian according to the parametrization
function [dYdX] = updateJacobianBlock(dYdX,parameterName, parametrizations,inputParams)
    % Get dimensions
    nPoints = size(dYdX,1);
    nLines = size(dYdX,2);
    secDim = size(dYdX,3);
    % Free parametrizations need to add secDim copies to the jacobian and
    % set partial derivatives to zero for e.g. df1/dph02 .
    if strcmp(parametrizations.(parameterName).type,'free')
        dYdX = squeeze(dYdX);                                               %Remove length 1 dims (e.g. for ph0, ph1, gaussLB)
        dYdX = repmat(dYdX,[1 1 1 secDim]);                                 %Create secDim copies 
        dYdX = squeeze(dYdX);                                               %Remove length 1 dims (needed for 3D dYdX case)
        factor = repmat(eye(secDim),[1 1 nPoints nLines]);                  %To delete partial derivatives not on diagonal
        factor = permute(factor,[3 4 1 2]);                                 %Reorder to match dYdX dimensions
        factor = squeeze(factor);                                           %Remove length 1 dims
        dYdX = dYdX .* factor;                                              %Delete partial derivatives not on diagonal
    end
    % Fixed parametrization needs to concatenate along secDim resulting in
    % a single line in the jacobian
    if strcmp(parametrizations.(parameterName).type,'fixed')
        dYdX = squeeze(dYdX);                                              %Remove length 1 dims (e.g. for ph0, ph1, gaussLB)
    end
    % Dynamic parametrization needs be updated according to external
    % function and has to include modified lines in te jacobian
    if strcmp(parametrizations.(parameterName).type,'dynamic')
        parameterEstimate = [];
        % Loop over new parameters and get estimates
        for rp = 1 : length(parametrizations.(parameterName).parameterNames)    
            parameterEstimate = cat(1,parameterEstimate,inputParams.([parameterName 'Reparametrization']).(parametrizations.(parameterName).parameterNames{rp}));
        end
        % Calculate the jacobian according to the external function, parameter estimates, and modulator 
        factor = parametrizations.metAmpl.fun.jac(parameterEstimate,parametrizations.(parameterName).modulator);
        factor = repmat(factor, [1 1 1 nPoints]);                          % Repeat nPoints times
        factor = permute(factor,[4 2 1 3]);                                % Dims have to be nPoints secDim nLines nPars  
        dYdXOrginal = dYdX;                                                % Backup original derivatives
        % Multiply the original derivatives with the derivatives from the
        % reparametrization
        dYdX = [];
        % Loop over new parameters
        for rp = 1 : length(parametrizations.(parameterName).parameterNames)
            dYdX = cat(4,dYdX,dYdXOrginal .* factor(:,:,:,rp));
        end
        dYdX = squeeze(dYdX);                                               %Remove length 1 dims
        if ndims(dYdX) ==3                                                  % Dims have to be nPoints secDim nLines*nPars 
            dYdX = permute(dYdX,[1 3 2]);
            secDim = size(dYdX,3);
        else
            secDim = size(dYdX,4);                                          %THIS MIGHT BE THE PROBLEM FOR CASE nBasiss > 1
        end
        
    end
    % Finally we have to concatenate along the indirect dimension
    switch ndims(dYdX)
        case 2   % For fixed parametrizations of ph0, ph1, gaussLB                                   
            dYdX = reshape(dYdX,[],1);
        case 3   % E.g. linked metAmpls or single basis function case
            if strcmp(parametrizations.(parameterName).type,'fixed')
                dYdX = permute(dYdX,[1 3 2]);
            end
            if strcmp(parameterName,'baseAmpl') || strcmp(parameterName,'metAmpl') || strcmp(parameterName,'freqShift') || strcmp(parameterName,'lorentzLB')
                secDim = nLines;
            end
            dYdX = reshape(dYdX,[],secDim);
        case  4 % E.g. free metAmpls 
            dYdX = permute(dYdX,[1 3 4 2]);
            dYdX = reshape(dYdX,[],secDim*nLines);
    end
end

% Create forward model
function [Y, baseline, metabs] = forwardModel(x, basisSet, baselineBasis, ppm, t, parametrizations)
    
    % Initialize output
    Y = [];
    baseline = [];
    metabs = [];
    
    fidsBasis = basisSet.fids;
    nBasisFcts = sum(basisSet.includeInFit(end,:));
    secDim = size(fidsBasis,3);
    
    inputParams = x2pars(x, secDim, parametrizations);

    % Loop over indirect dimension
    for sD = 1 : secDim
        fids        = squeeze(fidsBasis(:,:,sD));
        gaussLB     = squeeze(inputParams.gaussLB(sD,:));
        lorentzLB   = squeeze(inputParams.lorentzLB(sD,:));
        freqShift   = squeeze(inputParams.freqShift(sD,:));
        metAmpl     = squeeze(inputParams.metAmpl(sD,:))';
        baseAmpl    = squeeze(inputParams.baseAmpl(sD,:))';
        ph0         = squeeze(inputParams.ph0(sD));
        ph1         = squeeze(inputParams.ph1(sD));
    
        timeDomainMultiplier = zeros(size(fids));
        for ll = 1:nBasisFcts
            timeDomainMultiplier(:,ll) = exp(-(1i*freqShift(ll) + lorentzLB(ll) + gaussLB.^2.*t).*t)';    
        end
        
        Fl = timeDomainMultiplier .* fids;
        specs = fftshift(fft(Fl, [], 1),1);
        
        T_ph = exp(-1j .* (ph0 + ph1.*ppm)');
        mets = specs * metAmpl;
        bl = baselineBasis * baseAmpl;
    
        % Return full model, metabolites, and baseline for plotting
        if ~isempty(bl)
            Y = cat(2,Y,T_ph .* (mets + bl));
            baseline = cat(2,baseline,T_ph .* bl);
        else
            Y = cat(2,Y,T_ph .* mets);
            baseline = cat(2,baseline,zeros(size(specs)));
        end
        metabs = cat(3,metabs,repmat(T_ph, [1 size(specs,2)]) .* specs .* repmat(metAmpl', [size(specs,1) 1]));
    end   
end

function paramStruct = x2pars(x, secDim, parametrizations)      
    pars = fields(parametrizations);        
    % Converts a 1-D x vector into a parameter struct
    for ff = 1 : length(pars)
        % Start with a parameter struct of 1-D vectors according to the
        % indices
        paramStruct.(pars{ff}) = x(parametrizations.(pars{ff}).start:parametrizations.(pars{ff}).end);  

        % Reshape the 1-D vectors according to the number of basis functions,
        % number of subspectra, and parametrization
        switch pars{ff}
            case {'ph0','ph1','gaussLB'}
                if strcmp(parametrizations.(pars{ff}).type,'free')
                    paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),secDim,1);
                end
                if strcmp(parametrizations.(pars{ff}).type,'fixed')
                    paramStruct.(pars{ff}) = repmat(paramStruct.(pars{ff}),[secDim,1]);
                end
                if strcmp(parametrizations.(pars{ff}).type,'dynamic')
                    paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),size(parametrizations.(pars{ff}).lb));
                    paramStruct.(pars{ff}) = parametrizations.metAmpl.fun.fun(paramStruct.(pars{ff}),parametrizations.(pars{ff}).modulator);
                end
            case {'metAmpl', 'freqShift', 'lorentzLB','baseAmpl'}
                if strcmp(parametrizations.(pars{ff}).type,'free')
                    paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),secDim,[]);
                end
                if strcmp(parametrizations.(pars{ff}).type,'fixed')
                    paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),1,[]);
                    paramStruct.(pars{ff}) = repmat(paramStruct.(pars{ff}),[secDim,1]);
                end
                if strcmp(parametrizations.(pars{ff}).type,'dynamic')
                    paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),size(parametrizations.(pars{ff}).lb));
                    for rp = 1 : length(parametrizations.(pars{ff}).parameterNames)
                        paramStruct.([pars{ff} 'Reparametrization']).(parametrizations.(pars{ff}).parameterNames{rp}) = paramStruct.(pars{ff})(rp,:);
                    end
                    paramStruct.(pars{ff}) = parametrizations.(pars{ff}).fun.fun(paramStruct.(pars{ff}),parametrizations.(pars{ff}).modulator);
                end
        end
    end
   


end

function [x,indexStruct] = pars2x(paramStruct)
    % Converts a parameter struct into a 1-D x vector that can be
    % passed on to solvers. It also defines start and end indices for
    % easier identification
    pars = fields(paramStruct);
    x = [];
    for ff = 1 : length(pars)
        if ismember(pars{ff},{'ph0','ph1','gaussLB','lorentzLB','freqShift','metAmpl','baseAmpl'})
            % Skip new parameters from dynamic parameterization bc they
            % shoudl not turn up in the x vector
            if ~isempty(paramStruct.(pars{ff}))             
                if isempty(x)
                    indexStruct.(pars{ff}).start = 1;
                else
                    indexStruct.(pars{ff}).start = length(x)+1;
                end
                x = cat(2,x,reshape(paramStruct.(pars{ff}),1,[]));
                indexStruct.(pars{ff}).end = length(x);
            end
        end
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