%%  GeneralizedBasicPhysicsModel.m
%   This function contains the generalized physics model allowing for 2D
%   modeling of arbitray MRS data in different domains. It is desinged with
%   the highest possible felxibility.
%   There are additional functions that are not included as handles. They
%   are requried to construct the jacobian for 2D modeling and perform
%   parameter regularization:
%   
%   lossFunction    - calculates loss function for optimizer
%   forwardGradient - calculates the forward gradient for optimizer
%   forwardJacobian - calculates the forward jacobian for optimizer
%   forwardModel    - converts parameter struct into model prediction
%   x2pars          - converts x vector to parameter struct
%   pars2x          - converts parameter struct to x vector 
%   fminunc_wrapper - wrapper for MATLAB's fminunc
% 
%   It integreates into the OspreyFitObj environment and allows for easy
%   changes of solvers. You can also use this function as a template for
%   your own Physics model. For full functionalty you need to define the
%   same functions and handels described in the first section. This
%   includes the following 7 functions:
%   updateJacobianBlock - updates 2D jacobian blocks (only needed for 2D)
%   addParameterRegularization - adds regularization (only used for reg)
%   
%   USAGE:
%       Specify this in the module.ModelFunction
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Handle set up for optimizer
function fh = GeneralizedBasicPhysicsModel
    fh.lossFunction = @lossFunction;
    fh.forwardGradient = @forwardGradient;
    fh.forwardJacobian = @forwardJacobian;
    fh.forwardModel = @forwardModel;
    fh.x2pars = @x2pars;
    fh.pars2x = @pars2x;
    fh.fminunc_wrapper = @fminunc_wrapper;
end

%% Define lossfunction, forward gradient, forward jacobian, forward model, x2pars, pars2x, fminunc_wrapper

function sse = lossFunction(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange, SignalPart, Domain, SSE,Reg, parametrizations)
% This function generates the output for the solver accoring to the
% lossfunction. This includes all settings described in the model procedure
% json.
%
%   USAGE:
%       sse = lossFunction(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange, SignalPart, Domain, SSE,Reg, parametrizations)
%
%   INPUTS:
%       x                = x vector with parameters to optimize
%       data             = ppm axis
%       NormNoise        = 1-point noise estimate
%       basisSet         = basis set struct
%       baselineBasis    = baseline basis set
%       ppm              = ppm axis
%       t                = time vector
%       fitRange         = model range
%       SignalPart       = optimization signal part
%       Domain           = optimization domain
%       SSE              = lossfunction string
%       Reg              = regularizer flag
%       parametrizations = parameter struct
%
%   OUTPUTS:
%       sse       = return for solver (sse or residual vector)
%
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Calculate lossfunction
    % Run forward model to get prediction and regularizer output
    [prediction, ~, ~, regu]  = forwardModel(x, ...                         % x vector with parameters to optimize
                                             basisSet, ...                  % basis set struct
                                             baselineBasis, ...             % baseline basis set
                                             ppm, ...                       % ppm axis
                                             t, ...                         % time vector
                                             Reg, ...                       % regularizer flag
                                             parametrizations);             % parameter struct
    
    
    if strcmp(Domain,'FD')                                                  % frequency domain optimization  
        [indMin, indMax] = ppmToIndex(ppm, fitRange);                       % Get fit range ppm indices
        data        = fftshift(fft(data, [], 1),1);                         % Convert data frequency domain data
        data        = data(indMin:indMax,:);                                % Cut out fit range from data
        prediction  = prediction(indMin:indMax,:);                          % Cut out fit range from prediction
    end    
    if strcmp(Domain,'TD')                                                  % time domain optimization         
        prediction  = ifft(ifftshift(prediction,1), [], 1);                 % Convert prediction to time domain
        prediction  = ifft(ifftshift(regu,1), [], 1);                       % Convert regularizer to time domain
    end
    residual     = data - prediction;                                       % Calculate residual

    switch SignalPart                                                       % Switch for optimization signal part
        case 'R'
            residual = real(residual);                                      % Take real part of residual
            regu     = real(regu);                                          % Take real part of regularizer
        case 'I'
            residual = imag(residual);                                      % Take imaginary part
            regu     = imag(regu);                                          % Take imaginary part of regularizer
        case {'RI', 'IR'}
            residual = cat(1, real(residual), imag(residual));              % Concatenate real and imaginary part
            regu     = cat(1, real(regu), imag(regu));                      % Concatenate real and imaginary part
        case 'A'
            residual = abs(residual);                                       % Take magnitude residual
            regu     = abs(regu);                                           % Take magnitude regularizer
    end
    
    if ~isempty(NormNoise)                                                  % Don't normalize residual if GradientCheck is performed               
        SigmaSquared = (NormNoise * sqrt(size(data,1))).^2;                 % Normalization factor for residual
        SigmaSquared = repmat(SigmaSquared, [size(data,1) 1]);              % Repeat according to dimensions
        residual = residual ./ SigmaSquared;                                % Normalize residual
        SigmaSquared = (NormNoise * sqrt(size(data,1) + size(regu,1))).^2;  % Normalization factor for regularizer
        if size(regu,1) > 0                                                 % Has regularizer
            SigmaSquared = repmat(SigmaSquared, [size(regu,1) 1]);          % Repeat according to dimensions
            regu = regu ./ SigmaSquared;                                    % Normalize regularizer
        end
    end

    residual  = cat(1,residual,-regu);                                      % Add regularizer term
    residual = reshape(residual,[],1);                                      % Reshape the multidimenisonal case

    switch SSE                                                              % Switch for return to solver
        case 'res'
            sse = residual;                                                 % Return residual
        case 'sos'
            sse = sum(residual.^2);                                         % Return sum of squares
    end    
end

function grad = forwardGradient(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange, SignalPart,Reg, parametrizations)
% This function generates the gradient for the generalized physics model
%
%   USAGE:
%       grad = forwardGradient(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange,SignalPart,Reg, parametrizations)
%
%   INPUTS:
%       x                = x vector with parameters to optimize
%       data             = ppm axis
%       NormNoise        = 1-point noise estimate
%       basisSet         = basis set struct
%       baselineBasis    = baseline basis set
%       ppm              = ppm axis
%       t                = time vector
%       fitRange         = model range
%       SignalPart       = optimization signal part
%       Reg              = regularizer flag
%       parametrizations = parameter struct
%
%   OUTPUTS:
%       grad             = return for gradient matrix
%
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Calculate forward gradient  
    % Partial derivatives to allow loop across
    dYdph0          = [];                                                   % Initialize partial derivative wrt ph0
    dYdph1          = [];                                                   % Initialize partial derivative wrt ph1
    dYdgaussLB      = [];                                                   % Initialize partial derivative wrt gaussLB
    dYdlorentzLB    = [];                                                   % Initialize partial derivative wrt lorentzLB
    dYdfreqShift    = [];                                                   % Initialize partial derivative wrt freqShift
    dYdmetAmpl      = [];                                                   % Initialize partial derivative wrt metAmpl
    dYdbaseAmpl     = [];                                                   % Initialize partial derivative wrt baseAmpl

    fidsBasis = basisSet.fids;                                              % Get basis function
    nBasisFcts = size(fidsBasis,2);                                         % Get number basis functions
    nBaselineComps = size(baselineBasis, 2);                                % Get number baseline basis function                           
    secDim = size(data,2);                                                  % Get number spectra along indirect dimension

    inputParams = x2pars(x, secDim, parametrizations);                      % Convert x vector to parameter struct

    [indMin, indMax] = ppmToIndex(ppm, fitRange);                           % Get fit range ppm indices

    % Run forward model to get prediction and regularizer output
    [prediction,~,~,regu]  = forwardModel(x, ...                            % x vector with parameters to optimize
                                          basisSet, ...                     % basis set struct
                                          baselineBasis, ...                % baseline basis set
                                          ppm, ...                          % ppm axis
                                          t, ...                            % time vector
                                          Reg, ...                          % regularizer flag
                                          parametrizations);                % parameter struct


    

    % Construct the Jacobian matrix of partial derivatives. This is done in
    % a block-wise fashion for 2-D data.    
    for sD = 1 : secDim                                                     % Loop over indirect dimension
        
        fids        = squeeze(fidsBasis(:,:,sD));                           % Get time domain data
        gaussLB = squeeze(inputParams.gaussLB(sD,:));                       % Get gaussLB parameter
        lorentzLB = squeeze(inputParams.lorentzLB(sD,:));                   % Get lorentzLB parameter
        freqShift = squeeze(inputParams.freqShift(sD,:));                   % Get freqShift parameter
        metAmpl = squeeze(inputParams.metAmpl(sD,:))';                      % Get metAmpl parameter
        baseAmpl = squeeze(inputParams.baseAmpl(sD,:))';                    % Get baseAmpl parameter
        ph0 = squeeze(inputParams.ph0(sD));                                 % Get ph0 parameter
        ph1 = squeeze(inputParams.ph1(sD));                                 % Get ph1 parameter
        
        timeDomainMultiplier = zeros(size(fids));                           % Setup time domain multiplier
        for ll = 1:nBasisFcts                                               % Loop over basis functions
            timeDomainMultiplier(:,ll) = exp(-(1i*freqShift(ll) + lorentzLB(ll) + gaussLB.^2.*t).*t)';  % Apply freqshifts, lorentzianLB, and gaussLB  
        end                                                                 % End loop over basis functions
        
        T_ph = exp(-1j .* (ph0 + ph1.*ppm)');                               % Create phase evolution
        T_ph_basis = repmat(T_ph, [1, nBasisFcts]);                         % Repeat phase vector for all basis functions
        T_ph_baseline = repmat(T_ph, [1, nBaselineComps]);                  % Repeat phase vector for all baseline basis functions
        T_t = repmat(t', [1, nBasisFcts]);                                  % Create t multiplier
        T_tt = repmat(t'.*t', [1, nBasisFcts]);                             % Create t*t multiplier
        
        Fmet = timeDomainMultiplier .* fids;                                % Apply time domain multiplier                
        Fmett = fftshift(fft(-T_t .* Fmet, [], 1),1);                       % Apply t multiplier
        Fmett2gauss = fftshift(fft(-2 .* gaussLB .* T_tt .* Fmet, [], 1),1);% Apply t*t multiplier and gaussLB for partial derivaitve wrt gaussLB
        
        % Calculate and concatenate partial derivatives
        dYdmetAmpl      = cat(3,dYdmetAmpl,T_ph_basis .* fftshift(fft(Fmet, [], 1),1)); % Partial derivative wrt metAmpl
        dYdfreqShift    = cat(3,dYdfreqShift,T_ph_basis .* (-1j) .* Fmett .* metAmpl'); % Partial derivative wrt freqShift
        dYdlorentzLB    = cat(3,dYdlorentzLB,T_ph_basis  .* Fmett .* metAmpl');         % Partial derivative wrt lorentzLB
        dYdgaussLB      = cat(3,dYdgaussLB,T_ph .* Fmett2gauss * metAmpl);              % Partial derivative wrt gaussLB
        dYdph0          = cat(3,dYdph0,(-1j) .* prediction(:,sD));                      % Partial derivative wrt ph0
        dYdph1          = cat(3,dYdph1,(-1j) .* ppm' .* prediction(:,sD));              % Partial derivative wrt ph1
        if nBaselineComps ~= 0                                                          % Has baseline?
            dYdbaseAmpl           = cat(3,dYdbaseAmpl,T_ph_baseline .* baselineBasis);  % Partial derivative wrt baseAmpl
        else
            dYdbaseAmpl = cat(3,dYdbaseAmpl,[]);                                        % Empty partial derivative when no baseline is defined
        end
                           
    end                                                                     % End loop over indirect dimension
    
    % Reduce to fit range
    dYdph0          = dYdph0(indMin:indMax,:,:);                            % Cut out fit range
    dYdph1          = dYdph1(indMin:indMax,:,:);                            % Cut out fit range
    dYdgaussLB      = dYdgaussLB(indMin:indMax,:,:);                        % Cut out fit range
    dYdlorentzLB    = dYdlorentzLB(indMin:indMax,:,:,:);                    % Cut out fit range
    dYdfreqShift    = dYdfreqShift(indMin:indMax,:,:,:);                    % Cut out fit range
    dYdmetAmpl      = dYdmetAmpl(indMin:indMax,:,:,:);                      % Cut out fit range
    if nBaselineComps ~= 0                                                  % Has baseline?
        dYdbaseAmpl           = dYdbaseAmpl(indMin:indMax,:,:,:);           % Cut out fit range
    end
    
    
    if Reg                                                                  % Add parameter regularization
        for sD = 1 : secDim                                                 % Loop over indirect dimension 
            ph0 = squeeze(inputParams.ph0(sD));                             % Get ph0 parameter
            ph1 = squeeze(inputParams.ph1(sD));                             % Get ph1 parameter
            gaussLB = squeeze(inputParams.gaussLB(sD,:));                   % Get gaussLB parameter
            lorentzLB = squeeze(inputParams.lorentzLB(sD,:));               % Get lorentzLB parameter
            freqShift = squeeze(inputParams.freqShift(sD,:));               % Get freqShift parameter
            metAmpl = squeeze(inputParams.metAmpl(sD,:))';                  % Get metAmpl parameter
            baseAmpl = squeeze(inputParams.baseAmpl(sD,:))';                % Get baseAmpl parameter
            [dYdph0]        = addParameterRegularization(dYdph0,'ph0', parametrizations,ph0,1); % Add regularizer to ph0 parameter
            [dYdph1]        = addParameterRegularization(dYdph1,'ph1', parametrizations,ph1,1); % Add regularizer to ph1 parameter
            [dYdgaussLB]    = addParameterRegularization(dYdgaussLB,'gaussLB', parametrizations,gaussLB,1); % Add regularizer to gaussLB parameter
            [dYdlorentzLB]  = addParameterRegularization(dYdlorentzLB,'lorentzLB', parametrizations,lorentzLB,1); % Add regularizer to lorentzLB parameter
            [dYdfreqShift]  = addParameterRegularization(dYdfreqShift,'freqShift', parametrizations,freqShift,1); % Add regularizer to freqShift parameter
            [dYdmetAmpl]    = addParameterRegularization(dYdmetAmpl,'metAmpl', parametrizations,metAmpl,1); % Add regularizer to metAmpl parameter
            [dYdbaseAmpl]   = addParameterRegularization(dYdbaseAmpl,'baseAmpl', parametrizations,baseAmpl,1); % Add regularizer to baseAmpl parameter 
        end
    end                                                                     % End loop over indirect dimension


   
    if secDim > 1                                                           % update each block in the jacobian according to the 2-D parametrization
        dataPoints = size(data(indMin:indMax,:),1);                         % Get number of datapoints
        regularizerPoints = dYdph0 - dataPoints;                            % Get number of points in regularizer
        [dYdph0]        = updateJacobianBlock(dYdph0,'ph0', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for ph0 parameter
        [dYdph1]        = updateJacobianBlock(dYdph1,'ph1', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for ph1 parameter
        [dYdgaussLB]    = updateJacobianBlock(dYdgaussLB,'gaussLB', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for gaussLB parameter
        [dYdlorentzLB]  = updateJacobianBlock(dYdlorentzLB,'lorentzLB', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for lorentzLB parameter
        [dYdfreqShift]  = updateJacobianBlock(dYdfreqShift,'freqShift', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for freqShift parameter
        [dYdmetAmpl]    = updateJacobianBlock(dYdmetAmpl,'metAmpl', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for metAmpl parameter
        [dYdbaseAmpl]   = updateJacobianBlock(dYdbaseAmpl,'baseAmpl', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise);  % Update jacobian block for baseAmpl parameter       
    end

    if nBaselineComps ~= 0                                                  % Has baseline?
        jac = cat(2, dYdph0, dYdph1, dYdgaussLB, dYdlorentzLB, dYdfreqShift, dYdmetAmpl, dYdbaseAmpl); % Create final jacobian
    else
        jac = cat(2, dYdph0, dYdph1, dYdgaussLB, dYdlorentzLB, dYdfreqShift, dYdmetAmpl);   % Create final jacobian
    end
           
    residual     = data - prediction;                                       % Calculate residual
    residual = residual(indMin:indMax,:);                                   % Cut frequency range from residual
    dataPoints = size(data(indMin:indMax,:),1);                             % Get number of datapoints
    residual      = cat(1,residual,-regu);                                  % Concatenate regularizer
        
    if  ~isempty(NormNoise)                                                 % No normlaization for GradientCheck                                                      
        Sigma = (NormNoise * sqrt(dataPoints));                             % Calculate normalization factor for data
        residual(1:dataPoints,:) = residual(1:dataPoints,:) / Sigma;        % Normalize data part of residual
        if secDim == 1                                                      % 1D jacobian have not been normalized yet
            jac(1:dataPoints,:) = jac(1:dataPoints,:) / Sigma;              % Normalize data part jacobian
        end
        regularizerPoints = size(residual,1) - dataPoints;                  % Get number of points in regularizer
        Sigma = (NormNoise * sqrt(dataPoints + regularizerPoints));         % Calculate normalization factor for regularizer
        residual(dataPoints+1:end,:) = residual(dataPoints+1:end,:) / Sigma;% Normalize regularizer part
        if secDim == 1                                                      % 1D jacobian have not been normalized yet
            jac(dataPoints+1:end,:) = jac(dataPoints+1:end,:) / Sigma;          % Normalize regularizer part
        end
    end

    switch SignalPart
        case 'R'
            residual    = real(residual);                               % Take real part of residual
            jac         = real(jac);                                        % Take real part of jacobian
        case 'I'
            grad        = imag(grad);                                       % Take imaginary of residual
            jac         = imag(jac);                                        % Take imaginary of jacobian
        case {'RI', 'IR'} 
            grad        = cat(1, real(grad), imag(grad));                   % Concatenate real and imaginary part
            jac         = cat(1, real(jac), imag(jac));                     % Concatenate real and imaginary part
        case 'A'
            grad        = abs(grad);                                        % Take magnitude of residual
            jac         = abs(jac);                                         % Take magnitude of jacobian
        case 'C'
                                                                            % just return the complex
    end
    
    grad = sum((residual).*(-conj(jac)) + (-jac .* conj(residual)));    % Calculate gradient
    grad = grad';
end

function jac = forwardJacobian(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange,SignalPart,Reg, parametrizations)
% This function generates the jacobian for the generalized physics model
%
%   USAGE:
%       jac = forwardJacobian(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange,SignalPart,Reg, parametrizations)
%
%   INPUTS:
%       x                = x vector with parameters to optimize
%       data             = ppm axis
%       NormNoise        = 1-point noise estimate
%       basisSet         = basis set struct
%       baselineBasis    = baseline basis set
%       ppm              = ppm axis
%       t                = time vector
%       fitRange         = model range
%       SignalPart       = optimization signal part
%       Reg              = regularizer flag
%       parametrizations = parameter struct
%
%   OUTPUTS:
%       jac              = return for jacobian matrix
%
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Calculate forward jacobian matrix       
    % Partial derivatives to allow loop across
    dYdph0          = [];                                                   % Initialize partial derivative wrt ph0
    dYdph1          = [];                                                   % Initialize partial derivative wrt ph1
    dYdgaussLB      = [];                                                   % Initialize partial derivative wrt gaussLB
    dYdlorentzLB    = [];                                                   % Initialize partial derivative wrt lorentzLB
    dYdfreqShift    = [];                                                   % Initialize partial derivative wrt freqShift
    dYdmetAmpl      = [];                                                   % Initialize partial derivative wrt metAmpl
    dYdbaseAmpl     = [];                                                   % Initialize partial derivative wrt baseAmpl

    fidsBasis = basisSet.fids;                                              % Get basis function
    nBasisFcts = size(fidsBasis,2);                                         % Get number basis functions
    nBaselineComps = size(baselineBasis, 2);                                % Get number baseline basis function                           
    secDim = size(data,2);                                                  % Get number spectra along indirect dimension

    inputParams = x2pars(x, secDim, parametrizations);                      % Convert x vector to parameter struct

    [indMin, indMax] = ppmToIndex(ppm, fitRange);                           % Get fit range ppm indices

    % Run forward model to get prediction and regularizer output
    prediction  = forwardModel(x, ...                                       % x vector with parameters to optimize
                             basisSet, ...                                  % basis set struct
                             baselineBasis, ...                             % baseline basis set
                             ppm, ...                                       % ppm axis
                             t, ...                                         % time vector
                             Reg, ...                                       % regularizer flag
                             parametrizations);                             % parameter struct


    

    % Construct the Jacobian matrix of partial derivatives. This is done in
    % a block-wise fashion for 2-D data.    
    for sD = 1 : secDim                                                     % Loop over indirect dimension
        
        fids        = squeeze(fidsBasis(:,:,sD));                           % Get time domain data
        gaussLB = squeeze(inputParams.gaussLB(sD,:));                       % Get gaussLB parameter
        lorentzLB = squeeze(inputParams.lorentzLB(sD,:));                   % Get lorentzLB parameter
        freqShift = squeeze(inputParams.freqShift(sD,:));                   % Get freqShift parameter
        metAmpl = squeeze(inputParams.metAmpl(sD,:))';                      % Get metAmpl parameter
        baseAmpl = squeeze(inputParams.baseAmpl(sD,:))';                    % Get baseAmpl parameter
        ph0 = squeeze(inputParams.ph0(sD));                                 % Get ph0 parameter
        ph1 = squeeze(inputParams.ph1(sD));                                 % Get ph1 parameter
        
        timeDomainMultiplier = zeros(size(fids));                           % Setup time domain multiplier
        for ll = 1:nBasisFcts                                               % Loop over basis functions
            timeDomainMultiplier(:,ll) = exp(-(1i*freqShift(ll) + lorentzLB(ll) + gaussLB.^2.*t).*t)';  % Apply freqshifts, lorentzianLB, and gaussLB  
        end                                                                 % End loop over basis functions
        
        T_ph = exp(-1j .* (ph0 + ph1.*ppm)');                               % Create phase evolution
        T_ph_basis = repmat(T_ph, [1, nBasisFcts]);                         % Repeat phase vector for all basis functions
        T_ph_baseline = repmat(T_ph, [1, nBaselineComps]);                  % Repeat phase vector for all baseline basis functions
        T_t = repmat(t', [1, nBasisFcts]);                                  % Create t multiplier
        T_tt = repmat(t'.*t', [1, nBasisFcts]);                             % Create t*t multiplier
        
        Fmet = timeDomainMultiplier .* fids;                                % Apply time domain multiplier                
        Fmett = fftshift(fft(-T_t .* Fmet, [], 1),1);                       % Apply t multiplier
        Fmett2gauss = fftshift(fft(-2 .* gaussLB .* T_tt .* Fmet, [], 1),1);% Apply t*t multiplier and gaussLB for partial derivaitve wrt gaussLB
        
        % Calculate and concatenate partial derivatives
        dYdmetAmpl      = cat(3,dYdmetAmpl,T_ph_basis .* fftshift(fft(Fmet, [], 1),1)); % Partial derivative wrt metAmpl
        dYdfreqShift    = cat(3,dYdfreqShift,T_ph_basis .* (-1j) .* Fmett .* metAmpl'); % Partial derivative wrt freqShift
        dYdlorentzLB    = cat(3,dYdlorentzLB,T_ph_basis  .* Fmett .* metAmpl');         % Partial derivative wrt lorentzLB
        dYdgaussLB      = cat(3,dYdgaussLB,T_ph .* Fmett2gauss * metAmpl);              % Partial derivative wrt gaussLB
        dYdph0          = cat(3,dYdph0,(-1j) .* prediction(:,sD));                      % Partial derivative wrt ph0
        dYdph1          = cat(3,dYdph1,(-1j) .* ppm' .* prediction(:,sD));              % Partial derivative wrt ph1
        if nBaselineComps ~= 0                                                          % Has baseline?
            dYdbaseAmpl           = cat(3,dYdbaseAmpl,T_ph_baseline .* baselineBasis);  % Partial derivative wrt baseAmpl
        else
            dYdbaseAmpl = cat(3,dYdbaseAmpl,[]);                                        % Empty partial derivative when no baseline is defined
        end
                           
    end                                                                     % End loop over indirect dimension
    
    % Reduce to fit range
    dYdph0          = dYdph0(indMin:indMax,:,:);                            % Cut out fit range
    dYdph1          = dYdph1(indMin:indMax,:,:);                            % Cut out fit range
    dYdgaussLB      = dYdgaussLB(indMin:indMax,:,:);                        % Cut out fit range
    dYdlorentzLB    = dYdlorentzLB(indMin:indMax,:,:,:);                    % Cut out fit range
    dYdfreqShift    = dYdfreqShift(indMin:indMax,:,:,:);                    % Cut out fit range
    dYdmetAmpl      = dYdmetAmpl(indMin:indMax,:,:,:);                      % Cut out fit range
    if nBaselineComps ~= 0                                                  % Has baseline?
        dYdbaseAmpl           = dYdbaseAmpl(indMin:indMax,:,:,:);           % Cut out fit range
    end
    
    
    if Reg                                                                  % Add parameter regularization
        for sD = 1 : secDim                                                 % Loop over indirect dimension 
            ph0 = squeeze(inputParams.ph0(sD));                             % Get ph0 parameter
            ph1 = squeeze(inputParams.ph1(sD));                             % Get ph1 parameter
            gaussLB = squeeze(inputParams.gaussLB(sD,:));                   % Get gaussLB parameter
            lorentzLB = squeeze(inputParams.lorentzLB(sD,:));               % Get lorentzLB parameter
            freqShift = squeeze(inputParams.freqShift(sD,:));               % Get freqShift parameter
            metAmpl = squeeze(inputParams.metAmpl(sD,:))';                  % Get metAmpl parameter
            baseAmpl = squeeze(inputParams.baseAmpl(sD,:))';                % Get baseAmpl parameter
            [dYdph0]        = addParameterRegularization(dYdph0,'ph0', parametrizations,ph0,1); % Add regularizer to ph0 parameter
            [dYdph1]        = addParameterRegularization(dYdph1,'ph1', parametrizations,ph1,1); % Add regularizer to ph1 parameter
            [dYdgaussLB]    = addParameterRegularization(dYdgaussLB,'gaussLB', parametrizations,gaussLB,1); % Add regularizer to gaussLB parameter
            [dYdlorentzLB]  = addParameterRegularization(dYdlorentzLB,'lorentzLB', parametrizations,lorentzLB,1); % Add regularizer to lorentzLB parameter
            [dYdfreqShift]  = addParameterRegularization(dYdfreqShift,'freqShift', parametrizations,freqShift,1); % Add regularizer to freqShift parameter
            [dYdmetAmpl]    = addParameterRegularization(dYdmetAmpl,'metAmpl', parametrizations,metAmpl,1); % Add regularizer to metAmpl parameter
            [dYdbaseAmpl]   = addParameterRegularization(dYdbaseAmpl,'baseAmpl', parametrizations,baseAmpl,1); % Add regularizer to baseAmpl parameter 
        end
    end                                                                     % End loop over indirect dimension


   
    if secDim > 1                                                           % update each block in the jacobian according to the 2-D parametrization
        dataPoints = size(data(indMin:indMax,:),1);                         % Get number of datapoints
        regularizerPoints = dYdph0 - dataPoints;                            % Get number of points in regularizer
        [dYdph0]        = updateJacobianBlock(dYdph0,'ph0', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for ph0 parameter
        [dYdph1]        = updateJacobianBlock(dYdph1,'ph1', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for ph1 parameter
        [dYdgaussLB]    = updateJacobianBlock(dYdgaussLB,'gaussLB', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for gaussLB parameter
        [dYdlorentzLB]  = updateJacobianBlock(dYdlorentzLB,'lorentzLB', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for lorentzLB parameter
        [dYdfreqShift]  = updateJacobianBlock(dYdfreqShift,'freqShift', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for freqShift parameter
        [dYdmetAmpl]    = updateJacobianBlock(dYdmetAmpl,'metAmpl', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise); % Update jacobian block for metAmpl parameter
        [dYdbaseAmpl]   = updateJacobianBlock(dYdbaseAmpl,'baseAmpl', parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise);  % Update jacobian block for baseAmpl parameter       
    end

    if nBaselineComps ~= 0                                                  % Has baseline?
        jac = cat(2, dYdph0, dYdph1, dYdgaussLB, dYdlorentzLB, dYdfreqShift, dYdmetAmpl, dYdbaseAmpl); % Create final jacobian
    else
        jac = cat(2, dYdph0, dYdph1, dYdgaussLB, dYdlorentzLB, dYdfreqShift, dYdmetAmpl);   % Create final jacobian
    end

    if strcmp(SignalPart,'R') 
        jac         = real(jac);                                            % Take real part jacobian
    end
    if strcmp(SignalPart,'I') 
        jac         = imag(jac);                                            % Take imaginary part jacobian
    end
    if strcmp(SignalPart,'RI') 
        jac = cat(1, real(jac), imag(jac));                                 % Concatenate real and imaginary part
    end
    if strcmp(SignalPart,'A') 
        jac         = abs(jac);                                             % Take magnitude of jacobian
    end
    if strcmp(SignalPart,'C')                                               % just return the complex for CRLBs
       
    end
    jac = (-1) * jac;                                                       % Needed to match numerical jacobian

    if secDim == 1 && ~isempty(NormNoise)                                   % 1D jacobians have not been normalized yet
        dataPoints = size(data(indMin:indMax,:),1);                         % Get number of datapoints
        Sigma = (NormNoise * sqrt(dataPoints));                             % Calculate normalization factor for data
        if ~strcmp(SignalPart,'C')
            Sigma = Sigma^2;                                                % Square if not for CRLB
        end
        jac(1:dataPoints,:) = jac(1:dataPoints,:) / Sigma;                  % Normalize data part
        regularizerPoints = size(jac,1) - dataPoints;                       % Get number of points in regularizer
        Sigma = (NormNoise * sqrt(dataPoints + regularizerPoints));         % Calculate normalization factor for regularizer
        if ~strcmp(SignalPart,'C')
            Sigma = Sigma^2;                                                % Square if not for CRLB
        end
        jac(dataPoints+1:end,:) = jac(dataPoints+1:end,:) / Sigma;          % Normalize regularizer part
    end
end

function [Y, baseline, metabs,regu] = forwardModel(x, basisSet, baselineBasis, ppm, t, Reg, parametrizations)
% This function generates forward model and regularizer for the generalized physics model
%
%   USAGE:
%       [Y, baseline, metabs,regu] = forwardModel(x, basisSet, baselineBasis, ppm, t, Reg, parametrizations)
%
%   INPUTS:
%       x                = x vector with parameters to optimize
%       basisSet         = basis set struct
%       baselineBasis    = baseline basis set
%       ppm              = ppm axis
%       t                = time vector
%       Reg              = regularizer flag
%       parametrizations = parameter struct
%
%   OUTPUTS:
%       Y                = full prediction of forward model
%       baseline         = baseline prediction of forward model  
%       metabs           = basis function prediction of forward model
%       regu             = regularizer prediction of the forward model
%
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Calculate forward model and regularizer

    % Initialize output
    Y = [];                                                                 % Initialize model
    baseline = [];                                                          % Initialize baseline
    metabs = [];                                                            % Initialize basis function
    regu = [];                                                              % Initialize regularizer

    fidsBasis = basisSet.fids;                                              % Get basis function
    nBasisFcts = sum(basisSet.includeInFit(end,:));                         % Get number basis functions
    secDim = size(fidsBasis,3);                                             % Get number spectra along indirect dimension
    
    inputParams = x2pars(x, secDim, parametrizations);                      % Convert x vector to parameter struct

    % Loop over indirect dimension
    for sD = 1 : secDim % Loop over indirect dimension
        fids        = squeeze(fidsBasis(:,:,sD));                           % Get time domain data
        gaussLB = squeeze(inputParams.gaussLB(sD,:));                       % Get gaussLB parameter
        lorentzLB = squeeze(inputParams.lorentzLB(sD,:));                   % Get lorentzLB parameter
        freqShift = squeeze(inputParams.freqShift(sD,:));                   % Get freqShift parameter
        metAmpl = squeeze(inputParams.metAmpl(sD,:))';                      % Get metAmpl parameter
        baseAmpl = squeeze(inputParams.baseAmpl(sD,:))';                    % Get baseAmpl parameter
        ph0 = squeeze(inputParams.ph0(sD));                                 % Get ph0 parameter
        ph1 = squeeze(inputParams.ph1(sD));                                 % Get ph1 parameter
        
        timeDomainMultiplier = zeros(size(fids));                           % Setup time domain multiplier
        for ll = 1:nBasisFcts                                               % Loop over basis functions
            timeDomainMultiplier(:,ll) = exp(-(1i*freqShift(ll) + lorentzLB(ll) + gaussLB.^2.*t).*t)';  % Apply freqshifts, lorentzianLB, and gaussLB  
        end                                                                 % End loop over basis functions
        
        Fl = timeDomainMultiplier .* fids;                                  % Apply time domain multiplier to basis functions
        specs = fftshift(fft(Fl, [], 1),1);                                 % Convert to frequency domain             
        mets = specs * metAmpl;                                             % Multiply with basis function amplitude estimates                                                                                        
        bl = baselineBasis * baseAmpl;                                      % Multiply with baseline basis function amplitude estimates

        T_ph = exp(-1j .* (ph0 + ph1.*ppm)');                               % Create phase evolution

        if ~isempty(bl)                                                     % Has baseline?
            Y = cat(2,Y,T_ph .* (mets + bl));                               % Final model (phase evolution * estimates)
            baseline = cat(2,baseline,T_ph .* bl);                          % Baseline (phase evolution * estimates)
        else
            Y = cat(2,Y,T_ph .* mets);                                      % Final model
            baseline = cat(2,baseline,zeros(size(specs)));                  % zero baseline
        end
        metabs = cat(3,metabs,repmat(T_ph, [1 size(specs,2)]) .* specs .* repmat(metAmpl', [size(specs,1) 1]));  % Final basis functions (phase evolution * estimates)
    end                                                                     % End loop over indirect dimension

    if Reg                                                                  % Add parameter regularization
        for sD = 1 : secDim                                                 % Loop over indirect dimension
            gaussLB = squeeze(inputParams.gaussLB(sD,:));                   % Get gaussLB parameter
            lorentzLB = squeeze(inputParams.lorentzLB(sD,:));               % Get lorentzLB parameter
            freqShift = squeeze(inputParams.freqShift(sD,:));               % Get freqShift parameter
            metAmpl = squeeze(inputParams.metAmpl(sD,:));                   % Get metAmpl parameter
            baseAmpl = squeeze(inputParams.baseAmpl(sD,:));                 % Get baseAmpl parameter
            ph0 = squeeze(inputParams.ph0(sD));                             % Get ph0 parameter
            ph1 = squeeze(inputParams.ph1(sD));                             % Get ph1 parameter
            [ph0Reg]        = addParameterRegularization([],'ph0', parametrizations,ph0,0); % Calculate regularizer for ph0 parameter
            [ph1Reg]        = addParameterRegularization([],'ph1', parametrizations,ph1,0); % Calculate regularizer for ph1 parameter
            [gaussLBReg]    = addParameterRegularization([],'gaussLB', parametrizations,gaussLB,0); % Calculate regularizer for gaussLB parameter
            [lorentzLBReg]  = addParameterRegularization([],'lorentzLB', parametrizations,lorentzLB,0); % Calculate regularizer for lorentzLB parameter
            [freqShiftReg]  = addParameterRegularization([],'freqShift', parametrizations,freqShift,0); % Calculate regularizer for freqShift parameter
            [metAmplReg]    = addParameterRegularization([],'metAmpl', parametrizations,metAmpl,0); % Calculate regularizer for metAmpl parameter
            [baseAmplReg]   = addParameterRegularization([],'baseAmpl', parametrizations,baseAmpl',0); % Calculate regularizer for baseAmpl parameter
            regu = cat(2,regu,[ph0Reg , ph1Reg, gaussLBReg, lorentzLBReg, freqShiftReg, metAmplReg, baseAmplReg]); % Concatenate regularizer
        end
    end                                                                     % End loop over indirect dimension
end

function paramStruct = x2pars(x, secDim, parametrizations)      
% This function converts a 1-D x vector into a parameter struct
%
%   USAGE:
%       paramStruct = x2pars(x, secDim, parametrizations)   
%
%   INPUTS:
%       x                = x vector with parameters to optimize
%       secDim           = number of spectra in indirect dimension
%       parametrizations = parameterization options
%
%   OUTPUTS:
%       paramStruct      = parameter struct
%
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Converts a 1-D x vector
    pars = fields(parametrizations);                                        % Get parameter names       
    for ff = 1 : length(pars)                                               % Loop over parameters                
        if parametrizations.(pars{ff}).start                                % Start with a parameter struct of 1-D vectors according to the indices
            paramStruct.(pars{ff}) = x(parametrizations.(pars{ff}).start:parametrizations.(pars{ff}).end);  % Get parameters from x vector
        else
            paramStruct.(pars{ff})  =[];
        end

        % Reshape the 1-D vectors according to the number of basis functions,
        % number of subspectra, and parametrization
        switch pars{ff}                                                     % Switch for parameter names
            case {'ph0','ph1','gaussLB'}                                    % Parameters that appear once per spectrum
                if strcmp(parametrizations.(pars{ff}).type,'free')
                    paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),secDim,1);
                end
                if strcmp(parametrizations.(pars{ff}).type,'fixed')
                    paramStruct.(pars{ff}) = repmat(paramStruct.(pars{ff}),[secDim,1]);
                end
                if strcmp(parametrizations.(pars{ff}).type,'dynamic')
                    paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),size(parametrizations.(pars{ff}).lb));
                    for rp = 1 : length(parametrizations.(pars{ff}).parameterNames)
                        paramStruct.([pars{ff} 'Reparametrization']).(parametrizations.(pars{ff}).parameterNames{rp}) = paramStruct.(pars{ff})(rp,:);
                    end
                    paramStruct.(pars{ff}) = parametrizations.metAmpl.fun.fun(paramStruct.(pars{ff}),parametrizations.(pars{ff}).modulator);
                end
            case {'metAmpl', 'freqShift', 'lorentzLB','baseAmpl'}            % Parameters that appear once per basis function
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
                if strcmp(parametrizations.(pars{ff}).type,'none')
                    paramStruct.(pars{ff}) = zeros(secDim,1);
                end
        end
    end                                                                     % End loop over parameters          
end

function [x,indexStruct] = pars2x(paramStruct)
% This function converts a parameter struct into a 1-D x vector that can be
% passed on to solvers. It also defines start and end indices for
% easier identification
%
%   USAGE:
%       [x,indexStruct] = pars2x(paramStruct)
%
%   INPUTS:
%       paramStruct      = parameter struct
%
%   OUTPUTS:
%       x                = x vector with parameters to optimize
%       indexStruct      = struct with position indicies
%
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Converts parameter struct

    pars = fields(paramStruct);                                             % Get parameter names
    x = [];                                                                 % Initialize 1-D x vector
    for ff = 1 : length(pars)                                               % Loop over parameters
        if ismember(pars{ff},{'ph0','ph1','gaussLB','lorentzLB','freqShift','metAmpl','baseAmpl'}) % Skip new parameters from dynamic parameterization bc they should not turn up in the x vector
            if ~isempty(paramStruct.(pars{ff}))                             % baseAmpl is empty when no baseline included          
                if isempty(x)                                               % Setup index struct start value
                    indexStruct.(pars{ff}).start = 1;                       % Set index struct start value to 1
                else
                    indexStruct.(pars{ff}).start = length(x)+1;             % Set index struct start value to length + 1
                end
                x = cat(2,x,reshape(paramStruct.(pars{ff}),1,[]));          % Add new parameters to x vector
                indexStruct.(pars{ff}).end = length(x);                     % Set index struct end value to length
            end
        end
    end                                                                     % End loop over parameters   
end

function [f,g,h] = fminunc_wrapper(x,F,GJ,H)
% This function sets up the fminunc_wrapper for MATLAB
%%
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

%% Functions needed for 2D modeling and regularization 
function parameterMatrix = addParameterRegularization(parameterMatrix,parameterName, parametrizations,inputParams,jacobian)
% This function adds regularization terms to a parameter matrix
%
%   USAGE:
%       parameterMatrix = addParameterRegularization(parameterMatrix,parameterName, parametrizations,inputParams,jacobian)
%
%   INPUTS:
%       parameterMatrix  = input parameter matrix without regularizer term
%       parameterName    = name of parameter
%       parametrizations = paramterization options
%       inputParams      = parameter values
%       jacobian         = jacobian flag
%
%   OUTPUTS:
%       parameterMatrix  = updated parameter matrix (jacobian or penalty term)
%
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Add regularizer to parameter matrix
    if ~strcmp(parametrizations.(parameterName).RegFun,'')                  % Don't apply regularization to the parameter       
        if jacobian
           RegMatrix = parametrizations.(parameterName).RegFun.fun(parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1);
           RegMatrix = sqrt(parametrizations.(parameterName).RegPar)*RegMatrix;
           parameterMatrix = cat(1,parameterMatrix,RegMatrix);
        else
           RegMatrix = parametrizations.(parameterName).RegFun.fun(parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1);
           RegMatrix = sqrt(parametrizations.(parameterName).RegPar)*RegMatrix*inputParams;
           parameterMatrix =  RegMatrix;
        end
       
    end
    if jacobian                                                             % Parameter matrix is jacobian
        pars = fields(parametrizations);                                    % Get parameter names
        numberOfParameters = [];                                            % Initialize number of parameters
        for ff = 1 : length(pars)                                           % Loop over parameters       
            numberOfParameters(end+1) = parametrizations.(pars{ff}).end-parametrizations.(pars{ff}).start + 1; % Calculate number of parameters       
        end                                                                 % End loop over parameters       
        if (parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1) <= max(numberOfParameters) && strcmp(parametrizations.(parameterName).RegFun,'')  % Add zeros if jacobian is too short 
            parameterMatrix = cat(1,parameterMatrix,zeros(max(numberOfParameters),parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1));    %Add correct number of zeros to the end
        end
    end
end

function dYdX = updateJacobianBlock(dYdX,parameterName, parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise)
% This function updates jacobians according to the parametrization
%
%   USAGE:
%       dYdX = updateJacobianBlock(dYdX,parameterName, parametrizations,inputParams,SignalPart,regularizerPoints,NormNoise)
%
%   INPUTS:
%       dYdX             = jacobian 
%       parameterName    = name of parameter
%       parametrizations = paramterization options
%       inputParams      = parameter values
%       SignalPart       = optimization signal part
%       regularizerPoints= number of points in the regularizer
%       NormNoise        = 1-point noise estimate
%
%   OUTPUTS:
%       dYdX             = re-organized jacobian
%
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%% Update jacobian blocks
    % Get dimensions
    nPoints = size(dYdX,1);
    nLines = size(dYdX,2);
    secDim = size(dYdX,3);

    if ~isempty(NormNoise)
        Sigma = (NormNoise * sqrt(nPoints));                                % Calculate normalization factor for data
        if ~strcmp(SignalPart,'C')
            Sigma = Sigma.^2;                                               % Square if not for CRLBs
        end
        Sigma = repmat(Sigma, [nPoints nLines 1]);                          % Repeat according to dimensions
        dYdX = squeeze(dYdX);                                               % Remove zero dimensions 
        dYdX(1:nPoints,:) = dYdX(1:nPoints,:,:) ./ Sigma;                   % Normalize data part
        if regularizerPoints > 0                                            % Has regularizer
            Sigma = (NormNoise * sqrt(nPoints + regularizerPoints));        % Calculate normalization factor for regularizer
            if ~strcmp(SignalPart,'C')
                Sigma = Sigma.^2;                                           % Square if not for CRLBs
            end
            Sigma = repmat(Sigma, [regularizerPoints nLines 1]);            % Repeat according to dimensions
            dYdX(nPoints+1:end,:) = nPoints(nPoints+1:end,:,:) ./ Sigma;    % Normalize regularizer part
        end
    end

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
            secDim = size(dYdX,4);                                          
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
                nLines = secDim;
            end
            dYdX = reshape(dYdX,[],nLines);
        case  4 % E.g. free metAmpls 
            dYdX = permute(dYdX,[1 3 4 2]);
            dYdX = reshape(dYdX,[],secDim*nLines);
    end
end