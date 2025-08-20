%%  GeneralizedBasicPhysicsModel.m
%   This function contains the generalized physics model allowing for 2D
%   modeling of arbitray MRS data in different domains. It is desinged with
%   the highest possible felxibility.
%   It integreates into the OspreyFitObj environment and allows for easy
%   changes of solvers. You can also use this function as a template for
%   your own Physics model. For full functionalty you need to define the
%   same functions and handels described in the first section. This
%   includes the following 7 functions: 
%
%   lossFunction    - calculates loss function for optimizer
%   forwardGradient - calculates the forward gradient for optimizer
%   forwardJacobian - calculates the forward jacobian for optimizer
%   forwardModel    - converts parameter struct into model prediction
%   x2pars          - converts x vector to parameter struct
%   pars2x          - converts parameter struct to x vector
%   fminunc_wrapper - wrapper for MATLAB's fminunc
%
%   There are additional functions that are not included as handles. They
%   are requried to construct the jacobian for 2D modeling and perform
%   parameter regularization:
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
    fh.lossFunction     = @lossFunction;
    fh.forwardGradient  = @forwardGradient;
    fh.forwardJacobian  = @forwardJacobian;
    fh.forwardModel     = @forwardModel;
    fh.x2pars           = @x2pars;
    fh.pars2x           = @pars2x;
    fh.fminunc_wrapper  = @fminunc_wrapper;
end

%% Define lossfunction, forward gradient, forward jacobian, forward model, x2pars, pars2x, fminunc_wrapper

function sse = lossFunction(x, data, NoiseSD, basisSet, baselineBasis, ppm, t, fitRange, fitGap, SignalPart, Domain, SSE, Reg, parametrizations)
% This function generates the output for the solver according to the
% loss function. This includes all settings described in the model procedure
% json.
%
%   USAGE:
%       sse = lossFunction(x, data, NoiseSD, basisSet, baselineBasis, ppm, t, fitRange, SignalPart, Domain, SSE, Reg, parametrizations)
%
%   INPUTS:
%       x                = x vector with parameters to optimize
%       data             = ppm axis
%       NoiseSD          = standard deviation of the noise
%       basisSet         = basis set struct
%       baselineBasis    = baseline basis set
%       ppm              = ppm axis
%       t                = time vector
%       fitRange         = model range
%       fitGap           = GAP parameter
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
    [prediction, ~, ~, regu, penaltyTerm,penaltyTermSoftConstraint]  = forwardModel(x, ...                         % x vector with parameters to optimize
                                                                         basisSet, ...                  % basis set struct
                                                                         baselineBasis, ...             % baseline basis set
                                                                         ppm, ...                       % ppm axis
                                                                         t, ...                         % time vector
                                                                         Reg, ...                       % regularizer flag
                                                                         parametrizations);             % parameter struct


    if strcmp(Domain,'FD')                                                  % frequency domain optimization
        [indMin, indMax] = ppmToIndex(ppm, fitRange.FD);                       % Get fit range ppm indices
        data        = fftshift(fft(data, [], 1),1);                         % Convert data frequency domain data
        data        = data(indMin:indMax,:);                                % Cut out fit range from data
        prediction  = prediction(indMin:indMax,:);                          % Cut out fit range from prediction
    end
    if strcmp(Domain,'TD')                                                  % time domain optimization
        indMin = fitRange.TD(1);                                            % Get time domain indices
        indMax = fitRange.TD(2);                                            % Get time domain indices
        prediction  = ifft(ifftshift(prediction,1), [], 1);                 % Convert prediction to time domain
        regu        = ifft(ifftshift(regu,1), [], 1);                       % Convert regularizer to time domain
        data        = data(indMin:indMax,:);                                % Cut out fit range from data
        prediction  = prediction(indMin:indMax,:);                          % Cut out fit range from prediction
    end
    if strcmp(Domain,'FDTD')                                                % Simualtneous frequency and time domain optimization
        [indMinFD, indMaxFD] = ppmToIndex(ppm, fitRange.FD);                % Get fit range ppm indices
        dataFD        = fftshift(fft(data, [], 1),1);                       % Convert data frequency domain data
        dataFD        = dataFD(indMinFD:indMaxFD,:);                        % Cut out fit range from data
        predictionFD  = prediction(indMinFD:indMaxFD,:);                    % Cut out fit range from prediction

        indMinTD      = fitRange.TD(1);                                     % Get time domain indices
        indMaxTD      = fitRange.TD(2);                                     % Get time domain indices
        predictionTD  = ifft(ifftshift(prediction,1), [], 1);               % Convert prediction to time domain
        dataTD        = data(indMinTD:indMaxTD,:);                          % Cut out fit range from data
        predictionTD  = predictionTD(indMinTD:indMaxTD,:);                    % Cut out fit range from prediction
        data          = cat(1,dataFD,dataTD);                               % Combine FD and TD data
        prediction    = cat(1,predictionFD,predictionTD);                   % Combine FD and TD data
    end
    residual     = data - prediction;                                       % Calculate residual

    if ~isempty(fitGap) && strcmp(Domain,'FD')
        [GapindMin, GapindMax] = ppmToIndex(ppm, fitGap);                   % Get fit Gap range ppm indices
        residual((GapindMin-indMin):(GapindMax-indMin),:)=[];
    end

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

    if ~isempty(NoiseSD)                                                    % Don't normalize residual if GradientCheck is performed
        if strcmp(Domain,'FD')
            Sigma = NoiseSD.FD;                                             % Get sigma squared
            Sigma = repmat(Sigma, [size(residual,1) 1]);                    % Repeat according to dimensions
        end
        if strcmp(Domain,'TD')
            Sigma = NoiseSD.TD;                                             % Get sigma squared
            Sigma = repmat(Sigma, [size(residual,1) 1]);                    % Repeat according to dimensions
        end
        if strcmp(Domain,'FDTD')
            SigmaFD = NoiseSD.FD;                                           % Get sigma squared
            SigmaFD = repmat(SigmaFD, [size(predictionFD,1) 1]);            % Repeat according to dimensions
            SigmaTD = NoiseSD.TD;                                           % Get sigma squared
            SigmaTD = repmat(SigmaTD, [size(predictionTD,1) 1]);            % Repeat according to dimensions
            Sigma  = cat(1,SigmaFD,SigmaTD);                                % Combine frequency and time domain     
        end
        residual = residual ./ Sigma;                                       % Normalize residual
        if size(regu,1) > 0                                                 % Has regularizer
            if strcmp(Domain,'FD')
                Sigma = NoiseSD.FD;                                         % Get sigma squared
            end
            if strcmp(Domain,'TD')
                Sigma = NoiseSD.TD;                                         % Get sigma squared
            end
            if strcmp(Domain,'FDTD')
                Sigma = (NoiseSD.FD + NoiseSD.TD)/2;                        % Get sigma squared
                Sigma = NoiseSD.FD;                                         % Get sigma squared
            end
            Sigma = repmat(Sigma, [size(regu,1) 1]);                        % Repeat according to dimensions
            regu = regu ./ Sigma;                                           % Normalize regularizer
        end
    end

    penaltyTermSoftConstraintToAdd = zeros(size(penaltyTermSoftConstraint,1),size(residual,2)); % For 2D data we need to create a matrix
    penaltyTermSoftConstraintToAdd(:,1) = penaltyTermSoftConstraint;        % Update according to calculations 

    residual  = cat(1,residual,-regu);                                      % Add regularizer term
    residual  = cat(1,residual,penaltyTerm);                                % Add penalty term
    residual  = cat(1,residual,penaltyTermSoftConstraintToAdd);                  % Add soft constraint penalty term
    residual = reshape(residual,[],1);                                      % Reshape the multidimenisonal case

    switch SSE                                                              % Switch for return to solver
        case 'res'
            sse = residual;                                                 % Return residual
        case 'sos'
            sse = sum(residual.^2);                                         % Return sum of squares
    end
end

function grad = forwardGradient(x, data, NoiseSD, basisSet, baselineBasis, ppm, t, fitRange, SignalPart,Reg, parametrizations)
% This function generates the gradient for the generalized physics model
%
%   USAGE:
%       grad = forwardGradient(x, data, NoiseSD, basisSet, baselineBasis, ppm, t, fitRange,SignalPart,Reg, parametrizations)
%
%   INPUTS:
%       x                = x vector with parameters to optimize
%       data             = ppm axis
%       NoiseSD          = standard deviation of the noise
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

    [indMin, indMax] = ppmToIndex(ppm, fitRange.FD);                           % Get fit range ppm indices

    % Run forward model to get prediction and regularizer output
    [prediction,~,~,regu]  = forwardModel(x, ...                            % x vector with parameters to optimize
                                          basisSet, ...                     % basis set struct
                                          baselineBasis, ...                % baseline basis set
                                          ppm, ...                          % ppm axis
                                          t, ...                            % time vector
                                          Reg, ...                          % regularizer flag
                                          parametrizations);                % parameter struct

    residual     = data - prediction;                                       % Calculate residual

    if ~isempty(NoiseSD)                                                    % Don't normalize residual if GradientCheck is performed
        Sigma = NoiseSD;                                                    % Get sigma squared
        Sigma = repmat(Sigma, [size(data,1) 1]);              % Repeat according to dimensions
        residual = residual ./ Sigma;                                       % Normalize residual
        if size(regu,1) > 0                                                 % Has regularizer
            Sigma = NoiseSD;                                                % Get sigma squared
            Sigma = repmat(Sigma, [size(regu,1) 1]);                        % Repeat according to dimensions
            regu = regu ./ Sigma;                                           % Normalize regularizer
        end
    end

    residual      = residual(indMin:indMax,:);                              % Cut frequency range from residual
    residual      = cat(1,residual,-regu);                                  % Concatenate regularizer

    % Generate jacobian
    jac = forwardJacobian(x, NoiseSD, basisSet, baselineBasis, ppm, t, fitRange.FD,SignalPart,Reg, parametrizations);
    if strcmp(SignalPart,'R')
        jac         = real(jac);                                  % Take real part jacobian
        residual    = real(residual);
    end
    if strcmp(SignalPart,'I')
        jac         = imag(jac);                                  % Take imaginary part jacobian
        residual         = imag(residual);
    end
    if strcmp(SignalPart,'RI')
        jac = cat(1, real(jac), imag(jac));                  % Concatenate real and imaginary part
        residual = cat(1, real(residual), imag(residual));
    end
    if strcmp(SignalPart,'A')
        jac         = abs(jac);                                   % Take magnitude of residual
        residual         = abs(residual);
    end
    if strcmp(SignalPart,'C')                                               % just return the complex for CRLBs

    end

    grad = sum((residual).*(-conj(jac)) + (-jac .* conj(residual)));    % Calculate gradient

    grad = grad';
end

function jac = forwardJacobian(x, data, NoiseSD, basisSet, baselineBasis, ppm, t, fitRange, fitGap, SignalPart, Domain, Reg, parametrizations)
% This function generates the jacobian for the generalized physics model
%
%   USAGE:
%       jac = forwardJacobian(x, data, NoiseSD, basisSet, baselineBasis, ppm, t, fitRange,SignalPart,Reg, parametrizations)
%
%   INPUTS:
%       x                = x vector with parameters to optimize
%       NoiseSD          = standard deviation of the noise
%       basisSet         = basis set struct
%       baselineBasis    = baseline basis set
%       ppm              = ppm axis
%       t                = time vector
%       fitRange         = model range
%       SignalPart       = optimization signal part
%       Domain           = optimization domain
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
    secDim = size(fidsBasis,3);                                             % Get number spectra along indirect dimension

    inputParams = x2pars(x, secDim, parametrizations);                      % Convert x vector to parameter struct

    % Run forward model to get prediction and regularizer output
    prediction  = forwardModel(x, ...                                       % x vector with parameters to optimize
                             basisSet, ...                                  % basis set struct
                             baselineBasis, ...                             % baseline basis set
                             ppm, ...                                       % ppm axis
                             t, ...                                         % time vector
                             Reg, ...                                       % regularizer flag
                             parametrizations);                             % parameter struct


    % Initialize gaussian lw factor
    GaussianLWfactor = (1/2*(pi*sqrt(2*log(2)))^2);                         % Definition allows direct conversion to T2


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
            timeDomainMultiplier(:,ll) = exp(-(1i*2*pi*freqShift(ll) + ...       % Apply freqshifts
                                                (1/pi)*lorentzLB(ll) + ...         % Apply lorentzianLB
                                                GaussianLWfactor * gaussLB.^2.*t).*t)';        % Apply gaussLB
        end                                                                % End loop over basis functions

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
        dYdfreqShift    = cat(3,dYdfreqShift,T_ph_basis .* (-1j) .* (2*pi) .* Fmett .* metAmpl'); % Partial derivative wrt freqShift
        dYdlorentzLB    = cat(3,dYdlorentzLB,T_ph_basis .* (1/pi) .* Fmett .* metAmpl');         % Partial derivative wrt lorentzLB
        dYdgaussLB      = cat(3,dYdgaussLB,T_ph .* GaussianLWfactor .* Fmett2gauss * metAmpl);              % Partial derivative wrt gaussLB
        dYdph0          = cat(3,dYdph0,(-1j) .* prediction(:,sD));                      % Partial derivative wrt ph0
        dYdph1          = cat(3,dYdph1,(-1j) .* ppm' .* prediction(:,sD));              % Partial derivative wrt ph1
        if nBaselineComps ~= 0                                                          % Has baseline?
            dYdbaseAmpl           = cat(3,dYdbaseAmpl,T_ph_baseline .* baselineBasis);  % Partial derivative wrt baseAmpl
        else
            dYdbaseAmpl = cat(3,dYdbaseAmpl,[]);                                        % Empty partial derivative when no baseline is defined
        end

    end                                                                     % End loop over indirect dimension


    %Pick optimization domain
    if strcmp(Domain,'FD') || strcmp(Domain,'FDTD')                           % frequency domain optimization
        [indMin, indMax] = ppmToIndex(ppm, fitRange.FD);                    % Get fit range ppm indices
    end
    if strcmp(Domain,'TD')                                                  % time domain optimization
        indMin = fitRange.TD(1);                                            % Get time domain indices
        indMax = fitRange.TD(2);                                            % Get time domain indices
        dYdph0  = ifft(ifftshift(dYdph0,1), [], 1);                         % Convert partial derivative wrt ph0 to time domain
        dYdph1  = ifft(ifftshift(dYdph1,1), [], 1);                         % Convert partial derivative wrt ph1 to time domain
        dYdgaussLB  = ifft(ifftshift(dYdgaussLB,1), [], 1);                 % Convert partial derivative wrt gaussLB to time domain
        dYdlorentzLB  = ifft(ifftshift(dYdlorentzLB,1), [], 1);             % Convert partial derivative wrt lorentzLB to time domain
        dYdfreqShift  = ifft(ifftshift(dYdfreqShift,1), [], 1);             % Convert partial derivative wrt freqShift to time domain
        dYdmetAmpl  = ifft(ifftshift(dYdmetAmpl,1), [], 1);                 % Convert partial derivative wrt metAmpl to time domain
        if nBaselineComps ~= 0                                              % Has baseline?
            dYdbaseAmpl  = ifft(ifftshift(dYdbaseAmpl,1), [], 1);           % Convert partial derivative wrt baseAmpl to time domain
        end
    end

    if strcmp(Domain,'FD') || strcmp(Domain,'TD')                           % frequency domain optimization
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

        % Check if there is a gap in the fit range, if yes introduce it
        if ~isempty(fitGap) && strcmp(Domain,'FD')
            [GapindMin, GapindMax] = ppmToIndex(ppm, fitGap);                   % Get fit Gap range ppm indices
                
            dYdph0((GapindMin-indMin):(GapindMax-indMin),:,:)             = []; % Cut out the gap
            dYdph1((GapindMin-indMin):(GapindMax-indMin),:,:)             = []; % Cut out fit gap
            dYdgaussLB((GapindMin-indMin):(GapindMax-indMin),:,:)         = []; % Cut out fit gap
            dYdlorentzLB((GapindMin-indMin):(GapindMax-indMin),:,:,:)     = []; % Cut out fit gap
            dYdfreqShift((GapindMin-indMin):(GapindMax-indMin),:,:,:)     = []; % Cut out fit gap
            dYdmetAmpl((GapindMin-indMin):(GapindMax-indMin),:,:,:)       = []; % Cut out fit gap
            if nBaselineComps ~= 0                                              % Has baseline?
                dYdbaseAmpl((GapindMin-indMin):(GapindMax-indMin),:,:,:)  = []; % Cut out fit gap
            end
        end

    else                                                                    % Simulatneous frequency and time domain optimization             
        indMinTD = fitRange.TD(1);                                            % Get time domain indices
        indMaxTD = fitRange.TD(2);                                            % Get time domain indices
        dYdph0TD  = ifft(ifftshift(dYdph0,1), [], 1);                         % Convert partial derivative wrt ph0 to time domain
        dYdph1TD  = ifft(ifftshift(dYdph1,1), [], 1);                         % Convert partial derivative wrt ph1 to time domain
        dYdgaussLBTD  = ifft(ifftshift(dYdgaussLB,1), [], 1);                 % Convert partial derivative wrt gaussLB to time domain
        dYdlorentzLBTD  = ifft(ifftshift(dYdlorentzLB,1), [], 1);             % Convert partial derivative wrt lorentzLB to time domain
        dYdfreqShiftTD  = ifft(ifftshift(dYdfreqShift,1), [], 1);             % Convert partial derivative wrt freqShift to time domain
        dYdmetAmplTD  = ifft(ifftshift(dYdmetAmpl,1), [], 1);                 % Convert partial derivative wrt metAmpl to time domain
        if nBaselineComps ~= 0                                              % Has baseline?
            dYdbaseAmplTD  = ifft(ifftshift(dYdbaseAmpl,1), [], 1);           % Convert partial derivative wrt baseAmpl to time domain
        end
        % Combine frequency and time domain jacobians
        dYdph0          = cat(1,dYdph0(indMin:indMax,:,:),dYdph0TD(indMinTD:indMaxTD,:,:)); % Cut out fit ranges and combine domains
        dYdph1          = cat(1,dYdph1(indMin:indMax,:,:),dYdph1TD(indMinTD:indMaxTD,:,:));                            % Cut out fit range
        dYdgaussLB      = cat(1,dYdgaussLB(indMin:indMax,:,:),dYdgaussLBTD(indMinTD:indMaxTD,:,:));                        % Cut out fit range
        dYdlorentzLB    = cat(1,dYdlorentzLB(indMin:indMax,:,:,:),dYdlorentzLBTD(indMinTD:indMaxTD,:,:,:));                    % Cut out fit range
        dYdfreqShift    = cat(1,dYdfreqShift(indMin:indMax,:,:,:),dYdfreqShiftTD(indMinTD:indMaxTD,:,:,:));                    % Cut out fit range
        dYdmetAmpl      = cat(1,dYdmetAmpl(indMin:indMax,:,:,:),dYdmetAmplTD(indMinTD:indMaxTD,:,:,:));                      % Cut out fit range
        if nBaselineComps ~= 0                                                  % Has baseline?
            dYdbaseAmpl           = cat(1,dYdbaseAmpl(indMin:indMax,:,:,:),dYdbaseAmplTD(indMinTD:indMaxTD,:,:,:));           % Cut out fit range
        end
    end

    [dYdpen] = calcPenalty(1,inputParams, parametrizations, secDim);       % Calculate expectation value penalty jacobian  
    nParams = size(dYdpen,1);                                              % Store number of parameters before grouping for soft constraint jacobian
    dYdpen  = updateAccordingToGrouping('penaltyJacobian',dYdpen,[], parametrizations,x);   % Update expectation value penalty according to grouping

    
    tempTerm = [];
    [~,ph0PenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'ph0', 1, nParams);             % Calculate jacobian for soft constraint penalty for ph0           
    [~,ph1PenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'ph1', 1, nParams);             % Calculate jacobian for soft constraint penalty for ph1
    [~,gaussLBPenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'gaussLB', 1, nParams);     % Calculate jacobian for soft constraint penalty for gaussLB
    [~,lorentzLBPenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'lorentzLB', 1, nParams); % Calculate jacobian for soft constraint penalty for lorentzLB
    [~,freqShiftPenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'freqShift', 1, nParams); % Calculate jacobian for soft constraint penalty for freqShift
    [~,metAmplPenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'metAmpl', 1, nParams);     % Calculate jacobian for soft constraint penalty for metAmpl
    if nBaselineComps ~= 0 % If no baseline is included
        [~,baseAmplPenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'baseAmpl', 1, nParams); % Calculate jacobian for soft constraint penalty for baseAmpl
    else
        baseAmplPenaltyJac = zeros(size(ph0PenaltyJac));        % No baseline components means zeros only
    end

    % update each block in the jacobian according to the 2-D parametrization
    [dYdph0]        = updateJacobianBlock(dYdph0,'ph0', parametrizations,inputParams,Domain,NoiseSD,dYdpen,fitRange,ph0PenaltyJac,Reg); % Update jacobian block for ph0 parameter
    [dYdph1]        = updateJacobianBlock(dYdph1,'ph1', parametrizations,inputParams,Domain,NoiseSD,dYdpen,fitRange,ph1PenaltyJac,Reg); % Update jacobian block for ph1 parameter
    [dYdgaussLB]    = updateJacobianBlock(dYdgaussLB,'gaussLB', parametrizations,inputParams,Domain,NoiseSD,dYdpen,fitRange,gaussLBPenaltyJac,Reg); % Update jacobian block for gaussLB parameter
    [dYdlorentzLB]  = updateJacobianBlock(dYdlorentzLB,'lorentzLB', parametrizations,inputParams,Domain,NoiseSD,dYdpen,fitRange,lorentzLBPenaltyJac,Reg); % Update jacobian block for lorentzLB parameter
    [dYdfreqShift]  = updateJacobianBlock(dYdfreqShift,'freqShift', parametrizations,inputParams,Domain,NoiseSD,dYdpen,fitRange,freqShiftPenaltyJac,Reg); % Update jacobian block for freqShift parameter
    [dYdmetAmpl]    = updateJacobianBlock(dYdmetAmpl,'metAmpl', parametrizations,inputParams,Domain,NoiseSD,dYdpen,fitRange,metAmplPenaltyJac,Reg); % Update jacobian block for metAmpl parameter
    if nBaselineComps ~= 0
        [dYdbaseAmpl]   = updateJacobianBlock(dYdbaseAmpl,'baseAmpl', parametrizations,inputParams,Domain,NoiseSD,dYdpen,fitRange,baseAmplPenaltyJac,Reg);  % Update jacobian block for baseAmpl parameter
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

end


function [Y, baseline, metabs, regu, penaltyTerm,penaltyTermSoftConstraint] = forwardModel(x, basisSet, baselineBasis, ppm, t, Reg, parametrizations)
% This function generates forward model and regularizer for the generalized physics model
%
%   USAGE:
%       [Y, baseline, metabs, regu, penaltyTerm] = forwardModel(x, basisSet, baselineBasis, ppm, t, Reg, parametrizations)
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
%       penaltyTerm      = penalty term prediction of the forward model
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
    Y               = [];                                                   % Initialize model
    baseline        = [];                                                   % Initialize baseline
    metabs          = [];                                                   % Initialize basis function
    regu            = [];                                                   % Initialize regularizer
    penaltyTerm     = [];                                                   % Initialize penalty term
    penaltyTermSoftConstraint     = [];                                     % Initialize soft constraint penalty term

    fidsBasis = basisSet.fids;                                              % Get basis function
    nBasisFcts = sum(basisSet.includeInFit(end,:));                         % Get number basis functions
    secDim = size(fidsBasis,3);                                             % Get number spectra along indirect dimension

    inputParams = x2pars(x, secDim, parametrizations);                      % Convert x vector to parameter struct

    % Initialize gaussian lw factor
    GaussianLWfactor = (1/2*(pi*sqrt(2*log(2)))^2);                         % Definition allows direct conversion to T2

    % Loop over indirect dimension
    for sD = 1 : secDim % Loop over indirect dimension
        fids        = squeeze(fidsBasis(:,:,sD));                           % Get time domain data
        gaussLB     = squeeze(inputParams.gaussLB(sD,:));                   % Get gaussLB parameter
        lorentzLB   = squeeze(inputParams.lorentzLB(sD,:));                 % Get lorentzLB parameter
        freqShift   = squeeze(inputParams.freqShift(sD,:));                 % Get freqShift parameter
        metAmpl     = squeeze(inputParams.metAmpl(sD,:))';                  % Get metAmpl parameter
        baseAmpl    = squeeze(inputParams.baseAmpl(sD,:))';                 % Get baseAmpl parameter
        ph0         = squeeze(inputParams.ph0(sD));                         % Get ph0 parameter
        ph1         = squeeze(inputParams.ph1(sD));                         % Get ph1 parameter

        timeDomainMultiplier = zeros(size(fids));                           % Setup time domain multiplier
        for ll = 1:nBasisFcts                                               % Loop over basis functions
            timeDomainMultiplier(:,ll) = exp(-(1i*2*pi*freqShift(ll) + ...       % Apply freqshifts
                                                (1/pi)*lorentzLB(ll) + ...         % Apply lorentzianLB
                                                GaussianLWfactor * gaussLB.^2.*t).*t)';        % Apply gaussLB
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
            baseline = cat(2,baseline,zeros(size(specs,1),1));                  % zero baseline
        end
        metabs = cat(3,metabs,repmat(T_ph, [1 size(specs,2)]) .* specs .* repmat(metAmpl', [size(specs,1) 1]));  % Final basis functions (phase evolution * estimates)
    end                                                                     % End loop over indirect dimension

    if Reg                                                                  % Add parameter regularization
        for sD = 1 : secDim                                                 % Loop over indirect dimension
            gaussLB     = squeeze(inputParams.gaussLB(sD,:));               % Get gaussLB parameter
            lorentzLB   = squeeze(inputParams.lorentzLB(sD,:));             % Get lorentzLB parameter
            freqShift   = squeeze(inputParams.freqShift(sD,:));             % Get freqShift parameter
            metAmpl     = squeeze(inputParams.metAmpl(sD,:));               % Get metAmpl parameter
            baseAmpl    = squeeze(inputParams.baseAmpl(sD,:));              % Get baseAmpl parameter
            ph0         = squeeze(inputParams.ph0(sD));                     % Get ph0 parameter
            ph1         = squeeze(inputParams.ph1(sD));                     % Get ph1 parameter
            [ph0Reg]        = addParameterRegularization([],'ph0', parametrizations,ph0,0,secDim,sD); % Calculate regularizer for ph0 parameter
            [ph1Reg]        = addParameterRegularization([],'ph1', parametrizations,ph1,0,secDim,sD); % Calculate regularizer for ph1 parameter
            [gaussLBReg]    = addParameterRegularization([],'gaussLB', parametrizations,gaussLB,0,secDim,sD); % Calculate regularizer for gaussLB parameter
            [lorentzLBReg]  = addParameterRegularization([],'lorentzLB', parametrizations,lorentzLB,0,secDim,sD); % Calculate regularizer for lorentzLB parameter
            [freqShiftReg]  = addParameterRegularization([],'freqShift', parametrizations,freqShift,0,secDim,sD); % Calculate regularizer for freqShift parameter
            [metAmplReg]    = addParameterRegularization([],'metAmpl', parametrizations,metAmpl,0,secDim,sD); % Calculate regularizer for metAmpl parameter
            [baseAmplReg]   = addParameterRegularization([],'baseAmpl', parametrizations,baseAmpl',0,secDim,sD); % Calculate regularizer for baseAmpl parameter
            regu = cat(2,regu,[ph0Reg , ph1Reg, gaussLBReg, lorentzLBReg, freqShiftReg, metAmplReg, baseAmplReg]); % Concatenate regularizer
        end
    end                                                                     % End loop over indirect dimension

    penaltyTerm = calcPenalty(0,inputParams, parametrizations, secDim);    % Add penalty terms from expectation values
    penaltyTerm  = updateAccordingToGrouping('penalty',penaltyTerm,[], parametrizations,x); % Update expectation value penalty terms according to grouping

    ph0SoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'ph0', 1, size(inputParams.ph0,2));                   % Calculate soft constraint penalty term for ph0 parameter
    ph1SoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'ph1', 1, size(inputParams.ph1,2));                   % Calculate soft constraint penalty term for ph1 parameter
    gaussLBSoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'gaussLB', 1, size(inputParams.gaussLB,2));       % Calculate soft constraint penalty term for gaussLB parameter
    lorentzLBSoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'lorentzLB', 1, size(inputParams.lorentzLB,2)); % Calculate soft constraint penalty term for lorentzLB parameter
    freqShiftSoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'freqShift', 1, size(inputParams.freqShift,2)); % Calculate soft constraint penalty term for freqShift parameter
    metAmplSoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'metAmpl', 1, size(inputParams.metAmpl,2));       % Calculate soft constraint penalty term for metAmpl parameter
    if ~isempty(bl) 
        baseAmplSoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'baseAmpl', 1, size(inputParams.baseAmpl,2));% Calculate soft constraint penalty term for baseAmpl parameter
    else
        baseAmplSoftContstraint = [];   % No baseline parameter means no soft constraint penalty
    end
    penaltyTermSoftConstraint = cat(2,penaltyTermSoftConstraint,[ph0SoftContstraint , ph1SoftContstraint, gaussLBSoftContstraint, lorentzLBSoftContstraint, freqShiftSoftContstraint, metAmplSoftContstraint, baseAmplSoftContstraint]); % Concatenate penalty terms
    penaltyTermSoftConstraint  = updateAccordingToGrouping('penalty',penaltyTermSoftConstraint',[], parametrizations,x);    % Update soft constraint penalty terms according to grouping

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
                    if ~isempty(parametrizations.(pars{ff}).gr)  
                        idx=parametrizations.(pars{ff}).gr.idx;
                        paramStruct.(pars{ff}) = paramStruct.(pars{ff})(:,parametrizations.(pars{ff}).gr.idx); 
                    end
                end
                if strcmp(parametrizations.(pars{ff}).type,'fixed')
                    if ~isempty(parametrizations.(pars{ff}).gr)  
                        idx=parametrizations.(pars{ff}).gr.idx;
                        paramStruct.(pars{ff}) = paramStruct.(pars{ff})(idx); 
                    end
                    paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),1,[]);
                    paramStruct.(pars{ff}) = repmat(paramStruct.(pars{ff}),[secDim,1]);
                end
                if strcmp(parametrizations.(pars{ff}).type,'dynamic')                                        
                    if ~isempty(parametrizations.(pars{ff}).gr)  
                        idx=parametrizations.(pars{ff}).gr.idx;
                        nan_positions = isnan(reshape(parametrizations.(pars{ff}).gr.nan_marker,1,[]));
                        if sum(nan_positions) == 0
                            paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),size(parametrizations.(pars{ff}).lb));
                            paramStruct.(pars{ff})(:,:) = paramStruct.(pars{ff})(:,parametrizations.(pars{ff}).gr.idx);
                        else
                            temp_pars = nan(size(nan_positions));
                            temp_pars(~nan_positions) = paramStruct.(pars{ff});
                            paramStruct.(pars{ff}) = temp_pars;
                            paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),size(parametrizations.(pars{ff}).lb));
                             for gg_dyn = 1 : size(idx,1)
                                paramStruct.(pars{ff})(gg_dyn,:) = paramStruct.(pars{ff})(gg_dyn,idx(gg_dyn,:)); 
                             end
                        end 
                    else
                        paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),size(parametrizations.(pars{ff}).lb));
                    end
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
                temp_x =   reshape(paramStruct.(pars{ff}),1,[]);             % Remove nan values for unbalanced regrouping in dynamic models
                temp_x = temp_x(~isnan(temp_x));
                x = cat(2,x,temp_x);          % Add new parameters to x vector
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
function parameterMatrix = addParameterRegularization(parameterMatrix,parameterName, parametrizations,inputParams,jacobian,secDim,sD)
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
        if length(parametrizations.(parameterName).RegPar) ==1                  % Allows for regularization parameters for each sec dim entry
            RegPar = parametrizations.(parameterName).RegPar;
        else
            if ~isempty(sD)
                RegPar = parametrizations.(parameterName).RegPar(sD);
            else
                RegPar = parametrizations.(parameterName).RegPar;
            end
        end
       if strcmp(parametrizations.(parameterName).type,'free')
            RegMatrix = parametrizations.(parameterName).RegFun.fun((parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1)/secDim);
       else
            RegMatrix = parametrizations.(parameterName).RegFun.fun((parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1));
       end
        if jacobian
           if length(RegPar) == 1
                RegMatrix = repmat(sqrt(RegPar)*RegMatrix,[1 1 secDim]); 
           else
               tmpRegMatrix = RegMatrix;
               RegMatrix = sqrt(RegPar(1))*RegMatrix;
               for rp = 2 : length(RegPar)
                    RegMatrix = cat(3,RegMatrix,sqrt(RegPar(rp))*tmpRegMatrix); 
               end
           end
           parameterMatrix = cat(1,parameterMatrix,RegMatrix);
        else
           RegMatrix = sqrt(RegPar)*RegMatrix*inputParams; % removed a devision
           parameterMatrix =  RegMatrix;
        end

    end
    if jacobian                                                             % Parameter matrix is jacobian
        pars = fields(parametrizations);                                    % Get parameter names
        numberOfParameters = [];                                            % Initialize number of parameters
        for ff = 1 : length(pars)                                           % Loop over parameters
            numberOfParameters(end+1) = parametrizations.(pars{ff}).end-parametrizations.(pars{ff}).start + 1; % Calculate number of parameters
            if strcmp(parametrizations.(pars{ff}).type,'free')
                numberOfParameters(end) = numberOfParameters(end)/secDim;
            end    
        end                                                                 % End loop over parameters
        if strcmp(parametrizations.(parameterName).type,'free')
            nPars = (parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1)/secDim;
        else if strcmp(parametrizations.(parameterName).type,'dynamic')
                % nPars = (parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1)/length(parametrizations.(parameterName).parameterNames);
                nPars = (parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1);
            else
                nPars = (parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1);
            end
        end
       
        if  nPars <= max(numberOfParameters) && strcmp(parametrizations.(parameterName).RegFun,'')  % Add zeros if jacobian is too short
             if ~strcmp(parametrizations.(parameterName).type,'dynamic')
                parameterMatrix = cat(1,parameterMatrix,zeros(max(numberOfParameters),nPars,secDim));    %Add correct number of zeros to the end
             else
                parameterMatrix = cat(1,parameterMatrix,zeros(max(numberOfParameters),secDim,nPars));    %Add correct number of zeros to the end
             end
        end
    end
end

function dYdX = updateAccordingToGrouping(term, dYdX,parameterName, parametrizations,x)
% This function updates jacobians or panilty terms according to different parameter
% groupings
%
%   USAGE:
%       dYdX = updateAccordingToGrouping(term, dYdX,parameterName, parametrizations,x)
%
%   INPUTS:
%       term             = type of input (penalty, penaltyJacobian, Jacobian)
%       dYdX             = input matrix
%       parameterName    = name of parameter
%       parametrizations = paramterization options
%       x                = parameter vector
%
%   OUTPUTS:
%       dYdX             = re-organized output matrix
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
%% Update matrix lines
pars = fields(parametrizations);    % Get parameter names
run_idx = 1;                        % Setup run index 
switch term
    case 'penalty'                           
        if size(dYdX,2)~=1          % for 1D cases
            dYdXtemp = zeros(length(x),size(dYdX,2));   % Setup empty output matrix            
             for ff = 1 : length(pars)                                               % Loop over parameters
                par_idx = [parametrizations.(pars{ff}).start : parametrizations.(pars{ff}).end];    
                idx     = [1 : length(par_idx)];
                if ~isempty(parametrizations.(pars{ff}).gr)
                    if isempty(parametrizations.(pars{ff}).gr.idx_repar)
                        idx = parametrizations.(pars{ff}).gr.idx;
                    else
                        idx = parametrizations.(pars{ff}).gr.idx_repar;
                    end
                end
                if strcmp(parametrizations.(pars{ff}).type,'free') && ...
                    ~(strcmp(pars{ff},'ph0') || strcmp(pars{ff},'ph1') || strcmp(pars{ff},'gaussLB'))
                    for ll = 1 : length(idx)
                         if par_idx(1) ~= 0
                            dYdXtemp(par_idx(idx(ll)),:) = dYdXtemp(par_idx(idx(ll)),1) + dYdX(run_idx,:);
                            run_idx = run_idx + 1;
                         end
                    end
                else   
                    for ll = 1 : length(idx)
                        if par_idx(1) ~= 0         
                            dYdXtemp(par_idx(idx(ll)),:) =  dYdXtemp(par_idx(idx(ll)),:) + dYdX(run_idx,:);
                        end
                        run_idx = run_idx + 1;
                    end
                end
             end
        else
            dYdXtemp = zeros(length(x),1);         % Setup empty output matrix
             for ff = 1 : length(pars)                                               % Loop over parameters
                par_idx = [parametrizations.(pars{ff}).start : parametrizations.(pars{ff}).end]; 
                idx     = [1 : length(par_idx)];
                if ~isempty(parametrizations.(pars{ff}).gr)
                    if isempty(parametrizations.(pars{ff}).gr.idx_repar)
                        idx = parametrizations.(pars{ff}).gr.idx;
                    else
                        idx = parametrizations.(pars{ff}).gr.idx_repar;
                    end
                end            
                if strcmp(parametrizations.(pars{ff}).type,'free') && ...
                    ~(strcmp(pars{ff},'ph0') || strcmp(pars{ff},'ph1') || strcmp(pars{ff},'gaussLB'))
                    for ll = 1 : length(idx)
                         if par_idx(1) ~= 0 %&& ll == 1
                            dYdXtemp(par_idx(idx(ll)),1) = dYdXtemp(par_idx(idx(ll)),1) + dYdX(run_idx,1);
                            run_idx = run_idx + 1;
                         end
                    end
                else   
                    for ll = 1 : length(idx)
                        if par_idx(1) ~= 0         
                            dYdXtemp(par_idx(idx(ll)),1) =  dYdXtemp(par_idx(idx(ll)),1) + dYdX(run_idx,1);
                        end
                        run_idx = run_idx + 1;
                    end
                end
             end
        end
    case 'penaltyJacobian'
        if length(x) > 1
            dYdXtemp = zeros(length(x),length(x),size(dYdX,3));             % Setup empty output matrix
            for ff = 1 : length(pars)                                               % Loop over parameters
                par_idx = [parametrizations.(pars{ff}).start : parametrizations.(pars{ff}).end]; 
                idx     = [1 : length(par_idx)];
                if ~isempty(parametrizations.(pars{ff}).gr)
                    if isempty(parametrizations.(pars{ff}).gr.idx_repar)
                        idx = parametrizations.(pars{ff}).gr.idx;
                    else
                        idx = parametrizations.(pars{ff}).gr.idx_repar;
                    end
                    if strcmp(parametrizations.(pars{ff}).type,'free')
                        temp_idx = idx;
                        for sD = 2 : size(dYdX,3)
                            idx = cat(2,idx,temp_idx+ max(temp_idx)*(sD-1));
                        end
                    end
                end 
                    for ll = 1 : length(idx)
                        if par_idx(1) ~= 0         
                            dYdXtemp(par_idx(idx(ll)),par_idx(idx(ll)),:) = dYdXtemp(par_idx(idx(ll)),par_idx(idx(ll)),:) + dYdX(run_idx,run_idx,:);
                        end
                        run_idx = run_idx + 1;
                    end
             end
        else
            dYdXtemp = zeros(x,size(dYdX,2));                                   % Setup empty output matrix
            for ff = 1 : length(pars)                                               % Loop over parameters
                par_idx = [parametrizations.(pars{ff}).start : parametrizations.(pars{ff}).end]; 
                idx     = [1 : length(par_idx)];
                if ~isempty(parametrizations.(pars{ff}).gr)
                    if isempty(parametrizations.(pars{ff}).gr.idx_repar)
                        idx = parametrizations.(pars{ff}).gr.idx;
                    else
                        idx = parametrizations.(pars{ff}).gr.idx_repar;
                    end
                    if strcmp(parametrizations.(pars{ff}).type,'free')
                        temp_idx = idx;
                        for sD = 2 : size(dYdX,3)
                            idx = cat(2,idx,temp_idx+ max(temp_idx)*(sD-1));
                        end
                    end
                end 
                for ll = 1 : length(idx)
                    if par_idx(1) ~= 0         
                        dYdXtemp(par_idx(idx(ll)),:) = dYdXtemp(par_idx(idx(ll)),:) + dYdX(run_idx,:);
                    end
                    run_idx = run_idx + 1;
                end
            end
            dYdX = dYdXtemp;
            % Now we need to add up the columns (this is needed for soft
            % constraints which are not diagonal)
            dYdXtemp = zeros(x,x);                                          % Setup empty output matrix
            run_idx = 1;                                                    % Reset run index
            for ff = 1 : length(pars)                                               % Loop over parameters
                par_idx = [parametrizations.(pars{ff}).start : parametrizations.(pars{ff}).end]; 
                idx     = [1 : length(par_idx)];
                if ~isempty(parametrizations.(pars{ff}).gr)
                    if isempty(parametrizations.(pars{ff}).gr.idx_repar)
                        idx = parametrizations.(pars{ff}).gr.idx;
                    else
                        idx = parametrizations.(pars{ff}).gr.idx_repar;
                    end
                    if strcmp(parametrizations.(pars{ff}).type,'free')
                        temp_idx = idx;
                        for sD = 2 : size(dYdX,3)
                            idx = cat(2,idx,temp_idx+ max(temp_idx)*(sD-1));
                        end
                    end
                end 
                for ll = 1 : length(idx)
                    if par_idx(1) ~= 0         
                        dYdXtemp(:,par_idx(idx(ll))) = dYdXtemp(:,par_idx(idx(ll))) + dYdX(:,run_idx);
                    end
                    run_idx = run_idx + 1;
                end
            end
        end
                                
    case 'Jacobian'
        if strcmp(parametrizations.(parameterName).type,'free')            
            dYdXtemp = zeros(size(dYdX,1),(parametrizations.(parameterName).end - parametrizations.(parameterName).start + 1)/size(dYdX,3),size(dYdX,3));
            if ~isempty(parametrizations.(parameterName).gr)
                if isempty(parametrizations.(parameterName).gr.idx_repar)
                    idx = parametrizations.(parameterName).gr.idx;
                else
                    idx = parametrizations.(parameterName).gr.idx_repar;
                end
            end
            for ll = 1 : length(idx)  
                dYdXtemp(:,idx(ll),:) = dYdXtemp(:,idx(ll),:) +  dYdX(:,run_idx,:);
                run_idx = run_idx + 1;
            end
        end
        if strcmp(parametrizations.(parameterName).type,'fixed')
            idx     = [1 : parametrizations.(parameterName).end - parametrizations.(parameterName).start + 1];
            dYdXtemp = zeros(size(dYdX,1),length(idx),size(dYdX,3));
            if ~isempty(parametrizations.(parameterName).gr)
                if isempty(parametrizations.(parameterName).gr.idx_repar)
                    idx = parametrizations.(parameterName).gr.idx;
                else
                    idx = parametrizations.(parameterName).gr.idx_repar;
                end
            end
            for ll = 1 : length(idx)  
                dYdXtemp(:,idx(ll),:) = dYdXtemp(:,idx(ll),:) +  dYdX(:,run_idx,:);
                run_idx = run_idx + 1;
            end
        end
        if strcmp(parametrizations.(parameterName).type,'dynamic')
            idx     = [1 : parametrizations.(parameterName).end - parametrizations.(parameterName).start + 1];
            if ndims(dYdX)==3
                dYdXtemp = zeros(size(dYdX,1),size(dYdX,2),length(idx));
                if ~isempty(parametrizations.(parameterName).gr)
                    if isempty(parametrizations.(parameterName).gr.idx_repar)
                        idx = parametrizations.(parameterName).gr.idx;
                    else
                        idx = parametrizations.(parameterName).gr.idx_repar;
                    end
                end
                for ll = 1 : length(idx)  
                    dYdXtemp(:,:,idx(ll)) = dYdXtemp(:,:,idx(ll)) +  dYdX(:,:,run_idx);
                    run_idx = run_idx + 1;
                end
            else
                dYdXtemp = zeros(size(dYdX,1),length(idx));
                if ~isempty(parametrizations.(parameterName).gr)
                    if isempty(parametrizations.(parameterName).gr.idx_repar)
                        idx = parametrizations.(parameterName).gr.idx;
                    else
                        idx = parametrizations.(parameterName).gr.idx_repar;
                    end
                end
                for ll = 1 : length(idx)  
                    dYdXtemp(:,idx(ll)) = dYdXtemp(:,idx(ll)) +  dYdX(:,run_idx);
                    run_idx = run_idx + 1;
                end
            end
        end        
    end
dYdX = dYdXtemp;
end




function dYdX = updateJacobianBlock(dYdX,parameterName, parametrizations,inputParams,Domain,NoiseSD,dYdpenExpect,fitRange,dYdpenSoftCon,Reg)
% This function updates jacobians according to the parametrization
%
%   USAGE:
%       dYdX = updateJacobianBlock(dYdX,parameterName, parametrizations,inputParams,Domain,NoiseSD)
%
%   INPUTS:
%       dYdX             = jacobian
%       parameterName    = name of parameter
%       parametrizations = paramterization options
%       inputParams      = parameter values
%       Domain           = signal domain
%       NoiseSD          = standard deviation of the noise
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

    dYdpenSoftCon  = updateAccordingToGrouping('penaltyJacobian',dYdpenSoftCon,parameterName, parametrizations,size(dYdpenExpect,1));
    dYdpenSoftConToAdd = zeros([size(dYdpenSoftCon),secDim]);
    dYdpenSoftConToAdd(:,:,1) = dYdpenSoftCon;

    % We want to add the regularizer terms here
    if strcmp(parametrizations.(parameterName).type,'free') || ...
       strcmp(parametrizations.(parameterName).type,'fixed')  
        if ~isempty(parametrizations.(parameterName).gr)           
            dYdX  = updateAccordingToGrouping('Jacobian',dYdX,parameterName, parametrizations,[]);                    
        end
        if Reg                                                                  % Add parameter regularization
            [dYdX]        = addParameterRegularization(dYdX,parameterName, parametrizations,inputParams.(parameterName),1,secDim,[]);
        end  
        nLines = size(dYdX,2);
    end

    if ~isempty(NoiseSD)
        nPoints = size(dYdX,1);                                             % Update the number of points
        if strcmp(Domain,'FD')
            SigmaValue = NoiseSD.FD;                                             % Get sigma
        end
        if strcmp(Domain,'TD')
            SigmaValue = NoiseSD.TD;                                             % Get sigma
        end
        if strcmp(Domain,'FDTD')
            indMinTD = fitRange.TD(1);                                            % Get time domain indices
            indMaxTD = fitRange.TD(2);                                            % Get time domain indices
            SigmaFD = NoiseSD.FD;                                           % Get sigma squared
            SigmaFD = repmat(SigmaFD, [size(dYdX,1)-size(dYdX(indMinTD:indMaxTD,:,:),1) 1]);            % Repeat according to dimensions
            SigmaTD = NoiseSD.TD;                                           % Get sigma squared
            SigmaTD = repmat(SigmaTD, [size(dYdX(indMinTD:indMaxTD,:,:),1) 1]);            % Repeat according to dimensions
            SigmaValue  = cat(1,SigmaFD,SigmaTD);                                % Combine frequency and time domain     
        end
        Sigma = repmat(SigmaValue', [1 nPoints nLines]);                         % Repeat according to dimensions
        Sigma = permute(Sigma,[2 3 1]);                                     % Dims have to be nPoints nLines secDim
        Sigma = squeeze(Sigma);                                             % Remove zero dimensions
        dYdX = squeeze(dYdX);                                               % Remove zero dimensions
        dYdX(:,:,:) = dYdX(:,:,:) ./ Sigma;                 % Normalize jacobian
    end
    

    % Free parametrizations need to add secDim copies to the jacobian and
    % set partial derivatives to zero for e.g. df1/dph02 .
    if strcmp(parametrizations.(parameterName).type,'free')
        nPoints = size(dYdX,1);                                             % Update the number of points
        dYdX = squeeze(dYdX);                                               %Remove length 1 dims (e.g. for ph0, ph1, gaussLB)
        dYdX = repmat(dYdX,[1 1 1 secDim]);                                 %Create secDim copies
        dYdX = squeeze(dYdX);                                               %Remove length 1 dims (needed for 3D dYdX case)
        factor = repmat(eye(secDim),[1 1 nPoints nLines]);                  %To delete partial derivatives not on diagonal
        factor = permute(factor,[3 4 1 2]);                                 %Reorder to match dYdX dimensions
        factor = squeeze(factor);                                           %Remove length 1 dims
        dYdX(:,:,:,:) = dYdX(:,:,:,:) .* factor;            %Delete partial derivatives not on diagonal     
    end
    % Fixed parametrization needs to concatenate along secDim resulting in
    % a single line in the jacobian
    if strcmp(parametrizations.(parameterName).type,'fixed')
        dYdX = squeeze(dYdX);                                               % Remove length 1 dims (e.g. for ph0, ph1, gaussLB)
        nPoints = size(dYdX,1);                                             % Update the number of points
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
        factor = parametrizations.(parameterName).fun.jac(parameterEstimate,parametrizations.(parameterName).modulator);
        if ndims(factor) ==3                                               % For per metabolite cases with more than 1 entry       
            factor = repmat(factor, [1 1 1 nPoints]);                      % Repeat nPoints times
            factor = permute(factor,[4 2 1 3]);                            % Dims have to be nPoints secDim nLines nPars
        else
            factor = repmat(factor, [1 1 nPoints]);                      % Repeat nPoints times
            factor = permute(factor,[3 1 2]);                            % Dims have to be nPoints secDim nLines
        end
        dYdXOrginal = dYdX;                                                % Backup original derivatives
        % Multiply the original derivatives with the derivatives from the
        % reparametrization
        dYdX = [];
        % Loop over new parameters
        nPars = length(parametrizations.(parameterName).parameterNames);    % Parameters for indirect dynamic function
        for rp = 1 : nPars
            if ndims(factor) == 4                                               % For per metabolite cases with more than 1 entry 
                dYdX = cat(4,dYdX,dYdXOrginal .* factor(:,:,:,rp));
            else
                dYdX = cat(3,dYdX,dYdXOrginal .* factor(:,:,rp));
            end
        end 
        dYdX = squeeze(dYdX);                                               %Remove length 1 dims
        if ndims(dYdX) ==3                                                  
            secDim = size(dYdX,2);
        else
            secDim = size(dYdX,3);
        end

    end
    
    % We want to add the expectation value soft constraint regularizer here including
    % grouping of different variables
    dYdpenExpect = dYdpenExpect(:,parametrizations.(parameterName).start:parametrizations.(parameterName).end,:);
    dYdpenSoftConToAdd = dYdpenSoftConToAdd(:,parametrizations.(parameterName).start:parametrizations.(parameterName).end,:);
      

    % Finally we have to concatenate along the indirect dimension
    switch ndims(dYdX)
        case 2   % For fixed parametrizations of ph0, ph1, gaussLB or 1D data parameter per basis function
            dYdX  = cat(1,dYdX,-squeeze(dYdpenExpect)) ;                % Add expectation value regularizer to for global parameter
            dYdX  = cat(1,dYdX,-squeeze(dYdpenSoftConToAdd)) ;          % Add soft constraint value regularizer to for global parameter
            if secDim > 1
                dYdX = reshape(dYdX,[],1);
            end
        case 3   % E.g. Fixed parameter or single basis function case
            if strcmp(parametrizations.(parameterName).type,'fixed')
                dYdX = permute(dYdX,[1 3 2]);             
            end
            if strcmp(parametrizations.(parameterName).type,'free')
                if ~(strcmp(parameterName,'baseAmpl') || strcmp(parameterName,'metAmpl') || strcmp(parameterName,'freqShift') || strcmp(parameterName,'lorentzLB'))
                    nLines = secDim;
                else
                    if nLines == 1          % This is needed for single basis cases
                        nLines = secDim;
                    end
                end
            end
            if strcmp(parametrizations.(parameterName).type,'dynamic')
                nLines = nLines*nPars;
            end
            dYdpenExpect = permute(dYdpenExpect,[1 3 2]);
            dYdpenSoftConToAdd = permute(dYdpenSoftConToAdd,[1 3 2]);
            dYdX  = cat(1,dYdX,-dYdpenExpect);                % Add regularizer to parameter per basis function
            dYdX  = cat(1,dYdX,-dYdpenSoftConToAdd);          % Add regularizer to parameter per basis function
            dYdX = reshape(dYdX,[],nLines);
        case  4 % E.g. free parametrization of per metabolite parameters or dynamic cases
            if strcmp(parametrizations.(parameterName).type,'free')
                dYdX = permute(dYdX,[1 3 4 2]);
                dYdX = reshape(dYdX,nPoints,secDim,secDim*nLines);
                dYdpenExpect = reshape(permute(dYdpenExpect,[1 3 2]),size(dYdpenExpect,1),secDim,nLines,secDim);
                dYdpenExpect = permute(dYdpenExpect,[1 2 4 3]);
                dYdpenExpect = reshape(dYdpenExpect,size(dYdpenExpect,1),secDim,nLines*secDim);
                dYdpenSoftConToAdd = permute(dYdpenSoftConToAdd,[1 3 2]);
                dYdX  = cat(1,dYdX,-dYdpenExpect);
                dYdX  = cat(1,dYdX,-dYdpenSoftConToAdd);
                dYdX = reshape(dYdX,[],secDim*nLines);
            end
            if strcmp(parametrizations.(parameterName).type,'dynamic')
                dYdX = permute(dYdX,[1 3 4 2]);
                dYdX = reshape(dYdX,[],secDim,nPars*nLines); 
               if ~isempty(parametrizations.(parameterName).gr)
                    [dYdX]  = updateAccordingToGrouping('Jacobian',dYdX,parameterName, parametrizations,[]);                    
               end
               if Reg                                                                  % Add parameter regularization
                    [dYdX]        = addParameterRegularization(dYdX,parameterName, parametrizations,inputParams.(parameterName),1,secDim,[]);
                     Sigma = repmat(SigmaValue', [1 ((size(dYdX,1)-nPoints)) size(dYdX,3)]);  % Repeat according to dimensions
                     Sigma = permute(Sigma,[2 1 3]);                                     % Dims have to be nPoints nLines secDim
                     Sigma = squeeze(Sigma);                                             % Remove zero dimensions
                     dYdX(nPoints+1:end,:,:) = dYdX(nPoints+1:end,:,:) ./ Sigma;                 % Normalize jacobian
               end  
               dYdpenExpect = permute(dYdpenExpect,[1 3 2]);
               dYdpenSoftConToAdd = permute(dYdpenSoftConToAdd,[1 3 2]);
               dYdX  = cat(1,dYdX,-dYdpenExpect);                % Add regularizer to parameter per basis function
               dYdX  = cat(1,dYdX,-dYdpenSoftConToAdd);          % Add regularizer to parameter per basis function
               dYdX = reshape(dYdX,[],size(dYdX,3)); 
            end
   end
end

function [penaltyTerm] = calcPenalty(jacobian,inputParams, parametrizations, secDim);
param = fields(parametrizations);                                  % Loop over all model parameters

for sD = 1 : secDim                                                 % Loop over indirect dimension
    tempTerm = [];                                                      % Initialize empty vector to collect the penalties for each parameter            
    for pp = 1:length(param)
        
        % Calculates the penalty term for deviations from expectation values
        % relative to standard deviations       
        if strcmp(parametrizations.(param{pp}).type,'free')
            expectationValue    = parametrizations.(param{pp}).ex(sD,:);
            standardDeviation   = parametrizations.(param{pp}).sd(sD,:);
            actualValue         = squeeze(inputParams.(param{pp})(sD,:));
        end
        if strcmp(parametrizations.(param{pp}).type,'fixed')
            expectationValue    = parametrizations.(param{pp}).ex(1,:);
            standardDeviation   = parametrizations.(param{pp}).sd(1,:);
            actualValue         = squeeze(inputParams.(param{pp})(1,:));
        end
        if strcmp(parametrizations.(param{pp}).type,'dynamic')
            actualValue =[];
            expectationValue =[];
            standardDeviation =[];
            for dd = 1 : length(parametrizations.(param{pp}).parameterNames)
                actualValue         = [actualValue squeeze(inputParams.([param{pp} 'Reparametrization']).(parametrizations.(param{pp}).parameterNames{dd})(1,:))];
                expectationValue    = [expectationValue parametrizations.(param{pp}).ex(dd,:)];
                standardDeviation   = [standardDeviation parametrizations.(param{pp}).sd(dd,:)];
            end     
            nan_vector=isnan(actualValue);
            actualValue(nan_vector)=[];
            expectationValue(nan_vector)=[];
            standardDeviation(nan_vector)=[];
        end
        
        
        % In the LCModel objective functions, the terms are in the shape of
        % (difference-to-expectation-value)^2/(standard-deviation)^2. We are
        % here formulating the vector of penalties that is appended to the
        % residual, i.e., squared afterwards. We therefore give the penalty vector
        % in the shape of (difference-to-expectation)/(standard-deviation).
        diffVec     = actualValue-expectationValue;
        % Scale with regards to indirect dimension
        if strcmp(parametrizations.(param{pp}).type,'fixed') ||...
           strcmp(parametrizations.(param{pp}).type,'dynamic') 
           diffVec = diffVec/secDim; 
        end
        sdVec       = standardDeviation;
        penalty     = diffVec./sdVec;
        
        if ~jacobian                             
            if ~strcmp(parametrizations.(param{pp}).type,'none')                % The baseline parameter might be empty
                if ~strcmp(parametrizations.(param{pp}).type,'free')
                    tempTerm = cat(2,tempTerm,penalty);                         % Add penalty term
                else
                    for sD_free = 1 : secDim 
                        if sD_free == sD
                            tempTerm = cat(2,tempTerm,penalty);                 % Add penalty term
                        else
                            tempTerm = cat(2,tempTerm,zeros(size(penalty)));    % Add penalty term
                        end
                    end
                end
            end
        else
            % Jacobian
            if strcmp(parametrizations.(param{pp}).type,'fixed') ||...
               strcmp(parametrizations.(param{pp}).type,'dynamic') 
                penalty  = 1./(sdVec*secDim);                % Because the penalty is a linear sum of its parts
            else
                penalty  = 1./sdVec;                % Because the penalty is a linear sum of its parts
            end                           
            if ~strcmp(parametrizations.(param{pp}).type,'none')            % The baseline parameter might be empty
                if ~strcmp(parametrizations.(param{pp}).type,'free')
                    tempTerm = cat(2,tempTerm,penalty);                         % Add penalty term
                else
                    penalty = repmat(penalty,[1 secDim]);
                    factor = zeros(size(penalty));
                    factor(length(diffVec)*(sD-1)+1:length(diffVec)*(sD-1)+length(diffVec)) = 1;
                    penalty = penalty .* factor;
                    tempTerm = cat(2,tempTerm,penalty);
                end
            end
        end

    end
    if ~jacobian 
        penaltyTerm(:,sD) = tempTerm;                                        % Save as one concatenated vector of penalties to be appended to the residual
    else
        tempTerm = diag(tempTerm); 
        penaltyTerm(:,:,sD) = tempTerm;                                        % Save as one concatenated vector of penalties to be appended to the residual
    end
end  
end

function [penalty, penaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, param, sD, nParams)
penalty = [];
penaltyJac = [];
if ~isempty(parametrizations.(param).sc)
    if nargout == 1
        if strcmp(parametrizations.(param).type,'free')
            penalty = zeros(1,nParams*size(inputParams.(param),1));
        end
        if strcmp(parametrizations.(param).type,'fixed')
            penalty = zeros(1,nParams);
        end
        if strcmp(parametrizations.(param).type,'dynamic')
            n_repar = length(parametrizations.(param).parameterNames);
            penalty  = zeros(1,n_repar*length(squeeze(inputParams.(param)(sD,:))));
        end
    else
        penaltyJac  = zeros(nParams,nParams);
    end

    % if ~strcmp(parametrizations.(param).type,'dynamic')
    %     penalty  = zeros(1,length(squeeze(inputParams.(param)(sD,:))));
    % else    
    %     n_repar = length(parametrizations.(param).parameterNames);
    %     penalty  = zeros(1,n_repar*length(squeeze(inputParams.(param)(sD,:))));
    % end
    % penaltyJac  = zeros(nParams,nParams);
    for qq = 1 : length(parametrizations.(param).sc.fix_idx)
        fixValue=[];
        adjValue=[];
        % Calculates the penalty term for deviations from soft constraints
        for mm = 1 : size(parametrizations.(param).sc.fix_idx{qq},1)
            idx           = parametrizations.(param).sc.fix_idx{qq}(mm,:);
            idx(idx==0)   =[];
            tempfixValue  = squeeze(inputParams.(param)(sD,idx));
            fixValue(mm)  = sum(tempfixValue .* parametrizations.(param).sc.fix_factor{qq}{mm}');
        end
        for mm = 1 : size(parametrizations.(param).sc.adj_idx{qq},1)
            idx           = parametrizations.(param).sc.adj_idx{qq}(mm,:);
            idx(idx==0)   =[];
            tempadjValue  = squeeze(inputParams.(param)(sD,idx));
            adjValue(mm)  = sum(tempadjValue);
        end
        if strcmp(parametrizations.(param).sc.fun,'ratio')
            actualValue = adjValue ./fixValue;
            actualValue(isnan(actualValue)) = 0;
            actualValue(isinf(actualValue)) = 0;
        else
            actualValue = fixValue - adjValue;
        end
        
        expectationValue = parametrizations.(param).sc.ex{qq}{:}';
        standardDeviation   = parametrizations.(param).sc.sd{qq}{:}';
    
    
        diffVec     = actualValue-expectationValue;
        sdVec       = standardDeviation;
    
        if nargout == 1
            if strcmp(parametrizations.(param).type,'free')   
                factor=(reshape(repmat([1:size(inputParams.(param),1)]',[1 length(parametrizations.(param).sc.adj_idx{qq})])',1,[]))-1;
                factor = factor * nParams;
                adj_idx = repmat(parametrizations.(param).sc.adj_idx{qq}',[1 size(inputParams.(param),1)]);
                penalty(adj_idx+factor) = repmat(parametrizations.(param).sc.scaling*diffVec./sdVec,[1 size(inputParams.(param),1)]);
            end
            if strcmp(parametrizations.(param).type,'fixed')
                penalty(parametrizations.(param).sc.adj_idx{qq}) = parametrizations.(param).sc.scaling*diffVec./sdVec;
            end
            if strcmp(parametrizations.(param).type,'dynamic')
                n_repar = length(parametrizations.(param).parameterNames);
                penalty(parametrizations.(param).sc.adj_idx{qq}*n_repar-(n_repar-qq)) =  parametrizations.(param).sc.scaling*diffVec./sdVec;
            end
        end

        if strcmp(parametrizations.(param).type,'dynamic')
            n_repar = length(parametrizations.(param).parameterNames);
        end
          
        % Jacobian    
        if nargout ~=1            
            if sum(fixValue==0)==0 
                if strcmp(parametrizations.(param).type,'free')  
                    sdVec       = standardDeviation;
                    sdVec = repmat(sdVec,[1 size(inputParams.(param),1)]);
                end
                for m_adj = 1 : size(parametrizations.(param).sc.adj_idx{qq},1)
                    if fixValue(m_adj) ~= 0
                        fix_idx           = parametrizations.(param).sc.fix_idx{qq}(m_adj,:);
                        fix_idx(fix_idx==0)   =[];
                        adj_idx           = parametrizations.(param).sc.adj_idx{qq}(m_adj,:);
                        adj_idx(adj_idx==0)   =[];
                        if strcmp(parametrizations.(param).type,'dynamic')
                            fix_idx = n_repar*fix_idx+qq-n_repar;
                            adj_idx = n_repar*adj_idx+qq-n_repar;
                        end
                        for m_fix = 1 : length(fix_idx)
                            penaltyJac(parametrizations.(param).startNoGroup+adj_idx-1,parametrizations.(param).startNoGroup+fix_idx(m_fix)-1) = -parametrizations.(param).sc.scaling*parametrizations.(param).sc.fix_factor{qq}{m_adj}(m_fix)*adjValue(m_adj)./(sdVec(m_adj).*fixValue(m_adj).*fixValue(m_adj));
                        end
                        penaltyJac(parametrizations.(param).startNoGroup+adj_idx-1,parametrizations.(param).startNoGroup+adj_idx-1) = parametrizations.(param).sc.scaling*1./(sdVec(m_adj).*fixValue(m_adj)); 
                    end
                end
            end
        end
    end

else

    if strcmp(parametrizations.(param).type,'free')
        penalty = zeros(1,nParams*size(inputParams.(param),1));
    end
    if strcmp(parametrizations.(param).type,'fixed')
        penalty = zeros(1,nParams);
    end
    if strcmp(parametrizations.(param).type,'dynamic')
        penalty = zeros(1,nParams*length(parametrizations.(param).parameterNames));
    end
    penaltyJac  = zeros(nParams,nParams);    
end
end
