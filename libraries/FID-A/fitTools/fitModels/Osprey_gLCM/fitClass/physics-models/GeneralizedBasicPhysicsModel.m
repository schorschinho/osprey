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

    residual  = cat(1,residual,-regu);                                      % Add regularizer term
    residual  = cat(1,residual,penaltyTerm);                                % Add penalty term
    residual  = cat(1,residual,penaltyTermSoftConstraint);                  % Add soft constraint penalty term
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


    [dYdlorentzLB]  = updateAccordingToGrouping(dYdlorentzLB,'lorentzLB', parametrizations); % Update jacobian block for lorentzLB parameter
    [dYdfreqShift]  = updateAccordingToGrouping(dYdfreqShift,'freqShift', parametrizations); % Update jacobian block for freqShift parameter
    [dYdmetAmpl]    = updateAccordingToGrouping(dYdmetAmpl,'metAmpl', parametrizations); % Update jacobian block for metAmpl parameter

    if Reg                                                                  % Add parameter regularization
        [dYdph0]        = addParameterRegularization(dYdph0,'ph0', parametrizations,ph0,1,secDim); % Add regularizer to ph0 parameter
        [dYdph1]        = addParameterRegularization(dYdph1,'ph1', parametrizations,ph1,1,secDim); % Add regularizer to ph1 parameter
        [dYdgaussLB]    = addParameterRegularization(dYdgaussLB,'gaussLB', parametrizations,gaussLB,1,secDim); % Add regularizer to gaussLB parameter
        [dYdlorentzLB]  = addParameterRegularization(dYdlorentzLB,'lorentzLB', parametrizations,lorentzLB,1,secDim); % Add regularizer to lorentzLB parameter
        [dYdfreqShift]  = addParameterRegularization(dYdfreqShift,'freqShift', parametrizations,freqShift,1,secDim); % Add regularizer to freqShift parameter
        [dYdmetAmpl]    = addParameterRegularization(dYdmetAmpl,'metAmpl', parametrizations,metAmpl,1,secDim); % Add regularizer to metAmpl parameter
        if nBaselineComps ~= 0
            [dYdbaseAmpl]   = addParameterRegularization(dYdbaseAmpl,'baseAmpl', parametrizations,baseAmpl,1,secDim); % Add regularizer to baseAmpl parameter
        end
    end                                                                     % End loop over indirect dimension

    
    
    if secDim == 1 && ~isempty(NoiseSD)                                     % 1D jacobians have not been normalized yet
        if strcmp(Domain,'FD')
            Sigma = NoiseSD.FD;                                             % Get sigma
        end
        if strcmp(Domain,'TD')
            Sigma = NoiseSD.TD;                                             % Get sigma
        end
        if strcmp(Domain,'FDTD')
            SigmaFD = NoiseSD.FD;                                           % Get sigma squared
            SigmaFD = repmat(SigmaFD, [size(dYdph0,1)-size(dYdph0TD(indMinTD:indMaxTD,:,:),1) 1]);            % Repeat according to dimensions
            SigmaTD = NoiseSD.TD;                                           % Get sigma squared
            SigmaTD = repmat(SigmaTD, [size(dYdph0TD(indMinTD:indMaxTD,:,:),1) 1]);            % Repeat according to dimensions
            Sigma  = cat(1,SigmaFD,SigmaTD);                                % Combine frequency and time domain     
        end
        dYdph0          = dYdph0 ./ Sigma;                            % Cut out fit range
        dYdph1          = dYdph1 ./ Sigma;                            % Cut out fit range
        dYdgaussLB      = dYdgaussLB ./ Sigma;                        % Cut out fit range
        dYdlorentzLB    = dYdlorentzLB ./ Sigma;                    % Cut out fit range
        dYdfreqShift    = dYdfreqShift ./ Sigma;                    % Cut out fit range
        dYdmetAmpl      = dYdmetAmpl ./ Sigma;                      % Cut out fit range
        if nBaselineComps ~= 0                                                  % Has baseline?
            dYdbaseAmpl           = dYdbaseAmpl ./ Sigma;           % Cut out fit range
        end
    end

    if secDim > 1                                                           % update each block in the jacobian according to the 2-D parametrization
        [dYdph0]        = updateJacobianBlock(dYdph0,'ph0', parametrizations,inputParams,SignalPart,NoiseSD); % Update jacobian block for ph0 parameter
        [dYdph1]        = updateJacobianBlock(dYdph1,'ph1', parametrizations,inputParams,SignalPart,NoiseSD); % Update jacobian block for ph1 parameter
        [dYdgaussLB]    = updateJacobianBlock(dYdgaussLB,'gaussLB', parametrizations,inputParams,SignalPart,NoiseSD); % Update jacobian block for gaussLB parameter
        [dYdlorentzLB]  = updateJacobianBlock(dYdlorentzLB,'lorentzLB', parametrizations,inputParams,SignalPart,NoiseSD); % Update jacobian block for lorentzLB parameter
        [dYdfreqShift]  = updateJacobianBlock(dYdfreqShift,'freqShift', parametrizations,inputParams,SignalPart,NoiseSD); % Update jacobian block for freqShift parameter
        [dYdmetAmpl]    = updateJacobianBlock(dYdmetAmpl,'metAmpl', parametrizations,inputParams,SignalPart,NoiseSD); % Update jacobian block for metAmpl parameter
        if nBaselineComps ~= 0
            [dYdbaseAmpl]   = updateJacobianBlock(dYdbaseAmpl,'baseAmpl', parametrizations,inputParams,SignalPart,NoiseSD);  % Update jacobian block for baseAmpl parameter
        end
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

    if secDim == 1                                                          % Add penalty term for parameters deviating from expectation values
        tempTerm = [];

        [~,ph0PenaltyJac] = calcPenalty(inputParams, parametrizations, 'ph0', sD);
        [~,ph1PenaltyJac] = calcPenalty(inputParams, parametrizations, 'ph1', sD);
        [~,gaussLBPenaltyJac] = calcPenalty(inputParams, parametrizations, 'gaussLB', sD);
        [~,lorentzLBPenaltyJac] = calcPenalty(inputParams, parametrizations, 'lorentzLB', sD);
        [~,freqShiftPenaltyJac] = calcPenalty(inputParams, parametrizations, 'freqShift', sD);
        [~,metAmplPenaltyJac] = calcPenalty(inputParams, parametrizations, 'metAmpl', sD);
        if nBaselineComps ~= 0
            [~,baseAmplPenaltyJac] = calcPenalty(inputParams, parametrizations, 'baseAmpl', sD);
        else
            baseAmplPenaltyJac = [];
        end
        tempTerm = cat(2,tempTerm,[ph0PenaltyJac , ph1PenaltyJac, gaussLBPenaltyJac, lorentzLBPenaltyJac, freqShiftPenaltyJac, metAmplPenaltyJac, baseAmplPenaltyJac]); % Concatenate penalty terms

        dYdpen = diag(tempTerm);                                            % Copy to diagonal matrix (derivatives only affect each parameter itself, others are zero!)

        jac = cat(1, jac, dYdpen);                                          % Append to Jacobian

    end
    
    if secDim == 1                                                          % Add penalty term for parameters deviating from expectation values
        nParams = length(x);
        tempTerm = [];
        [~,ph0PenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'ph0', sD, nParams);
        [~,ph1PenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'ph1', sD, nParams);
        [~,gaussLBPenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'gaussLB', sD, nParams);
        [~,lorentzLBPenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'lorentzLB', sD, nParams);
        [~,freqShiftPenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'freqShift', sD, nParams);
        [~,metAmplPenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'metAmpl', sD, nParams);
        if nBaselineComps ~= 0 
            [~,baseAmplPenaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, 'baseAmpl', sD, nParams);
        else
            baseAmplPenaltyJac = zeros(size(ph0PenaltyJac));
        end
        tempTerm = ph0PenaltyJac + ph1PenaltyJac + gaussLBPenaltyJac + lorentzLBPenaltyJac + freqShiftPenaltyJac + metAmplPenaltyJac + baseAmplPenaltyJac; % Concatenate penalty terms

        jac = cat(1, jac, tempTerm);                                          % Append to Jacobian

    end

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
            [ph0Reg]        = addParameterRegularization([],'ph0', parametrizations,ph0,0,secDim); % Calculate regularizer for ph0 parameter
            [ph1Reg]        = addParameterRegularization([],'ph1', parametrizations,ph1,0,secDim); % Calculate regularizer for ph1 parameter
            [gaussLBReg]    = addParameterRegularization([],'gaussLB', parametrizations,gaussLB,0,secDim); % Calculate regularizer for gaussLB parameter
            [lorentzLBReg]  = addParameterRegularization([],'lorentzLB', parametrizations,lorentzLB,0,secDim); % Calculate regularizer for lorentzLB parameter
            [freqShiftReg]  = addParameterRegularization([],'freqShift', parametrizations,freqShift,0,secDim); % Calculate regularizer for freqShift parameter
            [metAmplReg]    = addParameterRegularization([],'metAmpl', parametrizations,metAmpl,0,secDim); % Calculate regularizer for metAmpl parameter
            [baseAmplReg]   = addParameterRegularization([],'baseAmpl', parametrizations,baseAmpl',0,secDim); % Calculate regularizer for baseAmpl parameter
            regu = cat(2,regu,[ph0Reg , ph1Reg, gaussLBReg, lorentzLBReg, freqShiftReg, metAmplReg, baseAmplReg]); % Concatenate regularizer
        end
    end                                                                     % End loop over indirect dimension

    if secDim == 1

        ph0Penalty = calcPenalty(inputParams, parametrizations, 'ph0', sD);
        ph1Penalty = calcPenalty(inputParams, parametrizations, 'ph1', sD);
        gaussLBPenalty = calcPenalty(inputParams, parametrizations, 'gaussLB', sD);
        lorentzLBPenalty = calcPenalty(inputParams, parametrizations, 'lorentzLB', sD);
        freqShiftPenalty = calcPenalty(inputParams, parametrizations, 'freqShift', sD);
        metAmplPenalty = calcPenalty(inputParams, parametrizations, 'metAmpl', sD);
        if ~isempty(bl) 
            baseAmplPenalty = calcPenalty(inputParams, parametrizations, 'baseAmpl', sD);
        else
            baseAmplPenalty = [];
        end
        penaltyTerm = cat(2,penaltyTerm,[ph0Penalty , ph1Penalty, gaussLBPenalty, lorentzLBPenalty, freqShiftPenalty, metAmplPenalty, baseAmplPenalty]); % Concatenate penalty terms
        penaltyTerm = penaltyTerm';
    end

    if secDim == 1
        nParams = length(x);
        ph0SoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'ph0', sD, nParams);
        ph1SoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'ph1', sD, nParams);
        gaussLBSoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'gaussLB', sD, nParams);
        lorentzLBSoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'lorentzLB', sD, nParams);
        freqShiftSoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'freqShift', sD, nParams);
        metAmplSoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'metAmpl', sD, nParams);
        if ~isempty(bl) 
            baseAmplSoftContstraint = calcSoftConstraintPenalty(inputParams, parametrizations, 'baseAmpl', sD, nParams);
        else
            baseAmplSoftContstraint = [];
        end
        penaltyTermSoftConstraint = cat(2,penaltyTermSoftConstraint,[ph0SoftContstraint , ph1SoftContstraint, gaussLBSoftContstraint, lorentzLBSoftContstraint, freqShiftSoftContstraint, metAmplSoftContstraint, baseAmplSoftContstraint]); % Concatenate penalty terms
        penaltyTermSoftConstraint = penaltyTermSoftConstraint';
    end    
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
                        paramStruct.(pars{ff}) = paramStruct.(pars{ff})(parametrizations.(pars{ff}).gr.idx); 
                    end
                end
                if strcmp(parametrizations.(pars{ff}).type,'fixed')
                    if ~isempty(parametrizations.(pars{ff}).gr)  
                        paramStruct.(pars{ff}) = paramStruct.(pars{ff})(parametrizations.(pars{ff}).gr.idx); 
                    end
                    paramStruct.(pars{ff}) = reshape(paramStruct.(pars{ff}),1,[]);
                    paramStruct.(pars{ff}) = repmat(paramStruct.(pars{ff}),[secDim,1]);
                end
                if strcmp(parametrizations.(pars{ff}).type,'dynamic')
                    if ~isempty(parametrizations.(pars{ff}).gr)  
                        paramStruct.(pars{ff}) = paramStruct.(pars{ff})(parametrizations.(pars{ff}).gr.idx); 
                    end
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
function parameterMatrix = addParameterRegularization(parameterMatrix,parameterName, parametrizations,inputParams,jacobian,secDim)
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
           RegMatrix = repmat(sqrt(parametrizations.(parameterName).RegPar)*RegMatrix,[1 1 secDim])/secDim;
           parameterMatrix = cat(1,parameterMatrix,RegMatrix);
        else
           RegMatrix = parametrizations.(parameterName).RegFun.fun(parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1);
           RegMatrix = sqrt(parametrizations.(parameterName).RegPar)*RegMatrix*inputParams/secDim;;
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
            parameterMatrix = cat(1,parameterMatrix,zeros(max(numberOfParameters),parametrizations.(parameterName).end-parametrizations.(parameterName).start + 1,secDim));    %Add correct number of zeros to the end
        end
    end
end

function dYdX = updateAccordingToGrouping(dYdX,parameterName, parametrizations)
% This function updates jacobians according to different parameter
% groupings
%
%   USAGE:
%       dYdX = updateAccordingToGrouping(dYdX,parameterName, parametrizations,inputParams,SignalPart,NoiseSD)
%
%   INPUTS:
%       dYdX             = jacobian
%       parameterName    = name of parameter
%       parametrizations = paramterization options
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
%% Update jacobian lines
if ~isempty(parametrizations.(parameterName).gr)                                                % Has grouping turned on
    dYdXtemp = zeros(size(dYdX,1),max(parametrizations.(parameterName).gr.idx),size(dYdX,3));
    for ll = 1 : length(parametrizations.(parameterName).gr.idx)
        dYdXtemp(:,parametrizations.(parameterName).gr.idx(ll),:) = dYdXtemp(:,parametrizations.(parameterName).gr.idx(ll),:) + dYdX(:,ll,:);
    end
    dYdX = dYdXtemp;
end

end



function dYdX = updateJacobianBlock(dYdX,parameterName, parametrizations,inputParams,SignalPart,NoiseSD)
% This function updates jacobians according to the parametrization
%
%   USAGE:
%       dYdX = updateJacobianBlock(dYdX,parameterName, parametrizations,inputParams,SignalPart,NoiseSD)
%
%   INPUTS:
%       dYdX             = jacobian
%       parameterName    = name of parameter
%       parametrizations = paramterization options
%       inputParams      = parameter values
%       SignalPart       = optimization signal part
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

    if ~isempty(NoiseSD)
        Sigma = NoiseSD.FD;                                                    % Get sigma
        Sigma = repmat(Sigma', [1 nPoints nLines]);                         % Repeat according to dimensions
        Sigma = permute(Sigma,[2 3 1]);                                     % Dims have to be nPoints nLines secDim
        Sigma = squeeze(Sigma);                                             % Remove zero dimensions
        dYdX = squeeze(dYdX);                                               % Remove zero dimensions
        dYdX(:,:,:) = dYdX(:,:,:) ./ Sigma;                                   % Normalize jacobian
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
        dYdX = squeeze(dYdX);                                               % Remove length 1 dims (e.g. for ph0, ph1, gaussLB)
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
        case 3   % E.g. Fixed metAmpls or single basis function case
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
            dYdX = reshape(dYdX,[],nLines);
        case  4 % E.g. free parametrization of per metabolite parameters
            dYdX = permute(dYdX,[1 3 4 2]);
            dYdX = reshape(dYdX,[],secDim*nLines);
    end
end

function [penalty, penaltyJac] = calcPenalty(inputParams, parametrizations, param, sD)
% Calculates the penalty term for deviations from expectation values
% relative to standard deviations
actualValue         = squeeze(inputParams.(param)(sD,:));
if ~isempty(parametrizations.(param).gr)
    idx = parametrizations.(param).gr.idx;
    actualValue =  actualValue(idx);
    actualValue =  actualValue(1:max(idx));
end
expectationValue    = parametrizations.(param).ex;
standardDeviation   = parametrizations.(param).sd;

% In the LCModel objective functions, the terms are in the shape of
% (difference-to-expectation-value)^2/(standard-deviation)^2. We are
% here formulating the vector of penalties that is appended to the
% residual, i.e., squared afterwards. We therefore give the penalty vector
% in the shape of (difference-to-expectation)/(standard-deviation).
diffVec     = actualValue-expectationValue;
sdVec       = standardDeviation;
penalty     = diffVec./sdVec;

% Jacobian
penaltyJac  = 1./sdVec;                % Because the penalty is a linear sum of its parts
end

function [penalty, penaltyJac] = calcSoftConstraintPenalty(inputParams, parametrizations, param, sD, nParams)
penalty = [];
penaltyJac = [];
if ~isempty(parametrizations.(param).sc)
    % Calculates the penalty term for deviations from soft constraints
    for mm = 1 : size(parametrizations.(param).sc.fix_idx,1)
        idx           = parametrizations.(param).sc.fix_idx(mm,:);
        idx(idx==0)   =[];
        tempfixValue  = squeeze(inputParams.(param)(sD,idx));
        fixValue(mm)  = sum(tempfixValue .* parametrizations.(param).sc.fix_factor{mm}');
    end
    for mm = 1 : size(parametrizations.(param).sc.adj_idx,1)
        idx           = parametrizations.(param).sc.adj_idx(mm,:);
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
    
    expectationValue = parametrizations.(param).sc.ex';
    standardDeviation   = parametrizations.(param).sc.sd';


    diffVec     = actualValue-expectationValue;
    sdVec       = standardDeviation;

    penalty  = zeros(1,length(squeeze(inputParams.(param)(sD,:))));


    penalty(parametrizations.(param).sc.adj_idx) = parametrizations.(param).sc.scaling*diffVec./sdVec;

    
    % Jacobian
    penaltyJac  = zeros(nParams,nParams);
    if sum(fixValue==0)==0 
        for m_adj = 1 : size(parametrizations.(param).sc.adj_idx,1)
            if fixValue(m_adj) ~= 0
                fix_idx           = parametrizations.(param).sc.fix_idx(m_adj,:);
                fix_idx(fix_idx==0)   =[];
                adj_idx           = parametrizations.(param).sc.adj_idx(m_adj,:);
                adj_idx(adj_idx==0)   =[];
                for m_fix = 1 : length(fix_idx)
                    penaltyJac(parametrizations.(param).start+adj_idx-1,parametrizations.(param).start+fix_idx(m_fix)-1) = -parametrizations.(param).sc.scaling*parametrizations.(param).sc.fix_factor{m_adj}(m_fix)*adjValue(m_adj)./(sdVec(m_adj).*fixValue(m_adj).*fixValue(m_adj));
                end
                penaltyJac(parametrizations.(param).start+adj_idx-1,parametrizations.(param).start+adj_idx-1) = parametrizations.(param).sc.scaling*1./(sdVec(m_adj).*fixValue(m_adj)); 
            end
        end
    end


else
    penalty  = zeros(1,length(squeeze(inputParams.(param)(sD,:))));
    if ~isempty(parametrizations.(param).gr)
        idx = parametrizations.(param).gr.idx;
        penalty =  penalty(idx);
        penalty =  penalty(1:max(idx));
    end 
    penaltyJac  = zeros(nParams,nParams);    
end
end
