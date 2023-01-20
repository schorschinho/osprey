function obj = createModel(obj)
            
    % This is the heart and soul of this class. This function takes
    % the parametrizations defined in the options, and translates
    % them into the appropriate loss and gradient functions.
    
    obj.step            = obj.step + 1;
    % Collect the basis functions
    basisSet            = obj.BasisSets;
    baselineBasis       = obj.BaselineBasis;
    data                = obj.Data.fids;
    ppm                 = obj.Data.ppm;
    t                   = obj.Data.t;
    Domain              = obj.Options{obj.step}.optimDomain;
    fitRange            = obj.Options{obj.step}.optimFreqFitRange;
    SignalPart          = obj.Options{obj.step}.optimSignalPart;
    solver              = obj.Options{obj.step}.solver;
    NormNoise           = obj.NormNoise;   
    
    
    % Only use basis functions that are included
    basisSet.fids   = basisSet.fids(:, logical(basisSet.includeInFit(obj.step,:)),:);
    
    % Create x0, lb, ub vectors by iteratively calling the
    % parameter-class-specific initialization
%     pars = {'ph0', 'ph1', 'gaussLB', 'metAmpl', 'freqShift', 'lorentzLB', 'baseAmpl'};
    pars = obj.Options{obj.step}.parameter;
    parsInit = [];
    parslb = [];
    parsub = [];
    parsfun = [];
    for pp = 1:length(pars)
        [parsInit, parslb, parsub, parsfun] = initializeParameters(obj, parsInit, parslb, parsub, parsfun, pars{pp});
    end

    % Set handles and create x0, lb, and ub vectors
    eval(['h = ' obj.Options{obj.step}.ModelFunction ';'])
    [x0, indexStruct] = h.pars2x(parsInit);
    [lb,~] = h.pars2x(parslb);
    [ub,~] = h.pars2x(parsub);
    
    % Update the parametrization according to the 2D json file. This is
    % crucial to define the type (free, fixed, dynamic) and new bounds
    for pp = 1:length(pars)
        obj.Options{obj.step}.parametrizations.(pars{pp}).start = indexStruct.(pars{pp}).start;
        obj.Options{obj.step}.parametrizations.(pars{pp}).end = indexStruct.(pars{pp}).end;
        obj.Options{obj.step}.parametrizations.(pars{pp}).init = parsInit.(pars{pp});
        obj.Options{obj.step}.parametrizations.(pars{pp}).lb = parslb.(pars{pp});
        obj.Options{obj.step}.parametrizations.(pars{pp}).ub = parsub.(pars{pp});
        obj.Options{obj.step}.parametrizations.(pars{pp}).type = parsfun.(pars{pp});
    end
    parametrizations = obj.Options{obj.step}.parametrizations;
    
    % Setup lossfunction outputs
    switch solver
        case {'lbfgsb', 'fminsearch'}
            sse = 'sos';
        case 'lsqnonlin'
            sse = 'res';
    end

    % Prepare the function wrapper
    fcn  = @(x) h.lossFunction(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange,SignalPart,Domain,sse,parametrizations);
   
    switch solver
        case 'lbfgsb'
            grad = @(x) h.forwardGradient(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange,SignalPart,parametrizations);
             % Request very high accuracy:
            opts            = struct('factr', 1e7, 'pgtol', 1e-2, 'm', 5);
            opts.printEvery = 1;
            % Run the algorithm:
            % Feed initial guess from the input parameters
            opts.x0 = x0;
            tstart = tic;
            [xk, ~, info] = lbfgsb({fcn,grad}, lb, ub, opts );
            time = toc(tstart);
        case 'lsqnonlin'
            jac = @(x) h.forwardJacobian(x, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange,SignalPart,parametrizations);
            fun  = @(x) h.fminunc_wrapper(x, fcn, jac);
            if obj.Options{obj.step}.CheckGradient
                CheckGrad = true;
            else
                CheckGrad = false;
            end
            opts = optimoptions('lsqnonlin','Display','iter','Algorithm','levenberg-marquardt','SpecifyObjectiveGradient',true,...
                                'CheckGradients',CheckGrad,'FiniteDifferenceType','central','MaxIterations',1000);
            tstart = tic;            
            [xk,info.resnorm,info.residual,info.exitflag,info.output,info.lambda,info.jacobian] = lsqnonlin(fun, x0, lb, ub, opts );
            time = toc(tstart);    
        case 'fminsearch'
            opts = optimset('fminsearch');
            opts.Display = 'iter';
            % Run the algorithm:
            % Feed initial guess from the input parameters
            tstart = tic;
            [xk, ~, info] = fminsearchbnd(fcn, x0, lb, ub, opts);
            time = toc(tstart);
    end
    

    
    % Save modeling results
    spec            = fftshift(fft(data, [], 1), 1);
    nBasisFcts      = size(basisSet.fids, 2);
    secDim          = size(data,2);      
    parsOut         = h.x2pars(xk, secDim, parametrizations);
    [fit, baseline, metabs] = h.forwardModel(xk, basisSet, baselineBasis, ppm, t, parametrizations);
    obj.Model{obj.step}.fit.fit       = fit;
    obj.Model{obj.step}.fit.baseline  = baseline;
    obj.Model{obj.step}.fit.residual  = spec - fit;
    obj.Model{obj.step}.fit.metabs    = metabs;
    obj.Model{obj.step}.time          = time;
    obj.Model{obj.step}.parsOut       = parsOut;
    obj.Model{obj.step}.info          = info;
    
    % Calculate CRLB
    [indMin, indMax] = h.ppmToIndex(ppm, fitRange);
    jac = h.forwardJacobian(xk, data, NormNoise, basisSet, baselineBasis, ppm, t, fitRange, 'C',parametrizations);
    SigmaSquared = (NormNoise * length(data(indMin:indMax)))^2;
    jac = -jac * SigmaSquared;
    
    % estimating the sigma based on the residual
    
    sigma   = std(real(obj.Model{obj.step}.fit.residual(indMin:indMax)));

    %calculate the fisher matrix
    fisher = (1./(sigma^2)) .* real(jac'*jac);
    
    %fisher = fisher + ones(size(fisher))*eps;
    invFisher = pinv(fisher);
    crlbs = sqrt(diag(invFisher));
    CRLB = h.x2pars(crlbs, secDim, parametrizations);
    
   
    %Relative CRLBs in percent
    relativeCRLB = CRLB.metAmpl./parsOut.metAmpl * 100;

    % Save table with basis function names and relative CRLB
    % Only use basis functions that are included
    obj.Model{obj.step}.rawCRLB = CRLB;
    try
        obj.Model{obj.step}.CRLB = table(basisSet.names(:, logical(basisSet.includeInFit(obj.step,:)))', relativeCRLB);
    catch
    end

end