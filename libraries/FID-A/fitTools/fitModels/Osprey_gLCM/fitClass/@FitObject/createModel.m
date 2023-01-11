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
    
    
    % Only use basis functions that are included
    basisSet.fids   = basisSet.fids(:, logical(basisSet.includeInFit(obj.step,:)),:);
    
    % Create x0, lb, ub vectors by iteratively calling the
    % parameter-class-specific initialization
%     pars = {'ph0', 'ph1', 'gaussLB', 'metAmpl', 'freqShift', 'lorentzLB', 'baseAmpl'};
    pars = obj.Options{obj.step}.parameter;
    parsInit = [];
    parslb = [];
    parsub = [];
    for pp = 1:length(pars)
        [parsInit, parslb, parsub] = initializeParameters(obj, parsInit, parslb, parsub, pars{pp});
    end

    eval(['h = ' obj.Options{obj.step}.ModelFunction ';'])
    x0 = h.pars2x(parsInit);
    lb = h.pars2x(parslb);
    ub = h.pars2x(parsub);
    
    % Setup lossfunction outputs
    switch solver
        case {'lbfgsb', 'fminsearch'}
            sse = 'sos';
        case 'lsqnonlin'
            sse = 'res';
    end

    % Prepare the function wrapper
    fcn  = @(x) h.lossFunction(x, data, basisSet, baselineBasis, ppm, t, fitRange,SignalPart,Domain,sse);
   
    switch solver
        case 'lbfgsb'
            grad = @(x) h.forwardGradient(x, data, basisSet, baselineBasis, ppm, t, fitRange,SignalPart);
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
            jac = @(x) h.forwardJacobian(x, data, basisSet, baselineBasis, ppm, t, fitRange,SignalPart);
            fun  = @(x) h.fminunc_wrapper(x, fcn, jac);
            if obj.Options{obj.step}.CheckGradient
                CheckGrad = true;
            else
                CheckGrad = false;
            end
            opts = optimoptions('lsqnonlin','Display','iter','Algorithm','levenberg-marquardt','SpecifyObjectiveGradient',true,'CheckGradients',CheckGrad,'FiniteDifferenceType','central');
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
    parsOut         = h.x2pars(xk, nBasisFcts);
    [fit, baseline, metabs] = h.forwardModel(xk, basisSet, baselineBasis, ppm, t);
    obj.Model{obj.step}.fit.fit       = fit;
    obj.Model{obj.step}.fit.baseline  = baseline;
    obj.Model{obj.step}.fit.residual  = spec - fit;
    obj.Model{obj.step}.fit.metabs    = metabs;
    obj.Model{obj.step}.time          = time;
    obj.Model{obj.step}.parsOut       = parsOut;
    obj.Model{obj.step}.info          = info;
    
    % Calculate CRLB
    jac = h.forwardJacobian(xk, data, basisSet, baselineBasis, ppm, t, fitRange, 'C');
    jac = -jac;
    
    % estimating the sigma based on the residual
    [indMin, indMax] = h.ppmToIndex(ppm, fitRange);
    sigma   = std(real(obj.Model{obj.step}.fit.residual(indMin:indMax)));

    %calculate the fisher matrix
    fisher = (1./(sigma^2)) .* real(jac'*jac);
    
    %fisher = fisher + ones(size(fisher))*eps;
    invFisher = pinv(fisher);
    crlbs = sqrt(diag(invFisher));
    CRLB = h.x2pars(crlbs, nBasisFcts);
    
   
    %Relative CRLBs in percent
    relativeCRLB = CRLB.metAmpl./parsOut.metAmpl * 100;

    % Save table with basis function names and relative CRLB
    % Only use basis functions that are included
    obj.Model{obj.step}.rawCRLB = CRLB;
    obj.Model{obj.step}.CRLB = table(basisSet.names(:, logical(basisSet.includeInFit(obj.step,:)))', relativeCRLB);

end