function obj = createModel(obj)
%%  obj = createModel(obj)
%   This is the heart and soul of this class. This method takes
%   the parametrizations defined in the options, and translates
%   them into the appropriate loss and gradient functions.   
%
%   USAGE:
%       obj.createModel()
%
%   OUTPUTS:
%       obj     = OspreyFitObj.
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
%%  Get object details
                
    obj.step            = obj.step + 1;                      % Update fit step counter
   
    if isfield(obj.Options{obj.step}, 'InitialPick')         % If an initial fit should be performed on a single spectrum we need to restore pick this single spectrum
        if obj.Options{obj.step}.InitialPick > 0             % Initial spectrum defined for 2D modeling
            temp.BasisSets.fids  = obj.BasisSets.fids;       % Backup full basis set
            temp.Data.fids       = obj.Data.fids;            % Backup full data
            obj.BasisSets.fids   = squeeze(obj.BasisSets.fids(:,:,obj.Options{obj.step}.InitialPick)); % Extract basis set for initial fit
            obj.Data.fids        = obj.Data.fids(:,obj.Options{obj.step}.InitialPick);  % Extract data for initial fit          
        end
    end
    % Basis set treatment
    basisSet            = obj.BasisSets;                                    % Get basis set
    % Subset the basis set to only include the basis functions included
    basisSet.fids       = basisSet.fids(:, logical(basisSet.includeInFit(obj.step,:)),:); % Get FIDs
    basisSet.names      = basisSet.names(:, logical(basisSet.includeInFit(obj.step,:)));  % Get names
    baselineBasis       = obj.BaselineBasis;                                % Get baseline basis
    data                = obj.Data.fids;                                    % Get time domain data
    ppm                 = obj.Data.ppm;                                     % Get ppm axis
    t                   = obj.Data.t;                                       % Get time vector
    Domain              = obj.Options{obj.step}.optimDomain;                % Get optimization domain
    if ~isempty(obj.Options{obj.step}.optimFreqFitRange)   
        fitRange.FD     = obj.Options{obj.step}.optimFreqFitRange;          % Get model range
    end
    if isfield(obj.Options{obj.step}, 'optimTimeFitRange') && ~isempty(obj.Options{obj.step}.optimTimeFitRange) 
        fitRange.TD     = obj.Options{obj.step}.optimTimeFitRange;          % Get model range
    end
    fitGap              = obj.Options{obj.step}.gap;                        % Get the gap in the model range
    SignalPart          = obj.Options{obj.step}.optimSignalPart;            % Get optimization signal part
    solver              = obj.Options{obj.step}.solver;                     % Get solver name
    NoiseSD             = obj.NoiseSD;                                      % Get standard deviation of the noise 
    % basisSet.fids       = basisSet.fids(:, logical(basisSet.includeInFit(obj.step,:)),:); % Only use basis functions that are included

%%  Update parameterizations according to model procedure step

    

    % Create x0, lb, ub vectors by iteratively calling the
    % parameter-class-specific initialization
    pars = obj.Options{obj.step}.parameter;                 % Get parameter names
    parsInit = [];                                          % Initialize parsInit struct
    parslb = [];                                            % Initialize parslb struct
    parsub = [];                                            % Initialize parsub struct
    parsex = [];                                            % Initialize parsex struct
    parssd = [];                                            % Initialize parssd struct
    parsfun = [];                                           % Initialize parsfun struct
    parsgr = [];                                            % Initialize parsgr struct
    parssc = [];                                            % Initialize parssc struct
    for pp = 1:length(pars)                                 % Loop over parameters
        [parsInit, parslb, parsub, parsex, parssd, parsfun, parsgr, parssc] = initializeParameters(obj, parsInit, parslb, parsub, parsex, parssd, parsfun, parsgr, parssc, pars{pp}); % Generate parameter structs
    end                                                     % End loop over parameters

    % Get indices for parameter soft constraints
    pars = fields(parssc);                                        % Get parameter names
    for ff = 1 : length(pars)                                     % Loop over parameters
        if ~isempty(parssc.(pars{ff}))                            % Grouping exists
            if isfield(parssc.(pars{ff}),'fix_idx')                         % Remove old index during optimization
               parssc.(pars{ff}) = rmfield(parssc.(pars{ff}),'fix_idx');
               parssc.(pars{ff}) = rmfield(parssc.(pars{ff}),'adj_idx'); 
            end
            basisNames = basisSet.names;
            
            for rr = 1 : length(parssc.(pars{ff}).fix)
                [metsToIncludeFixTemp, ~, ~] = intersect(basisNames, parssc.(pars{ff}).fix{rr}, 'stable');
                [metsToIncludeAdjTemp, ~, ~] = intersect(basisNames, parssc.(pars{ff}).adj{rr}, 'stable');
                
                for mm = 1 : length(metsToIncludeFixTemp)
                    idxToIncludeFix(rr,mm) = find(strcmp(metsToIncludeFixTemp{mm}, basisNames));    % Get index of basis function to include
                end
                for mm = 1 : length(metsToIncludeAdjTemp)
                    idxToIncludeAdj(rr,mm) = find(strcmp(metsToIncludeAdjTemp{mm}, basisNames));    % Get index of basis function to include
                end
            end
            parssc.(pars{ff}).fix_idx = idxToIncludeFix;
            parssc.(pars{ff}).adj_idx = idxToIncludeAdj;
            if ~iscell(parssc.(pars{ff}).fix_factor)
                for mm = 1:length(parssc.(pars{ff}).fix_factor)
                    temp(mm) = {parssc.(pars{ff}).fix_factor(mm)};                    
                end
                parssc.(pars{ff}).fix_factor = temp;
            end
        end
    end

    % Get indices for parameter mapping according to groups
    pars = fields(parsgr);                                        % Get parameter names
    for ff = 1 : length(pars)                                     % Loop over parameters
        if ~isempty(parsgr.(pars{ff}))                            % Grouping exists
            if isfield(parsgr.(pars{ff}),'idx')                         % Remove old index during optimization
               parsgr.(pars{ff}) = rmfield(parsgr.(pars{ff}),'idx'); 
            end
            groups = fields(parsgr.(pars{ff}));                   % Get group names
            
            basisNames = basisSet.names;
            idx = 1:sum(basisSet.includeInFit(obj.step,:));
            for gg = 1 : length(groups)            
                [metsToInclude, ~, ~] = intersect(basisNames, parsgr.(pars{ff}).(groups{gg}), 'stable');    % Get vector of logical indices
                firstIndex = find(strcmp(metsToInclude{1}, basisNames));    % Get index of basis function to include
                for rr = 2:length(metsToInclude)
                    idxToInclude = find(strcmp(metsToInclude{rr}, basisNames));    % Get index of basis function to include
                    idx(idxToInclude) = firstIndex;
                end 
                parsgr.(pars{ff}).idx = idx;
            end
            % Update paramters accoring to grouping
            parsInit.(pars{ff}) =  parsInit.(pars{ff})(idx);
            parslb.(pars{ff}) =  parslb.(pars{ff})(idx);
            parsub.(pars{ff}) =  parsub.(pars{ff})(idx);
            parsex.(pars{ff}) =  parsex.(pars{ff})(idx);
            parssd.(pars{ff}) =  parssd.(pars{ff})(idx);
            % Remove parameters that are not needed anymore
            parsInit.(pars{ff}) =  parsInit.(pars{ff})(1:max(idx));
            parslb.(pars{ff}) =  parslb.(pars{ff})(1:max(idx));
            parsub.(pars{ff}) =  parsub.(pars{ff})(1:max(idx));
            parsex.(pars{ff}) =  parsex.(pars{ff})(1:max(idx));
            parssd.(pars{ff}) =  parssd.(pars{ff})(1:max(idx));

        end
    end

    eval(['h = ' obj.Options{obj.step}.ModelFunction ';'])  % Set model function handle from ModelFunction field (e.g. GeneralizedPhysicsModel)
    [x0, indexStruct] = h.pars2x(parsInit);                 % Create x0 vector
    [lb,~] = h.pars2x(parslb);                              % Create lb vector
    [ub,~] = h.pars2x(parsub);                              % Create ub vector
    [ex,~] = h.pars2x(parsex);                              % Create ex vector
    [sd,~] = h.pars2x(parssd);                              % Create sd vector

    % Update the parametrization according to the 2D json file. This is
    % crucial to define the type (free, fixed, dynamic) and new bounds
    for pp = 1:length(pars)                                 % Loop over parameters
        if isfield(indexStruct,pars{pp})                    % Avoid any entries in the struct that are not parameters
            obj.Options{obj.step}.parametrizations.(pars{pp}).start = indexStruct.(pars{pp}).start; % Update start index for x vector
            obj.Options{obj.step}.parametrizations.(pars{pp}).end = indexStruct.(pars{pp}).end; % Update end index for x vector
            obj.Options{obj.step}.parametrizations.(pars{pp}).init = parsInit.(pars{pp}); % Update init values
            obj.Options{obj.step}.parametrizations.(pars{pp}).lb = parslb.(pars{pp});  % Update lb values
            obj.Options{obj.step}.parametrizations.(pars{pp}).ub = parsub.(pars{pp}); % Update ub values
            obj.Options{obj.step}.parametrizations.(pars{pp}).ex = parsex.(pars{pp}); % Update ex values
            obj.Options{obj.step}.parametrizations.(pars{pp}).sd = parssd.(pars{pp}); % Update sd values
            obj.Options{obj.step}.parametrizations.(pars{pp}).type = parsfun.(pars{pp}); % Update type strings
            obj.Options{obj.step}.parametrizations.(pars{pp}).gr = parsgr.(pars{pp}); % Update grouping information
            obj.Options{obj.step}.parametrizations.(pars{pp}).sc = parssc.(pars{pp}); % Update grouping information
        end
    end                                                     % End loop over parameters
    
    if isempty(obj.BaselineBasis)                                           % baseline setting none
        obj.Options{obj.step}.parametrizations.baseAmpl.start = 0;          % No baseline start index for x vector
        obj.Options{obj.step}.parametrizations.baseAmpl.end = 0;            % No baseline end index for x vector
        obj.Options{obj.step}.parametrizations.baseAmpl.type = 'none';      % Type none
    end
    
    parametrizations = obj.Options{obj.step}.parametrizations;              % Write parametrization in variable for solver
    
    if sum(cellfun(@isstruct,obj.returnParametrization(obj.step,'RegFun'))) > 0    % Apply regularizer?
        Reg = 1;                                                            % Set regularizer to yes
    else
        Reg = 0;                                                            % Set regularizer to no
    end

%% Setup inputs for optimizers 
   
    switch solver                                        % Switch to setup loss function outputs according to solver
        case {'lbfgsb', 'fminsearch'}
            sse = 'sos';                                 % Uses sum of squares
        case 'lsqnonlin'
            sse = 'res';                                 % Uses residual vector
    end

    % Set lossfunction handle
    fcn  = @(x) h.lossFunction(x, ...               % x vector with parameters to optimize
                               data, ...            % time domain data matrix
                               NoiseSD, ...         % standard deviation of the noise 
                               basisSet, ...        % basis set struct
                               baselineBasis, ...   % baseline basis set
                               ppm, ...             % ppm axis
                               t, ...               % time vector
                               fitRange, ...        % model range
                               fitGap, ...          % model gap
                               SignalPart, ...      % optimization signal part
                               Domain, ...          % optimization domain
                               sse, ...             % lossfunction string
                               Reg, ...             % regularizer flag
                               parametrizations);   % parameter struct
    
    switch solver                                        % Switch to pick solver
        case 'lbfgsb'
            % Set gradient function handle
            grad = @(x) h.forwardGradient(x, ...            % x vector with parameters to optimize
                                          data, ...         % time domain data matrix
                                          NoiseSD, ...      % standard deviation of the noise 
                                          basisSet, ...     % basis set struct
                                          baselineBasis, ...% baseline basis set
                                          ppm, ...          % ppm axis
                                          t, ...            % time vector
                                          fitRange, ...     % model range
                                          SignalPart, ...   % optimization signal part
                                          Reg, ...          % regularizer flag
                                          parametrizations);% parameter struct
            
            opts            = struct('factr', 1e7, 'pgtol', 1e-5, 'm', 5, 'printEvery', 0); % Set solver options
            opts.x0 = x0';

            tstart = tic;                                               % Start timer
            [xk, ~, info] = lbfgsb({fcn,grad}, lb', ub', opts );        % Run solver
            time = toc(tstart);                                         % End timer

        case 'lsqnonlin'               % Levenberg-Marquardt solver
            if obj.Options{obj.step}.NumericJacobian                    % Use numerically calculated jacobian (slow)
                SpecifyObjectiveGradient = false;
            else
                SpecifyObjectiveGradient = true;
            end
            if obj.Options{obj.step}.CheckGradient                      % Perform gradient check (for debugging)
                CheckGrad = true;
            else
                CheckGrad = false;
            end
            % Set jacobian handle
            jac = @(x) h.forwardJacobian(x, ...             % x vector with parameters to optimize
                                         data, ...          % time domain data matrix
                                         NoiseSD, ...       % standard deviation of the noise 
                                         basisSet, ...      % basis set struct
                                         baselineBasis, ... % baseline basis set
                                         ppm, ...           % ppm axis
                                         t, ...             % time vector
                                         fitRange, ...      % model range
                                         fitGap, ...          % model gap
                                         SignalPart, ...    % optimization signal part
                                         Domain,...         % optimization domain
                                         Reg, ...           % regularizer flag
                                         parametrizations); % parameter struct
            % Set fminunc wrappper handle
            fun  = @(x) h.fminunc_wrapper(x, fcn, jac);
             % Set solver options
            opts = optimoptions('lsqnonlin', ...
                                'Algorithm','levenberg-marquardt', ...      % Use LM
                                'SpecifyObjectiveGradient',SpecifyObjectiveGradient,... % Use analytic jacobian
                                'CheckGradients',CheckGrad, ...             % Check gradient
                                'FiniteDifferenceType','central', ...       % for numerically calculated jacobian only
                                'MaxIterations',1000, ...                   % Iterations
                                'Display','none');                       % Display no iterations

            % Add this if you want to plot per iteration
            % 'OutputFcn',@optimplotresidual,...);                          
           


            % Set tolerance this is useful for the regularization
            % optimization
            if isfield(obj.Options{obj.step},'FunctionTolerance')
                opts.FunctionTolerance = obj.Options{obj.step}.FunctionTolerance;
            end
            if isfield(obj.Options{obj.step},'StepTolerance')
                opts.StepTolerance = obj.Options{obj.step}.StepTolerance;
            end
            if isfield(obj.Options{obj.step},'OptimalityTolerance')
                opts.OptimalityTolerance = obj.Options{obj.step}.OptimalityTolerance;
            end

            tstart = tic;                                               % Start timer           
            [xk,info.resnorm,info.residual,info.exitflag,info.output,info.lambda,info.jacobian] = lsqnonlin(fun, x0, lb, ub, opts ); % Run solver
            time = toc(tstart);                                         % End timer 

        case 'fminsearch'               % Simplex algorithm uses no jacobian
            opts = optimset('fminsearch');                  % Set solver options

            tstart = tic;                                               % Start timer
            [xk, ~, info] = fminsearchbnd(fcn, x0, lb, ub, opts);       % Run solver
            time = toc(tstart);                                         % End timer

    end

%% Store outputs

    spec            = fftshift(fft(data, [], 1), 1);                % Get spectra
    secDim          = size(data,2);                                 % Number of entries along indirect dimension     
    parsOut         = h.x2pars(xk, secDim, parametrizations);       % Get final parametes

    % Generate final model for plots
    [fit, baseline, metabs] = h.forwardModel(xk, ...                        % x vector with final parameter estimates
                                             basisSet, ...                  % basis set struct
                                             baselineBasis, ...             % baseline basis set
                                             ppm, ...                       % ppm axis
                                             t, ...                         % time vector
                                             Reg, ...                       % regularizer flag
                                             parametrizations);             % parameter struct

    
    % Save modeling results
    obj.Model{obj.step}.fit.fit       = fit;                                % Store fit
    obj.Model{obj.step}.fit.baseline  = baseline;                           % Store baseline
    obj.Model{obj.step}.fit.residual  = spec - fit;                         % Store residual
    obj.Model{obj.step}.fit.metabs    = metabs;                             % Store metabolite results
    obj.Model{obj.step}.time          = time;                               % Store time vector
    obj.Model{obj.step}.parsOut       = parsOut;                            % Store final parameters
    obj.Model{obj.step}.info          = info;                               % Store info

    % Baseline model is 'residual' we have still generate the baseline
    if strcmp(obj.Options{obj.step}.baseline.type, 'residual')
        obj.residualSmoothing;
    end
%% Calculate CRLB

    jac = h.forwardJacobian(xk, ...            % x vector with final parameter estimates
                            data, ...          % time domain data matrix
                            NoiseSD, ...       % standard deviation of the noise 
                            basisSet, ...      % basis set struct
                            baselineBasis, ... % baseline basis set
                            ppm, ...           % ppm axis
                            t, ...             % time vector
                            fitRange, ...      % model range
                            fitGap, ...          % model gap
                            'R', ...           % get complex-valued jacobian
                            Domain,...         % optimization domain
                            Reg, ...           % regularizer flag
                            parametrizations); % parameter struct

    fisher = real(jac'*jac);                            % calculate the fisher matrix
    
    invFisher = pinv(fisher);                           % invert fisher matrix
    crlbs = sqrt(diag(invFisher));                      % get raw CRLBs values
    CRLB = h.x2pars(crlbs, secDim, parametrizations);   % convert CRLBs to parameter struct
    
    relativeCRLB = CRLB.metAmpl(1,:)./parsOut.metAmpl(1,:) * 100; % Relative CRLBs for amplitudes

    
    obj.Model{obj.step}.rawCRLB = CRLB;                 % Save raw CRLBs
    try
        obj.Model{obj.step}.CRLB = array2table(relativeCRLB,'VariableNames',basisSet.names); % Save table with basis function names and relative CRLBs
    catch
    end    
    
    % Calculate CRLBs for metabolite combinations
    metaboliteNames = basisSet.names;
    switch solver                                        % Switch to pick solver
        case 'lbfgsb'
            [obj] = calculateCombinedCRLB(obj,invFisher, xk', metaboliteNames, parametrizations);
        otherwise
            [obj] = calculateCombinedCRLB(obj,invFisher, xk, metaboliteNames, parametrizations);
    end
%% Do some clean-up    
   
    % If an inital fit was performed on a single spectrum we need to
    % restore the original dimensions
    if isfield(obj.Options{obj.step}, 'InitialPick')
        if obj.Options{obj.step}.InitialPick > 0
            obj.BasisSets.fids  = temp.BasisSets.fids;
            obj.Data.fids       = temp.Data.fids;
        end
    end


%% Calculate AICs & Fit quality number

    % Generate sum-of-squares for final model
    % Baseline model is 'residual' have a different residual vector
    if ~strcmp(obj.Options{obj.step}.baseline.type, 'residual')
        sse = 'res';                                 % Uses residual vector
        res  =      h.lossFunction(xk, ...               % x vector with final parameters estimates
                                   data, ...            % time domain data matrix
                                   NoiseSD, ...         % standard deviation of the noise 
                                   basisSet, ...        % basis set struct
                                   baselineBasis, ...   % baseline basis set
                                   ppm, ...             % ppm axis
                                   t, ...               % time vector
                                   fitRange, ...        % model range
                                   fitGap, ...          % model gap
                                   SignalPart, ...      % optimization signal part
                                   Domain, ...          % optimization domain
                                   sse, ...             % lossfunction string
                                   Reg, ...             % regularizer flag
                                   parametrizations);   % parameter struct
    else
        res = obj.Model{obj.step}.fit.residual;
    end
    sos = sum(res.^2); 
    % Calculate different AICs
    p       = length(xk);                                                       % Number of estimated parameter
    n       = length(res);                                                      % Number of points in model
    sigma   = sqrt(sos/n);                                                      % Sigma calculation for AIC
    obj.Model{obj.step}.AIC     = - 2 * n * log(sigma) - 2*p;                   % AIC
    obj.Model{obj.step}.AIC_c   = - 2 * n * log(sigma) - 2*(p+1)*(n/(n-p-2));   % AIC_c
    obj.Model{obj.step}.BIC     = - 2 * n * log(sigma) - log(n)*p;              % BIC
    obj.Model{obj.step}.fitQAnumber =  sos/length(res);
end