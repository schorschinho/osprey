function obj = optimizeRegularization(obj, opts)
%%  optimizeRegularization(obj, opts)
%   This is the method that runs the optimization loop top find the optimal
%   regularizer parameter lambda
%
%   USAGE:
%       obj.optimizeRegularization(opts)
%
%   INPUTS:
%       opts   = options struct for the regularizer  
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
    
    ppm                     = obj.Data.ppm;                                 % Get ppm range
    optimFreqFitRange       = obj.Options{1, end}.optimFreqFitRange;        % Get fit range
    SplineBasis             = obj.BaselineBasis;                            % Get spline basis 
    lambdaRange.ub          = opts.regularizer.ub;                          % Get upper bound of regularization parameter
    if isfield(opts.regularizer,'lb')
       lambdaRange.lb       = opts.regularizer.lb;                          % Get lower bound of regularization parameter
    else
        lambdaRange.lb      = [];                                           % Set to default value [] which will results in 2.01 ED per ppm
    end
    if isfield(opts.regularizer,'m')
       m       = opts.regularizer.m;                                        % Get m value for modified AIC
    else
       m      = 1;                                                          % Use normal AIC definition
    end
    OptimSteps                   = opts.regularizer.steps;                  % Get number of regularization steps    
  
%%  Set handles for the regularizer formalism

    eval(['h_Reg = ' obj.Options{1, end}.parametrizations.baseAmpl.RegFun ';']) % Get handle to regualrizer matrix for example second order differencing matrix SecondDerivativ
    eval(['h_RegOpt = ' opts.regularizer.fun ';'])                              % Get handle to the optimization function of the regualrizer for example LambdaEDFormalism

%% Run regularization parameter optimization

    % Set regularizer space
    if OptimSteps > 1   % Set up the regularizer parameter space according to effective dimension to lambda formalism 
        RegParSpace = h_RegOpt.SetLambdaSpace(ppm,optimFreqFitRange,SplineBasis,lambdaRange,h_Reg,OptimSteps); %Use SetLambdaSpace to get regularization parameter space
    else % If only one opimization step is performed we don't need a parameter space.
        RegParSpace = h_RegOpt.EDtoL(real(SplineBasis),real(h_Reg.fun(size(SplineBasis,2))),lambdaRange.ub*abs(optimFreqFitRange(1)-optimFreqFitRange(2))); %Use EDtoL to get regularization parameter 
    end
    residual = zeros(size(obj.Data.fids,1),OptimSteps);                     % Set residual matrix

    % Loop over optimizer steps
    for oo = 1 : OptimSteps                                                 % Loop over optimizer steps
        obj.Options{1, end}.parametrizations.(opts.regularizer.parameter).RegPar = RegParSpace(oo); % Update regularization parameter of OspreyFitObj according to regularizaiton paramter space    
        if oo > 1  || obj.step > 1                                          % If we are in the second step or have had an initial step get inital estiamtes for non amplitude paramters and fix them
            if isfield(opts.regularizer,'fixedParameter')                   % Get parameter names that should be fixed during the regularizer optimization          
                for pars = 1 : length(opts.regularizer.fixedParameter)      % Loop over parameters to fix
                    obj.Options{1, end}.parametrizations.(opts.regularizer.fixedParameter{pars}).init = obj.Model{1, 1}.parsOut.(opts.regularizer.fixedParameter{pars}); % Get initials from first step
                    obj.Options{1, end}.parametrizations.(opts.regularizer.fixedParameter{pars}).ub = obj.Model{1, 1}.parsOut.(opts.regularizer.fixedParameter{pars}); % Fix parameter
                    obj.Options{1, end}.parametrizations.(opts.regularizer.fixedParameter{pars}).lb = obj.Model{1, 1}.parsOut.(opts.regularizer.fixedParameter{pars}); % Fix parameter
                end                                                         % End loop over parameters to fix
            end
            obj.step            = obj.step - 1;                             % We have to subtract a step again because we are working on the same instance of the OspreyFitObj
        end
        obj.createModel();                                                  % Run model
        residual(:,oo)=obj.Model{1, obj.step}.fit.residual;
    end                                                                     % End loop over optimizer steps
    
    % Find the optimal regularization parameter using the mAIC formalism
    if OptimSteps>1                                                         % This is only needed when optimzation is actually ran
        [optimLambda, optimLambdaIndex,ed,AIC] = h_RegOpt.optimLambda(residual,ppm,optimFreqFitRange,SplineBasis,h_Reg,OptimSteps,m,RegParSpace);   % Use optimLambda to find optimal regularization parameter
        % Store optimal regularizer parameters in OspreyFitObj for next
        % step. Also store the full ed and AIC space to allow plotting
        % those later.
        obj.Model{obj.step}.Regularization.OptimalRegParIndex = optimLambdaIndex;   % Store index to optimal regularization parameter in OspreyFitObj
        obj.Model{obj.step}.Regularization.OptimalRegPar = optimLambda;             % Store optimal regularization parameter in OspreyFitObj
        obj.Model{obj.step}.Regularization.ed = ed;                                 % Store ED array in OspreyFitObj
        obj.Model{obj.step}.Regularization.AIC = AIC;                               % Store AIC array in OspreyFitObj
        obj.Options{1, end}.parametrizations.(opts.regularizer.parameter).RegPar = optimLambda; % Get optimal regularization parameter to re-run object for visualization
        if isfield(opts.regularizer,'fixedParameter')                       % Get parameter names that should be fixed during the regularizer optimization             
            for pars = 1 : length(opts.regularizer.fixedParameter)          % Loop over parameters to fix
                obj.Options{1, end}.parametrizations.(opts.regularizer.fixedParameter{pars}).init = obj.Model{1, 1}.parsOut.(opts.regularizer.fixedParameter{pars}); % Get initials from first step
                obj.Options{1, end}.parametrizations.(opts.regularizer.fixedParameter{pars}).ub = obj.Model{1, 1}.parsOut.(opts.regularizer.fixedParameter{pars}); % Fix parameter
                obj.Options{1, end}.parametrizations.(opts.regularizer.fixedParameter{pars}).lb = obj.Model{1, 1}.parsOut.(opts.regularizer.fixedParameter{pars}); % Fix parameter
            end                                                             % End loop over parameters to fix
        end
        obj.step            = obj.step - 1;                                 % We have to subtract a step again because we are working on the same instance of the OspreyFitObj
        obj.createModel();                                                  % Run model
    end

end