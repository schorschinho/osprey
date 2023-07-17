function updateOptsAccordingToStep(obj, options)
%%  updateOptsAccordingToStep(obj, options)
%   This is the method updates the parametrization options in the
%   OspreyFitObj according to the model procedure step
%
%   USAGE:
%       obj.updateOptsAccordingToStep(options)

%   INPUTS:
%       options            = struct with parametrization options   
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
%% Update parameters
    if ~isfield(options, 'optimDomain')
        options.optimDomain = 'FD'; % FD, TD, or FDTD
    else
        if ~(strcmpi(options.optimDomain,'FD') ||...
                strcmpi(options.optimDomain,'TD') ||...
                strcmpi(options.optimDomain,'FDTD'))
            error('Invalid optimization domain specification (options.optimDomain): (%s).', options.optimDomain)
        end
    end
    
    if ismember(options.optimDomain, {'FD', 'FDTD'})
        if ~isfield(options, 'optimFreqFitRange')
            options.optimFreqFitRange = [0.5 4.0];
        end
    end
    
    if ismember(options.optimDomain, {'TD', 'FDTD'})
        if ~isfield(options, 'optimTimeFitRange')
            options.optimTimeFitRange = [0 1];
        end
    end
    
    if ~isfield(options, 'optimSignalPart')
        options.optimSignalPart = 'R'; % R, I, or RI
    end
    
    
    if ~isfield(options, 'parametrizations')
        % Initialize phi0 as constant with value 0
        options.parametrizations.ph0.fun     = 'free';
        options.parametrizations.ph0.gradfun = 'free';
        options.parametrizations.ph0.lb      = -pi;
        options.parametrizations.ph0.ub      = pi;
        options.parametrizations.ph0.init    = 0;
        
        % Initialize phi1 as constant with value 0
        options.parametrizations.ph1.fun     = 'free';
        options.parametrizations.ph1.gradfun = 'free';
        options.parametrizations.ph1.lb      = -pi/30;
        options.parametrizations.ph1.ub      = pi/30;
        options.parametrizations.ph1.init    = 0;
        
        % Initialize Gaussian LB as constant with value [0.04 *
        % hz/ppm]
        options.parametrizations.gaussLB.fun     = 'free';
        options.parametrizations.gaussLB.gradfun = 'free';
        options.parametrizations.gaussLB.lb      = 0;
        options.parametrizations.gaussLB.ub      = sqrt(5000);
        options.parametrizations.gaussLB.init    = 0.04 * obj.Data.txfrq*1e-6;
        
        % Initialize Lorentzian LB as constant with value 2 Hz
        options.parametrizations.lorentzLB.fun     = 'free';
        options.parametrizations.lorentzLB.gradfun = 'free';
        options.parametrizations.lorentzLB.lb      = 0;
        options.parametrizations.lorentzLB.ub      = 10;
        options.parametrizations.lorentzLB.init    = 2;
        
        % Initialize frequency shifts as constant with value 0 Hz
        options.parametrizations.freqShift.fun     = 'free';
        options.parametrizations.freqShift.gradfun = 'free';
        options.parametrizations.freqShift.lb      = -5;
        options.parametrizations.freqShift.ub      = 5;
        options.parametrizations.freqShift.init    = 0;
        
        % Initialize metabolite amplitudes as free with value 0
        options.parametrizations.metAmpl.fun     = 'free';
        options.parametrizations.metAmpl.gradfun = 'free';
        options.parametrizations.metAmpl.lb      = 0;
        options.parametrizations.metAmpl.ub      = Inf;
        options.parametrizations.metAmpl.init    = 0;
        
        % Initialize baseline amplitudes as free with value 0
        options.parametrizations.baseAmpl.fun     = 'free';
        options.parametrizations.baseAmpl.gradfun = 'free';
        options.parametrizations.baseAmpl.lb      = -Inf;
        options.parametrizations.baseAmpl.ub      = Inf;
        options.parametrizations.baseAmpl.init    = 0;
        
        % Initialize x (the external dependency vector) as natural numbers
        options.parametrizations.x.values = [1:1:size(obj.Data.fids,2)];
        options.parametrizations.x.name   = 'independentVariable';
        
    end
    
    % Update baseline model 
    switch options.baseline.type                      % Switch for baseline type
        case 'spline'                                                   
            %%% CREATE BASELINE SPLINE BASIS %%%
            % Combine real and imaginary part to form a complex spline array.
            % Use the new, corrected function from here on                  
            dkntmn      = options.baseline.dkntmn;            % Get spline basis knot spacing
            [splineArray] = osp_gLCM_makeSplineBasis(obj.Data, options.optimFreqFitRange, dkntmn);   % Create spline baseline basis array     
            obj.BaselineBasis = splineArray;                                    % Store baseline array in object
        case 'poly'
            %%% CREATE BASELINE POLYNOMIAL BASIS %%%
            order      = options.baseline.order;              % Get order of the polynomial baseline
            [splineArray] = osp_gLCM_makePolyBasis(obj.Data, options.optimFreqFitRange, order);      % Create polynomial baseline basis array       
            obj.BaselineBasis = splineArray;                                    % Store baseline array in object
        case 'none'
            %%% NO BASELINE %%%
            obj.BaselineBasis = [];                                             % Store empty baseline array in object
    end

    obj.Options{obj.step+1} = options;                                      % Store options in OspreyFitObj
end