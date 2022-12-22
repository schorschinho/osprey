function updateOptsAccordingToStep(obj, options)
    % Initialize an empty container
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

    % Save the property struct
    obj.Options{obj.step+1} = options;
end