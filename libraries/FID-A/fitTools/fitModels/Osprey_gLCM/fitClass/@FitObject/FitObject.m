classdef FitObject < handle
%%  FitObject
%   This is the class defintion for the OspreyFitObj. This is the
%   center-piece of the new fit formalism in Osprey
%
%   USAGE:
%       obj = FitObject(DataToModel, basisSet, opts);
%
%   INPUTS:
%       DataToModel     = FID-A MRS struct (2D is modeled along extra dim)
%       opts            = model procedure step struct (see Osprey_gLCM)
%
%   OUTPUTS:
%       obj     = OspreyFitObj instance.
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
%%  Define OspreyFitObj properties  

    properties
        % Everything we need to perform a fit and store the results.
        step = [];                                                          % Initialize current step
        Data = struct('fids', [], ...                                       % Initialize time domain data
                      'DwellTime', [] , ...                                 % Initialize dwell time
                      'SpectralWidth', [], ...                              % Initialize spectral width
                      'txfrq', [], ...                                      % Initialize receiver frequency
                      't', [], ...                                          % Initialize time vector
                      'ppm', [], ...                                        % Initialize frequency axis in ppm
                      'nucleus', []);                                       % Initialize nucleus string, e.g. '1H'
        BasisSets = struct('fids', [], ...                                  % Initialize time domain basis functions
                           'names', [], ...                                 % Initialize basis function names
                           'includeInFit', []);                             % Initialize index vector of basis functions to be included
        BaselineBasis = struct('specs', []);                                % Initialize basline basis
        NormNoise = [];                                                     % Initialize 1-point noise estimate 
        Options = {struct};                                                 % Initialize options array
        Model = {struct};                                                   % Initialize model array
    end

%% Define OspreyFitObj nethods  
    methods      
        function obj = FitObject(data, basis, options)
            % class constructor
            if(nargin > 0) % don't have to initialize the fields
                
                obj.step                = 0;                                % Set current step to 1
                %%% DATA %%%
                % Copy the information necessary to create appropriate time-
                % and frequency-domain signals:
                obj.Data.fids            = data.fids;                       % Set time domain data
                obj.Data.DwellTime       = data.dwelltime(1);               % Set dwell time of first extra entry
                obj.Data.SpectralWidth   = data.spectralwidth(1);           % Set spectral width of first extra entry
                obj.Data.txfrq           = data.txfrq(1);                   % Set receiver frequency of first extra entry
                obj.Data.t               = data.t;                          % Set time vector
                obj.Data.nucleus         = data.nucleus;                    % Set nucleus string
                obj.NormNoise            = osp_gLCM_getOnePointNoise(data); % Estimate 1-point noise for normalization
                nptsData        = size(data.fids, 1);                       % Get number of datapoints
                obj.Data.ppm    = calculatePPMAxis(nptsData, data.spectralwidth(1), data.txfrq(1), data.nucleus);   % Create frequency axis in ppm
            
                %%% BASIS SET %%%
                basis.nucleus = data.nucleus;                               % Assume that the basis set nucleus matches the data nucleus                             
                basis.txfrq  = basis.Bo * lookUpGyromagRatioForNucleus(basis.nucleus) * 1e6;  % Calculate the receiver frequency

                % Check that basis set and data have the same resolution
                if abs(basis.dwelltime - obj.Data.DwellTime) > eps
                    warning('Dwell time does not agree between basis set (%5.2e) and data (%5.2e).', ...
                            obj.Data.DwellTime, basis.dwelltime);
                    fprintf('Resampling the basis set for you. \n');
                    basis = fit_resampleBasis(data, basis);                 % Resample basis set
                end
                if round(basis.spectralwidth) ~= round(obj.Data.SpectralWidth)
                    warning('Spectral width does not agree between basis set (%5.2e) and data (%5.2e).', ...
                            obj.Data.DwellTime, basis.spectralwidth);
                    fprintf('Resampling the basis set for you. \n');
                    basis = fit_resampleBasis(data, basis);                 % Resample basis set
                end
                              
                obj.BasisSets.fids  = basis.fids;                           % Set time domain basis functions
                obj.BasisSets.names = basis.name;                           % Set basis function names
                obj.BasisSets.includeInFit = ones(size(basis.name));        % Set index vector of basis functions to be included
                            
                
                %%% OPTIONS %%%
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
                
                obj.Options{1} = options;                                   % Save the property struct
                fitRange    = obj.Options{1}.optimFreqFitRange;              % Get fit range

                % Setup baseline model 
                switch obj.Options{1}.baseline.type                         % Switch for baseline type
                    case 'spline'
                        %%% CREATE BASELINE SPLINE BASIS %%%
                        % Combine real and imaginary part to form a complex spline array.
                        % Use the new, corrected function from here on                      
                        dkntmn      = obj.Options{1}.baseline.dkntmn;       % Get spline basis knot spacing
                        [splineArray] = osp_gLCM_makeSplineBasis(data, fitRange, dkntmn);   % Create spline baseline basis array
                        obj.BaselineBasis = splineArray;                    % Store baseline array in object
                    case 'poly'
                        %%% CREATE BASELINE POLYNOMIAL BASIS %%%
                        order      = obj.Options{1}.baseline.order;         % Get order of the polynomial baseline
                        [splineArray] = osp_gLCM_makePolyBasis(data, fitRange, order);  % Create polynomial baseline basis array    
                        obj.BaselineBasis = splineArray;                    % Store baseline array in object
                    case 'none'
                        %%% NO BASELINE %%%
                        obj.BaselineBasis = [];                             % Store empty baseline array in object 
                end
            end     
        end                                                    
    end
        
    % Static methods, helper functions
    methods (Static)
  
    end    
end