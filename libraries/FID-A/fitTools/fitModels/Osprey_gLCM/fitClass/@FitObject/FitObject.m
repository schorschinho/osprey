classdef FitObject < handle
    
    
    properties
        % Everything we need to perform a fit and store the results.
        step = [];
        Data = struct('fids', [], 'DwellTime', [] , 'SpectralWidth', [], 'txfrq', [], 't', [], 'ppm', []);
        BasisSets = struct('fids', [], 'names', [], 'includeInFit', []);
        BaselineBasis = struct('specs', []);
        Options = {struct};
        Model = {struct};
        Results = struct;
    end
    
    
    
    methods
        
        function obj = FitObject(data, basis, options)
            % class constructor
            if(nargin > 0) % don't have to initialize the fields
                
                obj.step                = 0;
                %%% DATA %%%
                % Copy the information necessary to create appropriate time-
                % and frequency-domain signals:
                obj.Data.fids            = data.fids;
                obj.Data.DwellTime       = data.dwelltime;
                obj.Data.SpectralWidth   = data.spectralwidth;
                obj.Data.txfrq           = data.txfrq;
                obj.Data.t               = data.t;
                
                % Calculate the ppm axis
                npts            = size(data.fids, 1);
                spectralwidth   = data.spectralwidth;
                f   = (-spectralwidth/2) + (spectralwidth/(2*npts)) : spectralwidth/npts : (spectralwidth/2) - (spectralwidth/(2*npts));
                obj.Data.ppm = f / (data.txfrq * 1e-6) + 4.68;
            
                
                %%% BASIS SET %%%
                % Check that basis set and data have the same resolution
                if abs(basis.dwelltime - obj.Data.DwellTime) > eps
                    warning('Dwell time does not agree between basis set (%5.2e) and data (%5.2e).', obj.Data.DwellTime, basis.dwelltime);
                    fprintf('Resampling the basis set for you. \n');
                    basis = fit_resampleBasis(data, basis);
                end
                if round(basis.spectralwidth) ~= round(obj.Data.SpectralWidth)
                    warning('Spectral width does not agree between basis set (%5.2e) and data (%5.2e).', obj.Data.DwellTime, basis.spectralwidth);
                    fprintf('Resampling the basis set for you. \n');
                    basis = fit_resampleBasis(data, basis);
                end
                
                % Copy the necessary information
                obj.BasisSets.fids  = basis.fids;
                obj.BasisSets.names = basis.name;
                obj.BasisSets.includeInFit = ones(size(basis.name));
                
                
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
                
                % Save the property struct
                obj.Options{1} = options;
                
                
                %%% CREATE BASELINE SPLINE BASIS %%%
                % Combine real and imaginary part to form a complex spline array.
                % Use the new, corrected function from here on
                fitRange    = obj.Options{1}.optimFreqFitRange;
                dkntmn      = obj.Options{1}.dkntmn;
                [splineArray] = osp_gLCM_makeSplineBasis(data, fitRange, dkntmn);
                B = splineArray;

                paddedBaseline = B;
                
                % Save into the property
                obj.BaselineBasis = paddedBaseline;
            
            end
            
        end     
            
        
                                                
    end
    
    
    
    % Static methods, helper functions
    methods (Static)

        sse = initialFitLossFunction(x, tdData, basisSet, baselineBasis, ppm, t, fitRange)
        
        specs = transformBasis(fids, gaussLB, lorentzLB, freqShift, t)        
        
    end
    

end