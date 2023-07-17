function obj = updateBaselineAccordingToStep(obj, options)
%%  updateBaselineAccordingToStep(obj)
%   This is the method updates the baseline basis according to the model
%   procedure step. This is used if for example the number of splines
%   changes for example to account for a new fit range.
%
%   USAGE:
%       obj.optimizeRegularization()
%
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
%%  Get object details
    
    data        = obj.Data;                                                 % Get data entry
    fitRange    = obj.Options{1,obj.step+1}.optimFreqFitRange;          % Get fit range
    step        = obj.step + 1;
%%  Update baseline basis
    if (obj.Options{1,step}.optimFreqFitRange(1) ~=  obj.Options{1,step}.optimFreqFitRange(1)) || ... % Did the fit range change?
       (obj.Options{1,step}.optimFreqFitRange(2) ~=  obj.Options{1,step}.optimFreqFitRange(2))    
        % Update baseline model 
        switch obj.Options{1,obj.step+1}.baseline.type                      % Switch for baseline type
            case 'spline'                                                   
                %%% CREATE BASELINE SPLINE BASIS %%%
                % Combine real and imaginary part to form a complex spline array.
                % Use the new, corrected function from here on                  
                dkntmn      = obj.Options{1,obj.step+1}.baseline.dkntmn;            % Get spline basis knot spacing
                [splineArray] = osp_gLCM_makeSplineBasis(data, fitRange, dkntmn);   % Create spline baseline basis array     
                obj.BaselineBasis = splineArray;                                    % Store baseline array in object
            case 'poly'
                %%% CREATE BASELINE POLYNOMIAL BASIS %%%
                order      = obj.Options{1,obj.step+1}.baseline.order;              % Get order of the polynomial baseline
                [splineArray] = osp_gLCM_makePolyBasis(data, fitRange, order);      % Create polynomial baseline basis array       
                obj.BaselineBasis = splineArray;                                    % Store baseline array in object
            case 'none'
                %%% NO BASELINE %%%
                obj.BaselineBasis = [];                                             % Store empty baseline array in object
        end
    end
end