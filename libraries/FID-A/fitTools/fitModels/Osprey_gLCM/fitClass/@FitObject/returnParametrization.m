function [value] = returnParametrization(obj,step,field)
%%  plotBasisSet(obj, step, plotRange)
%   This method generates a plot of the basis set 
%
%   USAGE:
%       obj.plotBaselineBasis(step, plotRange)
%
%   INPUTS:
%       step            = step to plot    
%       field           = parameterization to return, e.g. init
%
%   OUTPUTS:
%       figure
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
%%  Diverge to default options if required 

    if nargin < 3
        field = 'init';
        if nargin < 2
            step = obj.step;
        end
    end
    
%%  Get parametrization values
    pars = fields(obj(step).Options{1}.parametrizations);                   % Parameter names              
    value ={};                                                              % Initialize value cell arry
    for ff = 1 : length(pars)                                               % Loop over parameters
        value{end+1} = obj(step).Options{1}.parametrizations.(pars{ff}).(field);    % Get value 
    end                                                                     % End loop over parametes
end