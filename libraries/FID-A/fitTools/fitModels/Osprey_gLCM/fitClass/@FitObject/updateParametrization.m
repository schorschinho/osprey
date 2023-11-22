function updateParametrization(obj, parametrization)
%%  updateParametrization(obj, parametrization)
%   This method writes a parametrization into a FitObject Options property.
%
%   USAGE:
%       obj.updateParametrization(obj, parametrization)
%
%   INPUTS:
%       parametrization = struct describing parametrization
%       
%   OUTPUTS:
%       obj     = OspreyFitObj.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2023-10-18)
%       goeltzs1@jhmi.edu

obj.Options{obj.step+1}.parametrizations = parametrization;

end