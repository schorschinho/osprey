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

pars = fields(parametrization);                                             % Get parameter names to update
for ff = 1 : length(pars) 
    vals = fields(parametrization.(pars{ff}));
    for vv = 1 : length(vals)
        if ~isfield(obj.Options{obj.step+1}.parametrizations,(pars{ff}))
            obj.Options{obj.step+1}.parametrizations.(pars{ff}) = [];
        end
        if ~isfield(obj.Options{obj.step+1}.parametrizations.(pars{ff}),(vals{vv}))
            obj.Options{obj.step+1}.parametrizations.(pars{ff}).(vals{vv}) = parametrization.(pars{ff}).(vals{vv});
        end        
    end
end
end