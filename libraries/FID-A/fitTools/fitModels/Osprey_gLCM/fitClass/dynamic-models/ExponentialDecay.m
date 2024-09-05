%%  ExponentialDecay
%   This is a Osprey dynamic-models function. This function can be
%   specified in the DynamicModelJson field in the model procedure json. It
%   has to have a forward model and a jacobian with the reparameterized
%   model relation along the indirect dimension
%
%   USAGE:
%       Specify this in the extra.DynamicModelJson field in the model
%       procedure
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
%% Handle set up for optimizer

function fh = ExponentialDecay
    fh.fun = @ExpDecay;         % Forward model
    fh.jac = @ExpDecayJac;      % Jacobian
end

%% Dynamic forward models and jacobian functions

function prediction = ExpDecay(x, t)
% This describes an exponential decay releation of a parameter
% defined by an inital amplitude and a decay constant. This can be an
% amplitude in a TE-Series or diffusion behaviour. Can be applied to any
% parameter for example metabolite or baseline amplitudes.
%   USAGE:
%       prediction = ExpDecay(x, t)
%
%   INPUTS:
%       x          = parameter matrix
%       t          = modulator for example time
%
%   OUTPUTS:
%       prediction = evolution along modualtor axis
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
%%  Calculate forward model
 
    Ampl     = x(1,:);                      % Define amplitude from 2D vector
    Decay    = x(2,:);                      % Define decay from 2D vector
    prediction = Ampl .* exp(-repmat(t',[1 size(Ampl,2)])./Decay);    % Calculate amplitude evolution
    
end

function [jac] = ExpDecayJac(x, t)
% This describes an exponential decay releation of a parameter
% defined by an inital amplitude and a decay constant. This can be an
% amplitude in a TE-Series or diffusion behaviour. Can be applied to any
% parameter for example metabolite or baseline amplitudes.
%   USAGE:
%       jac        = ExpDecayJac(x, t)
%
%   INPUTS:
%       x          = parameter matrix
%       t          = modulator for example time
%
%   OUTPUTS:
%       jac        = jacobian
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
%%  Calculate jacobian
    Ampl     = x(1,:);                                          % Define amplitude from 2D vector
    Decay    = x(2,:);                                          % Define decay from 2D vector  
   
    dYdAmpl = exp(-repmat(t',[1 size(Ampl,2)])./Decay);                                   % Partial derivative wrt Ampl
    dYdDecay = repmat(t',[1 size(Ampl,2)]) .* Ampl .* (1./(Decay.^2)) .* exp(-repmat(t',[1 size(Ampl,2)])./Decay);  % Partial derivative wrt Decay

    jac = cat(3,dYdAmpl,dYdDecay);                              % Concatenate lines to generate jacobian
    jac = squeeze(jac);                                         % Remove zero dimensions
end

