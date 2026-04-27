%%  OppositeSign
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

function fh = OppositeSign
    fh.fun = @OppSign;         % Forward model
    fh.jac = @OppSignJac;      % Jacobian
end

%% Dynamic forward models and jacobian functions

function prediction = OppSign(x,m)
% This describes function uses opposite sign between a shared parameter. For example,
% Ph0 can be positive in the first spectrum and the exact opposite in the other.
%   USAGE:
%       prediction = OppSign(x)
%
%   INPUTS:
%       x          = parameter matrix
%
%   OUTPUTS:
%       prediction = evolution along modualtor axis
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2026-04-26)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%%  Calculate forward model

    prediction = x;
    x(2) = -x(2);     % Apply opposite sign

end

function [jac] = OppSignJac(x,m)
% This describes function uses opposite sign between a shared parameter. For example,
% Ph0 can be positive in the first spectrum and the exact opposite in the other.
%   USAGE:
%       jac        = OppSignJac(x)
%
%   INPUTS:
%       x          = parameter matrix
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
    dYdPar = [1,-1];

    jac = dYdPar;
end
