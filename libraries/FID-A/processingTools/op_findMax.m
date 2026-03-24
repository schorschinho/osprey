%% op_findMax.m
%
%   This function finds maximum value in the FD for all subspectra.
%
%
%   USAGE:
%           out=op_findMax(in,abs);
%
%   INPUT:  
%           in      = data in FID-A struct format
%           abs     = absolute values flag
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-02-08)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
% 

function out=op_findMax(in,absolute);
% Fall back into defaults
if nargin < 2
    absolute =0;
end

% Take absolute of real part of the spectrum if needed
if ~absolute
    specs=real(in.specs);
else
    specs=abs(real(in.specs));
end

    out = max(specs);
    while length(out) > 1
        out = max(out);
    end
end
