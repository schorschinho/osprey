% op_addNoise.m
% Helge Zoellner, Johns Hopkins University 2023.
% 
% USAGE:
% [out,noise]=op_addNoise(in,sdnoise,noise);
% 
% DESCRIPTION:
% Add correalted noise to a spectrum (useful for generating simulated data).  Normally
% distributed random noise is added to both the real and imaginary parts of
% the data.  Real and imaginary noise parts are uncorrelated.
% 
% INPUTS:
% in         = Input data in matlab structure format.
% sdnoise    = Standard deviation of the random noise to be added in the time domain.
% noise      = (optional)  Specific noise kernel to be added (if specified,
%               sdnoise variable is ignored).  
%
% OUTPUTS:
% out        = Output dataset with noise added.
% noise      = The noise vector that was added. 

function [out,noise]=op_addCorrelatedNoise(in,sdnoise,noise);

if nargin<3
    ampl = sqrt(2)*sdnoise*normrnd(0,1,[1, in.sz]); % Normal distribution
    phase = -pi + (pi + pi) * rand(1,in.sz); % Uniform distribution
    noise = amp .* exp(1i .* phase);
end

fids=in.fids + noise;

%re-calculate Specs using fft
specs=fftshift(fft(fids,[],in.dims.t),in.dims.t);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
end