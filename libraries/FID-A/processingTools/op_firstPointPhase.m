% op_firstPointPhase.m
% Georg Oeltzschner, Johns Hopkins University 2019.
% 
% USAGE:
% [out]=op_firstPointPhase(in);
% 
% DESCRIPTION:
% Performs a zero-order phase correction by multiplying each point of the
% FID with the normalized complex conjugate of the first point.
% 
% INPUTS:
% in         = input data in matlab structure format.  
%
% OUTPUTS:
% out        = Output following first-point phasing.

function [out]=op_firstPointPhase(in);

% Duplicate input structure
out = in;

% Determine phase correction
fids    = in.fids;
sz      = in.sz;
corrph = conj(fids(1,:))./abs(fids(1,:));
corrph = repmat(corrph, [sz(1) 1]);
if length(sz) == 2
    corrph = reshape(corrph, [sz(1) sz(2)]);
elseif length(sz) == 3
    corrph = reshape(corrph, [sz(1) sz(2) sz(3)]);
end

% Apply phase correction
fids = fids .* corrph;

% Produce specs
specs = fftshift(fft(fids,[],in.dims.t),in.dims.t);

%FILLING IN DATA STRUCTURE
out.specs   = specs;
out.fids    = fids;




