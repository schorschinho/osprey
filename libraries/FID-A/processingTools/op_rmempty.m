% op_rmempty.m
% Georg Oeltzschner, Johns Hopkins University 2018.
% 
% USAGE:
% [out]=op_rmempty(in);
% 
% DESCRIPTION:
% Removes empty FIDs (for example in Philips *_ref.sdat files).
% 
% INPUTS:
% in        = input data in matlab structure format.
%
% OUTPUTS:
% out       = Output dataset


function [out]=op_rmempty(in)

if in.flags.averaged
    error('ERROR:  Spectra have already been averaged!  Aborting!');
end

% If the standard deviation of an FID is ~0, this FID row is empty.
% Find the indices of the averages with non-zero standard deviations -
% these are the ones that will be 
idx_nonzeroFID = find(std(in.fids,0,1)>1e-40);
nonzero_fids = in.fids(:,idx_nonzeroFID);

%re-calculate Specs using fft
nonzero_specs=fftshift(fft(nonzero_fids,[],in.dims.t),in.dims.t);

%recalculate the sz vector
sz=size(nonzero_fids);

% 
%FILLING IN DATA STRUCTURES
out=in;
out.fids=nonzero_fids;
out.specs=nonzero_specs;
out.sz=sz;
out.averages=sz(2);

%FILLING IN THE FLAGS
out.flags=in.flags;