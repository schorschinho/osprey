function [basisSet]=osp_recalculate_basis_specs(basisSet)
%% [basisSet]=osp_recalculate_basis_specs(basisSet)
% This function recalculates the basis spectra and ppm-axis of the basis
% set. 
%   USAGE:
%       [basisSet]=osp_recalculate_basis_specs(basisSet);
%
%   INPUTS:
%       basisSet     = Basis set loaded in FID format.
%
%   OUTPUTS:
%       MRSCont     = Basis set, now including frequency domain spectra and
%                     ppm axis
%
%   AUTHOR:
%       Helge Zoellner, Johns Hopkins University 2023.

basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% Calcualte ppm-axis
f = [(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
basisSet.ppm = f/(basisSet.Bo*42.577);
basisSet.ppm=basisSet.ppm + basisSet.centerFreq;
end