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

% GO 6/25/2025
% There have been a few changes in the MRSCloud simulation code that I have
% been using internally and for collaborators (this will be obsolete in the
% eventual Osprey 3.0 with the new model. Basically, the metabolite basis
% functions are shifted in the time domain (by the difference between the
% slice-selection center frequency and the water resonance) - this resolves
% a couple not-so-nice effects with the old definitions when the output
% data were zero-filled.
% However, this new way of doing things needs to be accommodated here
% (since Osprey expects .mat basis sets to be in the 'old' convention).
% I will simply intercept 'new versions' of the MRSCloud sims and shift 
% the metabolites by the same amount.

% If the field basisSet.mrsCloudVersion does not exist (or is v1), we can
% safely assume that this is an 'old-convention' .mat basis set
if isfield(basisSet, 'mrsCloudVersion')
    if strcmpi(basisSet.mrsCloudVersion, 'v2')

        % Shift ppm axis by centreFreq
        txfrq       = basisSet.Bo * 42577000;
        f_shift     = -(4.68-basisSet.centerFreq)*txfrq*1e-6;
        t           = repmat(basisSet.t',[1 basisSet.sz(2:end)]);
        % Only shift the metabolites, not the MMs
        t           = t(:,1:basisSet.nMets,:);
        basisSet.fids(:,1:basisSet.nMets,:) = basisSet.fids(:,1:basisSet.nMets,:).*exp(-1i*t*f_shift*2*pi);

    end
end

basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% Calcualte ppm-axis
f = [(-basisSet.spectralwidth/2)+(basisSet.spectralwidth/(2*basisSet.sz(1))):basisSet.spectralwidth/(basisSet.sz(1)):(basisSet.spectralwidth/2)-(basisSet.spectralwidth/(2*basisSet.sz(1)))];
basisSet.ppm = f/(basisSet.Bo*42.577);
basisSet.ppm=basisSet.ppm + basisSet.centerFreq;
end