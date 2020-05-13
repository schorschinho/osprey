function [resBasisSet] = fit_resampleBasis(dataToFit, basisSet)
%% [resBasisSet] = fit_resampleBasis(dataToFit, basisSet)
%   This function downsamples basis functions (which are usually given in
%   higher resolution than the data that is to be fitted) to match the
%   resolution and frequency range of the data.
%
%   USAGE:
%       [resBasisSet] = fit_resampleBasis(dataToFit, basisSet);
%
%   INPUTS:
%       dataToFit   = Individual FID-A data container.
%       basisSet    = FID-A basis set container.
%
%   OUTPUTS:
%       resBasisSet = FID-A basis set container containing resampled basis functions.
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-09-05)
%       goeltzs1@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2019-09-05: First version of the code.


% Determine the ppm ranges of both the data and the basis functions.
ppmRangeData        = dataToFit.ppm';
ppmRangeBasis       = basisSet.ppm;
ppmIsInDataRange    = (ppmRangeBasis < ppmRangeData(1)) & (ppmRangeBasis > ppmRangeData(end));
if sum(ppmIsInDataRange) == 0
    ppmIsInDataRange    = (ppmRangeBasis > ppmRangeData(1)) & (ppmRangeBasis < ppmRangeData(end));
end

% Now resample the basis functions to match the resolution and frequency
% range (ppm axis) of the data.
fids_interp     = zeros(length(ppmRangeData), (basisSet.nMets + basisSet.nMM), size(basisSet.fids, 3)); % allocate memory
specs_interp    = zeros(length(ppmRangeData), (basisSet.nMets + basisSet.nMM), size(basisSet.fids, 3)); % allocate memory
% Loop over the number of basis functions (i.e. metabolites and
% MM/lipids)
for ll=1:(basisSet.nMets + basisSet.nMM)
    % Loop over the number of sub-spectra (MEGA, HERMES, HERCULES)
    for rr = 1:size(basisSet.fids, 3)
        specs_interp(:,ll,rr)      = interp1(ppmRangeBasis(ppmIsInDataRange), basisSet.specs(ppmIsInDataRange,ll,rr), ppmRangeData, 'pchip', 'extrap');
        %convert back to time domain
        %if the length of Fids is odd, then you have to do a circshift of one to
        %make sure that you don't introduce a small frequency shift into the fids
        %vector.
        if mod(length(basisSet.specs),2)==0
            %disp('Length of vector is even.  Doing normal conversion');
            fids_interp(:,ll,rr)   = ifft(fftshift(specs_interp(:,ll,rr),1),[],1);
        else
            %disp('Length of vector is odd.  Doing circshift by 1');
            fids_interp(:,ll,rr)   = ifft(circshift(fftshift(specs_interp(:,ll,rr),1)),[],1);
        end
    end
end

% Create the output resampled basis set container
resBasisSet         = basisSet;
resBasisSet.ppm     = ppmRangeData;
resBasisSet.specs   = specs_interp;
resBasisSet.fids    = fids_interp;
resBasisSet.sz      = size(resBasisSet.fids);
resBasisSet.n       = length(resBasisSet.fids);

% Calculate the new spectral width and dwelltime:
dppm                        = abs(resBasisSet.ppm(2)-resBasisSet.ppm(1));
ppmrange                    = abs((resBasisSet.ppm(end)-resBasisSet.ppm(1)))+dppm;
resBasisSet.spectralwidth   = ppmrange*dataToFit.Bo*42.577;
resBasisSet.dwelltime       = 1/resBasisSet.spectralwidth;
% Calculate the time scale
resBasisSet.t = (0:resBasisSet.dwelltime:(resBasisSet.sz(1)-1)*resBasisSet.dwelltime);


end