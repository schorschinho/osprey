function ppm = calculatePPMAxis(npts, spectralwidth, txfrq, nucleus)
%% ppm = calculatePPMAxis(npts, spectralwidth, txfrq, nucleus)
% Calculates the ppm axis for a given set of: npnts, spectralwidth,
% txfrq, and nucleus
%
%   USAGE:
%       ppm = calculatePPMAxis(npts, spectralwidth, txfrq, nucleus)
%
%   INPUTS:
%       npts          = number of data points in an FID [1]
%       spectralwidth = spectral width of the FID [Hz]
%       txfrq         = transmitter frequency [Hz]
%       nucleus       = string describing the nucleus
%
%   OUTPUTS:
%       ppm          = ppm axis
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu

%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017) 
%%  Calculate ppm axis

    centerFreq = lookUpCenterFreqForNucleus(nucleus);   % Get center frequency

f   = (-spectralwidth/2) + (spectralwidth/(2*npts)) : spectralwidth/npts : (spectralwidth/2) - (spectralwidth/(2*npts)); % Get axis in Hz
ppm = f / (txfrq * 1e-6) + centerFreq;  % Convert to ppm and shift by center frequency 

end
