function ppm = calculatePPMAxis(npts, spectralwidth, txfrq, nucleus)
% Calculates the ppm axis for a given set of:
% ppm           = ppm axis
% npts          = number of data points in an FID [1]
% spectralwidth = spectral width of the FID [Hz]
% txfrq         = transmitter frequency [Hz]
% nucleus       = string describing the nucleus

centerFreq = lookUpCenterFreqForNucleus(nucleus);

f   = (-spectralwidth/2) + (spectralwidth/(2*npts)) : spectralwidth/npts : (spectralwidth/2) - (spectralwidth/(2*npts));
ppm = f / (txfrq * 1e-6) + centerFreq;

end
