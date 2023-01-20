% osp_gLCM_getOnePointNoise.m
% Helge Zoellner, The Johns Hopkins UNiversity 2023.
%
% USAGE:
% [Noise]=osp_gLCM_getNoise(in,noiseppmmin,noiseppmmax);
%
% DESCRIPTION:
% Find the noise in a spectrum normalized to the number of points.
%
% INPUTS:
% in             = specs vector
% noiseppmmin    = min of frequency range in which to measure noise.
%                  (Optional.  Default = -2 ppm);
% noiseppmmax    = max of frequency range in which to measure noise.
%                  (Optional.  Default = 0 ppm);
%
% OUTPUTS:
% Noise            = Estimated noise of the input spectrum.

function [Noise]=osp_gLCM_getOnePointNoise(in,noiseppmmin,noiseppmmax);


if nargin<3
    noiseppmmax=0;
    if nargin<2
        noiseppmmin=-2;
    end
end

%NOW FIND THE STANDARD DEVIATION OF THE NOISE:
noisewindow=in.specs(in.ppm>noiseppmmin & in.ppm<noiseppmmax,1);
ppmwindow2=in.ppm(in.ppm>noiseppmmin & in.ppm<noiseppmmax)';

P=polyfit(ppmwindow2,noisewindow,2);
noise=noisewindow-polyval(P,ppmwindow2);

noisesd=std(real(noise));

Noise = noisesd/length(noisewindow);


