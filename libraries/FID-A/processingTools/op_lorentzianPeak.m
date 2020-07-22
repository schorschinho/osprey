% op_lorentzianPeak.m
% Helge Zollner, Johns Hopkins University 2020.
% 
% USAGE:
% out=op_lorentzianPeak(n,sw,Bo,centerFreq,lw,ppm0,amp);
% 
% DESCRIPTION:
% Generate a noiseless spectrum containing a single lorentzian peak with
% desired parameters (frequency, amplitude, linewidth, etc.).
% 
% INPUTS:
% n             = Number of points in spectrum.
% sw            = spectral width of spectrum (Hz).
% Bo            = Magnetic field strength (Tesla).
% centerFreq    = center frequency of the spectrum (ppm).
% lw            = Linewidth of gaussian peak (Hz).
% ppm0          = Frequency of gaussian peak (ppm).
% amp           = Amplitude of gaussian peak.  
%
% OUTPUTS:
% out    = Lorentzian lineshape peak in FID-A structure format 


function out=op_lorentzianPeak(n,sw,Bo,centerFreq,lw,ppm0,amp);

dt=1/sw;
df=sw/n;
f0=(-centerFreq+ppm0)*(Bo*42.577)* 2 * pi;

t=[0:dt:(n-1)*dt];
f=[(-sw/2+(df/2)):df:(sw/2-(df/2))];

ppm=f/(Bo*42.577);
ppm=ppm+4.68;
ppm=ppm-(4.68-centerFreq);

fids=amp * exp(-t*pi*lw) .* exp(-1j*f0*t);
fids=fids';

specs=fftshift(fft(fids,[],1),1);

out.t=t;
out.fids=fids;
out.ppm=ppm;
out.specs=specs;

%FILLING IN DATA STRUCTURE

out.sz=size(out.fids);
out.spectralwidth=sw;
out.dwelltime=dt;
out.txfrq=Bo*42.577;
out.centerFreq=centerFreq;
out.date=date;
out.dims.t=1;
out.dims.averages=0;
out.dims.subSpecs=0;
out.dims.coils=0;
out.Bo=Bo;
out.averages=1;
out.rawAverages=1;
out.subspecs=1;
out.rawSubspecs=1;

%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=0;
out.flags.addedrcvrs=0;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
if out.dims.subSpecs==0
    out.flags.isISIS=0;
else
    out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
end