% op_creChoFit.m
% Georg Oeltzschner, Johns Hopkins University 2019.
% 
% USAGE:
% parsFit=op_creChoFit(in, suppressPlot);
% 
% DESCRIPTION:
% Perform a Lorentzian lineshape fit to the creatine resonance in a brain 
% proton MRS dataset.
% 
% INPUTS:
% in            = input data in matlab stucture format.
% suppressPlot  = (optional) Boolian to suppress plots.  Default = 0;
%
% OUTPUTS:
% parsFit  = Fit parameters for the Creatine peak fit.  
%               parsFit(1) = Amplitude (in arbitrary units);
%               parsFit(2) = Linewidth (in Hz);
%               parsFit(3) = Frequency (in ppm);
%               parsFit(4) = Baseline slope;
%               parsFit(5) = Baseline DC Offset;

function parsFit=op_creChoFit(in, suppressPlot);

if nargin < 2
    echo = 0;
end

if in.flags.isISIS
    error('ERROR:  must have combined subspecs in order to do this!  ABORTING');
end

if ~in.flags.averaged
    error('ERROR:  must have averaged in order to do this!  ABORTING');
end

if ~in.flags.addedrcvrs
    error('ERROR:  must have added receivers in order to do this!  ABORTING');
end

% Extract the frequency axis and spectrum
ppm = in.ppm;
spec = in.specs;

% Define the spectral range around the Cr-Cho peaks
freqLim = ppm <= 3.02+0.15 & ppm >= 3.02-0.15;
[~,i] = max(abs(spec(freqLim)));
freq2 = ppm(freqLim);
maxFreq = freq2(i);
freqLim = ppm <= maxFreq+0.58 & ppm >= maxFreq-0.42;
specRange = real(spec(freqLim));

% Create initial guesses for the fit parameters
Baseline = (specRange(1) + specRange(end))/2;
Width = 0.05;
Area = (max(specRange) - min(specRange)) * Width * 4;
x0 = [Area Width maxFreq 0 Baseline 0 1];

% Set up the least-squares and non-linear solvers
lsqopts = optimset('lsqcurvefit');
lsqopts = optimset(lsqopts,'MaxIter',800,'TolX',1e-4,'TolFun',1e-4,'Display','off');
lb = [0 x0(2)*0.00001 x0(3)-0.42];
ub = [x0(1)*10000 x0(2)*2 x0(3)+0.58];
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts,'MaxIter',400,'TolX',1e-6,'TolFun',1e-6);

% Set up the fit problem and solve it
x0 = lsqcurvefit(@TwoLorentzModel, x0, ppm(freqLim)', real(specRange), lb, ub, lsqopts);
[parsFit, residCr] = nlinfit(ppm(freqLim)', specRange, @TwoLorentzModel, x0, nlinopts);

% Optional plotting of data and results
if ~suppressPlot
    % Plot original data
    figure(101);
    plot(ppm(freqLim),specRange);
    hold on
    % Plot the fit
    plot(ppm(freqLim),TwoLorentzModel(parsFit,ppm(freqLim)));
    % Plot the fit with phase 0
    parsFit2 = parsFit;
    parsFit2(4) = 0;
    plot(ppm(freqLim),TwoLorentzModel(parsFit2,ppm(freqLim)));
    legend('data','fit','phased fit');
    hold off;
end

end

function Lorentz = TwoLorentzModel(x,freq)
% CJE LorentzModel
% Lorentzian  = (1/pi) * (hwhm) / (deltaf^2 + hwhm^2) (Wolfram)
% Peak height of Lorentzian = 4 / (pi*hwhm)
% This defnition of the Lorentzian has Area = 1
%
% This model is adapted from Gannet3.1

area = x(1);
hwhm = x(2);
f0 = x(3);
phase = x(4);
baseline0 = x(5);
baseline1 = x(6);

Absorption = 1/(2*pi) * area * ones(size(freq)) .* hwhm ./ ((freq-f0).^2 + hwhm.^2)+ ...
    1/(2*pi) * area * x(7) * ones(size(freq)) .* hwhm ./ ((freq-f0-0.18).^2 + hwhm.^2);

Dispersion = 1/(2*pi) * area * (freq-f0) ./ ((freq-f0).^2 + hwhm.^2) + ...
    1/(2*pi) * area * x(7) * (freq-f0-0.18) ./ ((freq-f0-0.18).^2 + hwhm.^2);

Lorentz = Absorption*cos(phase) + Dispersion*sin(phase) ...
    + baseline0 + baseline1*(freq-f0);
end
