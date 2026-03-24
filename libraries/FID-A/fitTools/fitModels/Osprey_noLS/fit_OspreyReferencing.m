function [refShift, refFWHM] = fit_OspreyReferencing(dataToFit)
%% [refShift, refFWHM] = fit_OspreyReferencing(dataToFit)
%   Calculates the reference shift and full-width half-maximum (FWHM) of a
%   spectrum according to the LCModel algorithm. The algorithm is described in:
%       S.W. Provencher, "Estimation of metabolite concentrations from
%       localized in vivo NMR spectra", Magn Reson Med 30(6):672-679 (1993)
%
%   The algorithm first determines the power spectrum of the input
%   spectrum, before calculating the cross-correlation with a function
%   containing delta functions at 2.01 ppm (NAA), 3.03 ppm (Cr), and 3.22
%   ppm (Cho).
%   The frequency of the largest peak of the cross-correlation function is
%   returned as the initial referencing shift of the input spectrum.
%   The FWHM of the largest peak of the cross-correlation function is 
%   returned as an estimate of the FWHM of the input spectrum.
%
%   Input:
%       dataToFit   = FID-A data structure
%
%   Output:
%       refShift    = Reference frequency shift (in Hz)
%       FWHM        = full-width half-maximum of the input spectrum (in ppm)
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-10-12)
%       goeltzs1@jhmi.edu
%
%   History:
%       2018-10-12: First version of the code.
%

% Calculate power spectrum to estimate reference shift and FWHM
dataToFit=op_freqrange(dataToFit,1.85,4.2);
spec = dataToFit.specs;
ppm  = dataToFit.ppm;
powerspec = spec .* conj(spec) / length(spec);
x = ppm;
y = zeros(size(x));

% Set up delta functions
[a,b] = min((abs(x-2.01)));
[c,d] = min((abs(x-3.03)));
[e,f] = min((abs(x-3.22)));
y(b) = 1;
y(d) = 1;
y(f) = 1;

% Plot cross-correlation function, normalize it to its maximum
r = crosscorr(powerspec,y)';
r2 = crosscorr(real(spec),y)';
r = r./(max(r));
r2 = r2./(max(r2));
% Set up Lorentzian fit to the central peak of the cross-correlation
% function
% x-axis of cross-correlation function
newx = [1:1:(2*length(spec)-1)];

% fit the center peak
limits      = newx >= 0.94*length(spec) & newx <= 1.06*length(spec);
[~,max_ind] = max(r(limits));
tempx = newx(limits);
nlinopts    = statset('nlinfit');
nlinopts    = statset(nlinopts,'MaxIter',1e8,'MaxFunEvals',1e8,'TolX',1e-10,'TolFun',1e-10);
LorentzModelInit = [1 1/pi tempx(max_ind) 0];
[LorentzModelParam, ~] = nlinfit(newx(limits), real(r(limits)), @LorentzModel, LorentzModelInit, nlinopts);

% Return FWHM and reference shift
refFWHM     = 2 * LorentzModelParam(2) * abs(ppm(1)-ppm(2));
refShift    = (LorentzModelParam(3) - length(newx)/2) * (ppm(1)-ppm(2)) * dataToFit.txfrq*1e-6;


%%% embedded Lorentzian model function
function Lorentz = LorentzModel(x,freq)
% CJE LorentzModel
% Lorentzian  = (1/pi) * (hwhm) / (deltaf^2 + hwhm^2) (Wolfram)
% Peak height of Lorentzian = 4 / (pi*hwhm)
% This defnition of the Lorentzian has Area = 1

area = x(1);
hwhm = x(2);
f0 = x(3);
phase = x(4);

Absorption = 1/(2*pi) * area * ones(size(freq)) .* hwhm ./ ((freq-f0).^2 + hwhm.^2);
Dispersion = 1/(2*pi) * area * (freq-f0) ./ ((freq-f0).^2 + hwhm.^2);

Lorentz = Absorption*cos(phase) + Dispersion * sin(phase);
end

end