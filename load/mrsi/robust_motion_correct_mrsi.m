function [spec_fs, spec_phs, spec_w] = robust_motion_correct_mrsi(k_spec,first,Larmor)

spec_w = ones(size(k_spec,1),size(k_spec,2),size(k_spec,3));
k_spec_shift = zeros(size(k_spec));
refShift = zeros(size(k_spec,1),size(k_spec,2),size(k_spec,3),size(k_spec,5),size(k_spec,6));
SimMetricOn = ones(size(k_spec,1),size(k_spec,2),size(k_spec,3),size(k_spec,5))*nan;
SimMetricOff = ones(size(k_spec,1),size(k_spec,2),size(k_spec,3),size(k_spec,5))*nan;

dims = size(k_spec);
n_points = max(dims);
    
ppm_axis = linspace(-1000,1000,n_points)/Larmor + 4.68;


% figure, plot(ppm_axis,squeeze(abs(k_spec(first(1),first(2),first(3),:,1,1)))), hold on
% plot(ppm_axis,squeeze(abs(k_spec(first(1),first(2),first(3),:,2,1))))
% plot(ppm_axis,squeeze(abs(k_spec(first(1),first(2),first(3),:,1,2))))
% plot(ppm_axis,squeeze(abs(k_spec(first(1),first(2),first(3),:,2,2))))

%Frequency Shift
for kz = 1 : size(k_spec,1)
    for kx = 1 : size(k_spec,2)
        for ky = 1 : size(k_spec,3)
            for av = 1 : size(k_spec,5)
                for ss = 1 : size(k_spec,6)
                    if sum(squeeze(k_spec(kz,kx,ky,:,av,ss)))>0
                        [refShiftPoints(kz,kx,ky,av,ss)] = CrChoReferencingKSpace(squeeze(k_spec(kz,kx,ky,:,av,ss)),ppm_axis);  
                        refShiftPointsCirc(kz,kx,ky,av,ss) = -round(refShiftPoints(kz,kx,ky,av,ss));
                        k_spec_shift(kz,kx,ky,:,av,ss) = circshift(squeeze(k_spec(kz,kx,ky,:,av,ss)),refShiftPointsCirc(kz,kx,ky,av,ss));
                    end
                end
            end
        end
    end
end

spec_fs = refShiftPoints * (ppm_axis(1)-ppm_axis(2)) *Larmor;

%Similarity Metric
OnFirst = mean(squeeze(k_spec_shift(first(1),first(2),first(3),:,:,1)),2);
OffFirst = mean(squeeze(k_spec_shift(first(1),first(2),first(3),:,:,2)),2);
for kz = 1 : size(k_spec,1)
    for kx = 1 : size(k_spec,2)
        for ky = 1 : size(k_spec,3)
            for av = 1 : size(k_spec,5)
                if sum(squeeze(k_spec(kz,kx,ky,:,av,1)))>0                    
                    SimMetricOff(kz,kx,ky,av) = sum((abs(OffFirst(ppm_axis <= 6 & ppm_axis >= 1.95)) - squeeze(abs(k_spec_shift(kz,kx,ky,ppm_axis <= 6 & ppm_axis >= 1.95,av,2)))).^2) / (length(OffFirst));
                end
                 if sum(squeeze(k_spec(kz,kx,ky,:,av,2)))>0                    
                    SimMetricOn(kz,kx,ky,av) = sum((abs(OnFirst(ppm_axis <= 6 & ppm_axis >= 1.95)) - squeeze(abs(k_spec_shift(kz,kx,ky,ppm_axis <= 6 & ppm_axis >= 1.95,av,1)))).^2) / (length(OnFirst));
                end
            end
        end
    end
end

% Calculate weights  
SimMetricOff(isnan(SimMetricOff)) = 1;
WOff = 1./SimMetricOff.^2;        
WOff = WOff./sum(WOff(:,:,:,:),4);
SimMetricOn(isnan(SimMetricOn)) = 1;
WOn = 1./SimMetricOn.^2;        
WOn = WOn./sum(WOn(:,:,:,:),4);

for kz = 1 : size(k_spec,1)
    for kx = 1 : size(k_spec,2)
        for ky = 1 : size(k_spec,3)
            spec_w(kz,kx,ky) = min([WOff(kz,kx,ky,1),WOff(kz,kx,ky,2), WOn(kz,kx,ky,1), WOn(kz,kx,ky,2)]);
        end
    end
end

% Calculate phases
spec_phs = zeros(size(k_spec,1),size(k_spec,2),size(k_spec,3),size(k_spec,5),size(k_spec,6));
for kz = 1 : size(k_spec,1)
    for kx = 1 : size(k_spec,2)
        for ky = 1 : size(k_spec,3)
            for av = 1 : size(k_spec,5)
                for ss = 1 : size(k_spec,6)
                    noise = std(squeeze(abs(k_spec(kz,kx,ky,ppm_axis <= 11 & ppm_axis >= 10,av,ss))));
                    sig = max(squeeze(abs(k_spec(kz,kx,ky,ppm_axis <= 3.15 & ppm_axis >= 2.9,av,ss))));
                    SNR = sig/noise;
                    if sum(squeeze(k_spec(kz,kx,ky,:,av,ss)))>0 && ~isnan(SNR) && SNR >45
                        temp_spec = squeeze(k_spec(kz,kx,ky,ppm_axis <= 3.15 & ppm_axis >= 2.9,av,ss));
                        ppmindex=find(abs(temp_spec)==max(abs(temp_spec)));
                        spec_phs(kz,kx,ky,av,ss)=-phase(squeeze(temp_spec(ppmindex)))*180/pi; 
                    end
                end
            end
        end
    end
end



end



function [refShift] = CrChoReferencingKSpace(specs,ppm_axis)
%% [refShift, refFWHM] = osp_CrChoReferencing(dataToFit)
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
spec = specs(ppm_axis <= 6 & ppm_axis >= 1.85);
ppm  = ppm_axis(ppm_axis <= 6 & ppm_axis >= 1.85);
powerspec = spec .* conj(spec) / length(spec);
x = ppm;
y = zeros(size(x));

% Set up delta functions
[a,b] = min((abs(x-4.68)));
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
[max_value,max_ind] = max(r(limits));
tempx = newx(limits);
tempr = r(limits);
try
    gtHalfMax=find(tempr >= 0.5 * max_value);
    HWHM=abs(tempx(gtHalfMax(1)) - tempx(gtHalfMax(end)))/2 * 0.8;
    if HWHM == 0
        gtHalfMax=find(tempr >= 0.3 * max_value);
        HWHM=abs(tempx(gtHalfMax(1)) - tempx(gtHalfMax(end)))/2 * 0.8;
    end
    
    nlinopts    = statset('nlinfit');
    % nlinopts    = statset(nlinopts,'MaxIter',1e8,'MaxFunEvals',1e8,'TolX',1e-10,'TolFun',1e-10);
    LorentzModelInit = [1 HWHM tempx(max_ind) 0];


    LorentzModelParam = lsqcurvefit(@LorentzModel,LorentzModelInit,newx(limits),real(r(limits)),[1 0 tempx(1) -180],[1 2*HWHM tempx(end) 180],nlinopts);
    % Return FWHM and reference shift
    refFWHM     = 2 * LorentzModelParam(2) * abs(ppm(1)-ppm(2));
    refShift    = (LorentzModelParam(3) - length(newx)/2);
catch
    refFWHM     = nan;
    refShift    = 0;
end

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