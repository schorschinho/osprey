% fit_OspreyMM2.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% [fitParams] = fit_Osprey(dataToFit, resBasisSet, fitOpts);
% 
% DESCRIPTION:
% THIS RETURNS A CLEANED MM SPECTRUM - need to fix up comments and maybe
% change name to somethign useful

function [out] = fit_OspreyCleanMM(dataToFit, basisSet, fitOpts,scale,fitParams)


%%% 0. PREPARE DATA AND BASIS SET %%%
dataToFit               = op_zeropad(dataToFit, 2);

% switch fitOpts.sequence
%     case 'unedited'
%         %Invert the NAA, Cr and CrCH2 peaks
%         basisSet.specs(:,1)=-basisSet.specs(:,1);
%         basisSet.specs(:,2)=-basisSet.specs(:,2);
%         basisSet.specs(:,4)=-basisSet.specs(:,4);
%     case 'MEGA'
%         %Invert the NAA/NAAG peak
%         basisSet.specs(:,1)=-basisSet.specs(:,1);
%         basisSet.specs(:,2)=-basisSet.specs(:,2);
% end

basisSet             = fit_resampleBasis(dataToFit, basisSet);

out=dataToFit;


%%% 1. EXTRACT OPTIONS AND PREPARE FIT %%%

% ... fit parameters
nMets       = basisSet.nMets;
nMM         = basisSet.nMM;
nBasisFcts  = nMets + nMM; % number of basis functions
lineShape   = fitParams.lineShape;
ph0         = fitParams.ph0 * pi/180; % zero-order phase correction [convert from deg to rad]
ph1         = fitParams.ph1 * pi/180; % first-order phase correction [convert from deg/ppm to rad/ppm]
gaussLB     = fitParams.gaussLB; % Gaussian damping [Hz]
lorentzLB   = fitParams.lorentzLB; % Lorentzian damping [Hz] for each basis function
freqShift   = fitParams.freqShift; % Frequency shift [Hz] for each basis function
ampl        = fitParams.ampl; % Amplitudes for metabolite/MM/lipid basis functions
refShift    = fitParams.refShift; % Reference shift applied to the data during first step of fitting

% Normalize
lineShape = lineShape/sum(lineShape);

%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
% Run the time-domain operations on the metabolite basis functions
% (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
t = basisSet.t;
%ph0=-ph0;
for ii=1:nBasisFcts
    basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp(-1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)';    
     basisSet.fids(:,ii) = basisSet.fids(:,ii) * exp(1i*ph0);
end
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% Run the frequency-domain operations on the basis functions
% (first order phase correction)
% Cut out the frequency range of the basis set
%basisSet = op_freqrange(basisSet,fitRangePPM(1),fitRangePPM(end),length(splineArray(:,1,1)));
% Create a ppm vector around a pivot point (water)
ppm_ax = basisSet.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
%ph1=-ph1;
% Apply the linear phase correction
for ii=1:nBasisFcts
    basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp(1i*ph1*multiplier);
end
basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);


%%% 3. APPLY THE LINEAR PARAMETERS %%%
% Convolve the lineshape with the metabolite basis functions only
% (NOT the macromolecules or lipids or baseline splines).
A = basisSet.specs;
for kk = 1:basisSet.nMets
    A(:,kk) = conv(A(:,kk), lineShape, 'same');    
end

% Calculate the final metabolite residual

MetabFit = A(:,1:basisSet.nMets)  * ampl(1:basisSet.nMets);
% MetabFit = MetabFit *scale;


% Calculate the residual
% Cut out the frequency range of the spectrum to be fit
% Apply initial referencing shift
dataToFit   = op_freqshift(dataToFit, -refShift);
%dataToFit.fids = dataToFit.fids * exp(-1i*ph0);
dataToFit.specs = fftshift(fft(dataToFit.fids,[],1),1);
% dataToFit   = op_ampScale(dataToFit, scale);
%dataToFit   = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end),length(splineArray(:,1,1)));
%dataToFit.specs   = dataToFit.specs .* exp(-1i*ph1*multiplier);
%dataToFit.fids = ifft(fftshift(dataToFit.specs,1),[],1);
out.specs=dataToFit.specs-MetabFit;
out.fids = ifft(fftshift(out.specs,1),[],1);

out = op_truncate(out,2);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
