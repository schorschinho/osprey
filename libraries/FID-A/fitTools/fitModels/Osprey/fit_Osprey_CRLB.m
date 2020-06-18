function [J,Jfd,CRLB] = fit_Osprey_CRLB(dataToFit, basisSet, minKnotSpacingPPM, fitRangePPM,fitParams,refShift);

j = sqrt(-1); % i

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
beta_j      = fitParams.beta_j; % Amplitudes for baseline spline basis functions

% Normalize
lineShape = lineShape/sum(lineShape);

% Model function
% Y(f) = exp(i(ph0+ph1(f))) * (sum(betaj * B) + sum(amp conv(M, lineshape))
% with M = fft(basisSet.fids .* exp(-1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)')
% The derivatives are 'easier' in the time-domain, so we will stick with
% the FIDS for the calculations

% Create M for each metabolite function by applying the freqShifts,
% gaussian and lorentzian linebroadening, and the phase corrections to each
% basis function

t = basisSet.t;
for ii=1: nBasisFcts
    basisSet.fids(:,ii) = basisSet.fids(:,ii).* exp(-1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)';
     basisSet.fids(:,ii) = basisSet.fids(:,ii) * exp(1i*ph0);
end
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);
basisSet     = op_freqrange(basisSet,fitRangePPM(1),fitRangePPM(2));

% Create a ppm vector around a pivot point (water)
ppm_ax = basisSet.ppm;
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nBasisFcts
      basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp(1i*ph1*multiplier);
end
basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);
% basisSet.fids = basisSet.fids * ampl;

% dimensions and number of matrix entries
[npoints,~] = size(basisSet.fids); %number of points and number of basis functions
t = basisSet.t;

%Computation of the Jacobian
% The size is MxN with M (number of estimated parameters) and N (number of
% points). The parameters are (ampl, freqShift, lorentzLB, gaussLB, lineshape
% ,ph0 and ph1). In this case M = 7 * nMets + 6 * nMM.

J = zeros(npoints,7 * nMets + 6 * nMM);
Jfd = zeros(npoints,7 * nMets + 6 * nMM);

%Store the indices of the diferent partial derivatives 
amplCol = 1:nMets+nMM;
freqShiftCol = nMets+nMM+1 : 2*(nMets+nMM);
lorentzLBCol = 2*(nMets+nMM)+1 : 3*(nMets+nMM);
gaussLBCol = 3*(nMets+nMM)+1 : 4*(nMets+nMM);
lineshapeCol = 4*(nMets+nMM)+1 : (4*(nMets+nMM) + nMets);
ph0Col = ones(1,nMets+nMM)*lineshapeCol(end) + 1;
ph1Col = ones(1,nMets+nMM)*lineshapeCol(end) + 2;

% derivative  wrt amplitude 
for ii=1:nBasisFcts
    col = amplCol(ii);
    J(:,col) = J(:,col)+basisSet.fids(:,ii);
    Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);
    Jfd(:,col) = conv(Jfd(:,col), lineShape, 'same')*ampl(ii);
end

%derivative wrt freqShift
for ii=1:nBasisFcts
    col = freqShiftCol(ii);
    J(:,col) = J(:,col)+j*basisSet.fids(:,ii)*ampl(ii).*t';
    Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);
    Jfd(:,col) = conv(Jfd(:,col), lineShape, 'same');
end

%derivative wrt lorentzLB
for ii = 1:nBasisFcts
    col = lorentzLBCol(ii);
    J(:,col) = J(:,col)-basisSet.fids(:,ii)*ampl(ii).*t';
    Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);
    Jfd(:,col) = conv(Jfd(:,col), lineShape, 'same');    
end

%derivative wrt gaussLB
for ii = 1:nBasisFcts
    col = gaussLBCol(ii);
    J(:,col) = J(:,col)-basisSet.fids(:,ii)*ampl(ii).*(t.^2)';
    Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);
    Jfd(:,col) = conv(Jfd(:,col), lineShape, 'same');    
end

%derivative wrt lineshape S (not sure how to handle these)
for ii = 1:nBasisFcts
    col = lineshapeCol(ii);
    J(:,col) = J(:,col)+0;
    Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);
    Jfd(:,col) = conv(Jfd(:,col), lineShape, 'same');    
end

%derivative wrt ph0
for ii = 1:nBasisFcts
    col = ph0Col(ii);
    J(:,col) = J(:,col)+j*basisSet.fids(:,ii)*ampl(ii);
    Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);
    Jfd(:,col) = conv(Jfd(:,col), lineShape, 'same');    
end

%derivative wrt ph1
for ii = 1:nBasisFcts
    col = ph1Col(ii);
    J(:,col) = J(:,col)+j*basisSet.fids(:,ii)*ampl(ii);
    Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);
    Jfd(:,col) = conv(Jfd(:,col), lineShape, 'same');    
end


% Convolve the lineshape with the metabolite basis functions only
% (NOT the macromolecules or lipids or baseline splines).
A = basisSet.specs;
for kk = 1:nMets
    A(:,kk) = conv(A(:,kk), lineShape, 'same');
end


% Select fit range
dataToFit   = op_freqrange(dataToFit,fitRangePPM(1),fitRangePPM(2));
% Create an array of normalized cubic baseline spline basis functions.
[splineArray, ~]    = fit_makeSplineBasis(dataToFit, fitRangePPM, minKnotSpacingPPM);
if length(fitParams.beta_j)>size(splineArray,2)
[splineArray, ~]    = fit_makeSplineBasis(dataToFit, fitRangePPM, 0.1);
end

% Apply phasing to the spline basis functions
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
B = B  * exp(1i*ph0);
B = B .* exp(1i*ph1*multiplier);
% B = [real(B)];

% Calculate the final baseline
baseline    = B * beta_j;
completeFit = A * ampl + baseline;


% Calculate the residual
% Cut out the frequency range of the spectrum to be fit
% Apply initial referencing shift
dataToFit   = op_freqshift(dataToFit, -refShift);


% estimating the sigma based on the residual
sigma   = std(real(dataToFit.specs)-real(completeFit));

%calculate the fisher matrix
fisher = (1./(sigma^2)) * real(Jfd(:,1:nBasisFcts).'*Jfd(:,1:nBasisFcts));

% Get non zero values from the fisher matrix
non_zero_values = find(sum(fisher,1)~=0);
non_zero_fisher = fisher(non_zero_values,non_zero_values);

%Calculate the invers fisher matrix
non_zero_CRLB   = sqrt(inv(non_zero_fisher));

%Fill the non zero CRLBs back into the matrix
CRLB = zeros(nBasisFcts);
for ii = 1:length(non_zero_values)
    for jj = 1:length(non_zero_values)
        CRLB(non_zero_values(jj),non_zero_values(ii)) = non_zero_CRLB(jj,ii);
    end
end

%Calculate relativ error of the amplitude estimates
for ii = 1:nBasisFcts
   relative_CRLB(ii) = 100./(ampl(ii)/CRLB(ii,ii));
end

%Relative CRLBs in percent
relative_CRLB = real(relative_CRLB)';


end