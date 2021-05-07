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

% Create an array of normalized cubic baseline spline basis functions.
[splineArray]    = fit_makeSplineBasis(dataToFit, fitRangePPM, minKnotSpacingPPM);
if length(fitParams.beta_j)>size(splineArray,2)
[splineArray]    = fit_makeSplineBasis(dataToFit, fitRangePPM, 0.1);
end

% Apply phasing to the spline basis functions
B = [splineArray(:,:,1) + 1i*splineArray(:,:,2)];
B = B  * exp(1i*ph0);
B = B .* exp(1i*ph1*multiplier);
% B = [real(B)];

% Convolve the lineshape with the metabolite basis functions only
% (NOT the macromolecules or lipids or baseline splines).
A = basisSet.specs;
for kk = 1:nMets
    A(:,kk) = conv(A(:,kk), lineShape, 'same');
end

% Calculate the final baseline
baseline    = B * beta_j;
completeFit = A * ampl + baseline;


%Computation of the Jacobian
% The size is MxN with M (number of estimated parameters) and N (number of
% points).

%Store the indices of the different partial derivatives 
ph0Col = 1;
ph1Col = 2;
gaussLBCol = 3 * ones(1,nMets + nMM);
lorentzLBCol = gaussLBCol(end) + 1 : (gaussLBCol(1) + nMets + nMM);
freqShiftCol = lorentzLBCol(end) + 1 : (lorentzLBCol(end) + nMets + nMM);
AmplCol = freqShiftCol(end) + 1 : (freqShiftCol(end) + nMets + nMM);
SplineAmplCol = AmplCol(end) + 1 : (AmplCol(end) + size(splineArray,2));
lineshapeCol = SplineAmplCol(end)+1 : SplineAmplCol(end) + length(lineShape);

nparams = 3 + length(lorentzLBCol) + length(freqShiftCol) + length(AmplCol) + length(SplineAmplCol) + length(lineshapeCol);

J = zeros(npoints,nparams);
Jfd = zeros(npoints,nparams);


%derivative wrt ph0
% for ii = 1:nMets
%     col = ph0Col(ii);
%     J(:,col) = J(:,col)+j*basisSet.fids(:,ii)*ampl(ii);
%     Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);
%     Jfd(:,col) = conv(Jfd(:,col), lineShape, 'same') + B * beta_j;    
% end
Jfd(:,1) = j * completeFit;

%derivative wrt ph1
% for ii = 1:nBasisFcts
%     col = ph1Col(ii);
%     J(:,col) = J(:,col)+j*basisSet.fids(:,ii)*ampl(ii);
%     Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);
%     Jfd(:,col) = conv(Jfd(:,col), lineShape, 'same') + B * beta_j;    
% end
Jfd(:,2) = j *  completeFit .* multiplier;

%derivative wrt gaussLB
for ii = 1:nBasisFcts
    if ii <= nMets % Sum up derivarives of all metabolite functions first
        col = gaussLBCol(ii);
        J(:,col) = J(:,col)-basisSet.fids(:,ii).*(t.^2)';
        Jfd(:,col) = Jfd(:,col) +  conv(fftshift(fft(-basisSet.fids(:,ii).*(t.^2)',[],1),1), lineShape, 'same')*ampl(ii); 
    else  % No convolution is applied to the MM functions
        J(:,col) = J(:,col)-basisSet.fids(:,ii).*(t.^2)';
        Jfd(:,col) = Jfd(:,col) +  fftshift(fft(-basisSet.fids(:,ii).*(t.^2)',[],1),1)*ampl(ii);             
    end
end


%derivative wrt lorentzLB
for ii = 1:nBasisFcts
    if ii <= nMets % Sum up derivarives of all metabolite functions first
        col = lorentzLBCol(ii);
        J(:,col) = J(:,col)-basisSet.fids(:,ii).*t';
        Jfd(:,col) = Jfd(:,col) +  conv(fftshift(fft(J(:,col),[],1),1), lineShape, 'same')*ampl(ii);  
    else % No convolution is applied to the MM functions
        col = lorentzLBCol(ii);
        J(:,col) = J(:,col)-basisSet.fids(:,ii).*t';
        Jfd(:,col) = Jfd(:,col) +  fftshift(fft(J(:,col),[],1),1)*ampl(ii);         
    end
end

%derivative wrt freqShift
for ii=1:nBasisFcts
    if ii <= nMets  % Sum up derivarives of all metabolite functions first
        col = freqShiftCol(ii);
        J(:,col) = J(:,col)-j*basisSet.fids(:,ii).*t';
        Jfd(:,col) = Jfd(:,col) +  conv(fftshift(fft(J(:,col),[],1),1), lineShape, 'same')*ampl(ii);  
    else % No convolution is applied to the MM functions
        col = freqShiftCol(ii);
        J(:,col) = J(:,col)-j*basisSet.fids(:,ii).*t';
        Jfd(:,col) = Jfd(:,col) + fftshift(fft(J(:,col),[],1),1)*ampl(ii);          
    end
end

% derivative  wrt basis set  amplitudes 
for ii=1:nBasisFcts
    if ii <= nMets
        col = AmplCol(ii);
        J(:,col) = basisSet.fids(:,ii);
        Jfd(:,col) = conv(fftshift(fft(J(:,col),[],1),1), lineShape, 'same');
    else
        col = AmplCol(ii);
        J(:,col) = basisSet.fids(:,ii);
        Jfd(:,col) = fftshift(fft(J(:,col),[],1),1);    
    end
end


% derivative wrt spline amplitudes
for ii=1:length(SplineAmplCol)
    col = SplineAmplCol(ii);
    Jfd(:,col) = Jfd(:,col)+B(:,ii);
end


%derivative wrt lineshape 
% We will do a discrete convolution of S'*M by using the Toeplitz matrix
% form of the partial derivative of the lineshape vector S and M beeing 
% the metabolite basis functions. This is allowed as convolutions are 
% commutative and (M*S)' = M'*S = M*S' -> (M*S)' = S'*M.

%We need to create S' as a Toeplitz matrix. The derivatives will
%essantially be ones on the diagonal (first lineshape coeff) or the upper off
%diagonal (all other lineshape coeff). 

nLineShape = length(lineShape);
nPoints = length(basisSet.fids(:,1));

%Set up the first rows of each partial derivative
Toep1row = zeros(nLineShape,nPoints);
for ii = 1 : nLineShape
    Toep1row(ii,ii) = 1; 
end

%We will set it up as a 3D vector with the partial derivatives in the third dimensions and the
%Toeplitz matrix spanning the first two dimensions.
ToepLineShape = zeros(nPoints,nPoints,nLineShape);
for ii = 1 : nLineShape
    if ii <= 1
        ToepLineShape(:,:,ii) = toeplitz(Toep1row(ii,:));
    else
        ToepLineShape(:,:,ii) = tril(toeplitz(Toep1row(ii,:))); % extract lower triangle of the Toeplitz matrix
    end
end

for ii=1:nLineShape % Loop over the lineshape derivatives
    col = lineshapeCol(ii);
    for kk = 1 : nMets % Convolute all basis functions to create the full model function
        J(:,col) = J(:,col)+basisSet.fids(:,kk);
        Jfd(:,col) = Jfd(:,col) + ampl(ii)*(ToepLineShape(:,:,ii) * fftshift(fft(basisSet.fids(:,kk),[],1),1));
    end
        
end


% Select fit range
dataToFit   = op_freqrange(dataToFit,fitRangePPM(1),fitRangePPM(2));


% Calculate the residual
% Cut out the frequency range of the spectrum to be fit
% Apply initial referencing shift
dataToFit   = op_freqshift(dataToFit, -refShift);


% estimating the sigma based on the residual
sigma   = std(real(dataToFit.specs)-real(completeFit));

%calculate the fisher matrix
fisher = (1./(sigma^2)) * real(Jfd.'*Jfd);

% % Get non zero values from the fisher matrix
% non_zero_values = find(sum(fisher,1)~=0);
% non_zero_fisher = fisher(non_zero_values,non_zero_values);
% 
% %Calculate the invers fisher matrix
% non_zero_CRLB   = sqrt(inv(non_zero_fisher));
% 
% %Fill the non zero CRLBs back into the matrix
% CRLB = zeros(nBasisFcts);
% for ii = 1:length(non_zero_values)
%     for jj = 1:length(non_zero_values)
%         CRLB(non_zero_values(jj),non_zero_values(ii)) = non_zero_CRLB(jj,ii);
%     end
% end

 CRLB   = sqrt(pinv(fisher));
%Calculate relativ error of the amplitude estimates
for ii = 1:nBasisFcts
    col = AmplCol(ii);
   relative_CRLB(ii) = 100./(ampl(ii)/CRLB(col,col));
end

%Relative CRLBs in percent
relative_CRLB = real(relative_CRLB)';


end