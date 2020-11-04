% fit_makeSplineBasis.m
% Georg Oeltzschner, Johns Hopkins University 2020.
%
% USAGE:
% [splineArray] = fit_makeSplineBasis(dataToFit, fitRangePPM, dkntmn);
% 
% DESCRIPTION:
% This function creates an array of cubic baseline spline functions, which
% serve as a basis set for the LCModel baseline computation.
%
% The array contains the real part of the cubic splines in the (:,:,1)
% dimension, and the imaginary part in the (:,:,2) dimension.
%
% The matrix G is the regularizor matrix specified in Eq (3.13) in the
% original publication of the CONTIN algorithm behind LCModel:
% Provencher, Comp Phys Comm 27:213-227 (1982)
% and the LCModel paper itself:
% Provencher, Magn Reson Med 30:672-679 (1993).
% 
% OUTPUTS:
% splineArray = Matrix containing all cubic baseline spline functions.
% G           = Matrix containing the integral over the squared second
%               derivatives of the baseline functions
%
% INPUTS:
% dataToFit   = FID-A data structure
% fitRangePPM = Vector containing the ppm range [PPMEND PPMST], 
%               e.g. [0.5 4.0]
% dkntmn      = Scalar denoting the minimal baseline knot spacing in ppm
%               (this is the hidden DKNTMN parameter in LCModel)
% plot        = Boolean whether spline generation is supposed to be plotted
%               (Default = 0)

function [splineArray, G] = fit_makeSplineBasis(dataToFit, fitRangePPM, dkntmn, plot);

% Parse input arguments
if nargin < 4
    plot = 0;
end

% Cut out the frequency range of the spectrum to be fit
dataToFit       = op_freqrange(dataToFit, fitRangePPM(1), fitRangePPM(end));
% Create vector with the same number of points as the ppm axis
baselineRange   = linspace(1, length(dataToFit.ppm), length(dataToFit.ppm));

% Calculate number of total spline knots:
% First, calculate the designated knot spacing by rounding how many times
% the minimum knot spacing is going to fit into the fit range.
% Then, add 1 to arrive at the number of knots that will create these segments.
% Then, add 2 knots on either side of the fit range, as LCModel does.
% This nominal number of knots is the number of baseline parameters that 
% enter the LCModel.
nKnotsNominal = round((fitRangePPM(end) - fitRangePPM(1))/dkntmn) + 1 + 2;

% In order to create the splines that we ACTUALLY NEED (including the
% ones sitting on the edges), we need to add two additional (temporary) knots
% on either side of the fit range.
nKnotsActual = nKnotsNominal + 4;

% Generate an extended range so that we can create all the splines,
% including the one that lap over the edges of the actual fit range
nPtsBaseLineExtended = (nKnotsActual-1)*length(baselineRange)/(nKnotsNominal-1);
baselineRangeExtended = linspace(1, round(nPtsBaseLineExtended), round(nPtsBaseLineExtended));

% Generate a B-form cubic B-spline such that we have nSegments of cubic
% B-spline basis functions in the fit range
knotLocations = linspace(1, length(baselineRangeExtended), nKnotsActual);
nSegmentsActual = nKnotsActual - 4;
beta_j = zeros(1,nSegmentsActual);

% Loop over number of segments to generate all basis functions
for rr = 1:length(beta_j)
    beta_j(:) = 0;
    beta_j(rr) = 1;
    bFormSpline{rr} = spmak(knotLocations, beta_j);
    % Evaluate spline
    splineBasisFunctionMatrix(:,rr) = fnval(bFormSpline{rr}, baselineRangeExtended);
end
% visualize
if plot
    figure
    hold;
    for rr = 1:length(beta_j)
        plot(baselineRangeExtended, rr + splineBasisFunctionMatrix(:,rr));
        plot(knotLocations, rr + fnval(bFormSpline{rr}, knotLocations), 'x');
    end
    hold off;
end

% Get imaginary part through Hilbert transform
for rr = 1:size(splineBasisFunctionMatrix,2)
    splineBasisFunctionMatrix_Hilb(:,rr) = hilbert(splineBasisFunctionMatrix(:,rr));
end
splineBasisFunctionMatrix(:,:,1) = real(splineBasisFunctionMatrix_Hilb);
splineBasisFunctionMatrix(:,:,2) = -imag(splineBasisFunctionMatrix_Hilb);

warning off
% Now cut the functions out only between the nominal knots
for rr = 1:length(beta_j)
    splineInnerBasisFunctionMatrix(:,rr,:) = splineBasisFunctionMatrix(knotLocations(3):1:knotLocations(end-2),rr,:);
end
warning on
% visualize
if plot
    figure
    hold;
    for rr = 1:length(beta_j)
        plot(knotLocations(3):1:knotLocations(end-2), rr + splineInnerBasisFunctionMatrix(:,rr));
        plot(knotLocations, rr + fnval(bFormSpline{rr}, knotLocations), 'x');
    end
    hold off;
end

% Now generate the second derivatives of the real parts
for rr = 1:length(beta_j)
    bFormSpline1stDeriv{rr} = fnder(bFormSpline{rr});
    bFormSpline2ndDeriv{rr} = fnder(bFormSpline1stDeriv{rr});
    % Evaluate spline
    spline2ndDerivMatrix(:,rr) = fnval(bFormSpline2ndDeriv{rr}, knotLocations(3):1:knotLocations(end-2));
end
% visualize
if plot
    figure
    hold;
    for rr = 1:length(beta_j)
        plot(knotLocations(3):1:knotLocations(end-2), spline2ndDerivMatrix(:,rr));
        plot(knotLocations, fnval(bFormSpline2ndDeriv{rr}, knotLocations), 'x');
    end
    hold off;
end

% Now generate the G matrix
G = zeros(length(beta_j));
for ii = 1:length(beta_j)
    for jj = 1:length(beta_j)
        y = spline2ndDerivMatrix(:,ii) .* spline2ndDerivMatrix(:,jj);
        G(ii,jj) = trapz(knotLocations(3):1:knotLocations(end-2), y);
    end
end

% % Now test with an example
% beta_j = [2.083E-01   1.770E-01   1.455E-01   1.131E-01   8.091E-02   5.232E-02   2.836E-02   9.195E-03  -2.912E-03  -9.074E-03  -9.486E-03  -6.413E-03  -7.256E-04   7.255E-03   1.698E-02   2.800E-02   4.057E-02   5.314E-02];
% alphaB = 7.289E-01;
% regB = norm(alphaB*sqrtm(G)*sqrt((length(baselineRange)/length(knotLocations))^3)*beta_j')^2;
% calculate with actual number of points!!! this will lead you to the way
% that the number of points and knots is absorbed into the regularization
% parameter


% Return variables
splineArray = splineInnerBasisFunctionMatrix;

end
