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
% OUTPUTS:
% splineArray = Matrix containing all cubic baseline spline functions.
%
% INPUTS:
% dataToFit   = FID-A data structure
% fitRangePPM = Vector containing the ppm range [PPMEND PPMST], 
%               e.g. [0.5 4.0]
% dkntmn      = Scalar denoting the minimal baseline knot spacing in ppm
%               (this is the hidden DKNTMN parameter in LCModel)
% plotFlag    = Boolean whether spline generation is supposed to be plotted
%               (Default = 0)

function [splineArray] = fit_makeSplineBasis(dataToFit, fitRangePPM, dkntmn, plotFlag)

% Parse input arguments
if nargin < 4
    plotFlag = 0;
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

% Previously, we used specific functions only available in the Curve
% Fitting Toolbox (spmak, fnval) to generate the splines. We can get around
% the need for the curve fitting toolbox if we generate them separately.
% % Loop over number of segments to generate all basis functions
% for rr = 1:length(beta_j)
%     beta_j(:) = 0;
%     beta_j(rr) = 1;
%     bFormSpline{rr} = spmak(knotLocations, beta_j);
%     % Evaluate spline
%     splineBasisFunctionMatrix(:,rr) = fnval(bFormSpline{rr}, baselineRangeExtended);
% end

% In the following, we recursively create the cubic (n = 4) B-spline basis
% functions. 
% See "Definition" under https://en.wikipedia.org/wiki/B-spline.
% This has been double-checked against the output from the MATLAB 'spmak'
% and 'fnval'functions that were previously used.

% First, for k = 1. We could do this probably more elegantly, but we only
% need to do it for four recursive iterations.
for ii = 1:length(knotLocations)-1
    for xx = 1:length(baselineRangeExtended) 
        if xx >= knotLocations(ii) && xx < knotLocations(ii+1)
            Bi1(xx, ii) = 1;
        else
            Bi1(xx, ii) = 0;
        end
    end
end

for ii = 1:length(knotLocations)-1
    for xx = 1:length(baselineRangeExtended)
        if knotLocations(ii) ~= knotLocations(ii+1)
            wi1(xx,ii) = (xx-knotLocations(ii))/(knotLocations(ii+1) - knotLocations(ii));
        else
            wi1(xx,ii) = 0;
        end
    end
end

% k = 2
for ii = 1:length(knotLocations)-2
    for xx = 1:length(baselineRangeExtended)
        Bi2(xx,ii) = wi1(xx,ii).*Bi1(xx,ii) + (1-wi1(xx,ii+1)).*Bi1(xx,ii+1);
    end
end

for ii = 1:length(knotLocations)-2
    for xx = 1:length(baselineRangeExtended)
        if knotLocations(ii) ~= knotLocations(ii+2)
            wi2(xx,ii) = (xx-knotLocations(ii))/(knotLocations(ii+2) - knotLocations(ii));
        else
            wi2(xx,ii) = 0;
        end
    end
end

% k = 3
for ii = 1:length(knotLocations)-3
    for xx = 1:length(baselineRangeExtended)
        Bi3(xx,ii) = wi2(xx,ii).*Bi2(xx,ii) + (1-wi2(xx,ii+1)).*Bi2(xx,ii+1);
    end
end

% n = 2
for ii = 1:length(knotLocations)-3
    for xx = 1:length(baselineRangeExtended)
        if knotLocations(ii) ~= knotLocations(ii+3)
            wi3(xx,ii) = (xx-knotLocations(ii))/(knotLocations(ii+3) - knotLocations(ii));
        else
            wi3(xx,ii) = 0;
        end
    end
end

% k = 4
for ii = 1:length(knotLocations)-4
    for xx = 1:length(baselineRangeExtended)
        Bi4(xx,ii) = wi3(xx,ii).*Bi3(xx,ii) + (1-wi3(xx,ii+1)).*Bi3(xx,ii+1);
    end
end

% Copy over into the array we previously used to store the splines in.
splineBasisFunctionMatrix = Bi4;

% visualize
if plotFlag
    figure
    hold;
    for rr = 1:length(beta_j)
        plot(baselineRangeExtended, rr + splineBasisFunctionMatrix(:,rr));
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
if plotFlag
    figure
    hold;
    for rr = 1:length(beta_j)
        plot(knotLocations(3):1:knotLocations(end-2), rr + splineInnerBasisFunctionMatrix(:,rr));
    end
    hold off;
end

% Return variables
splineArray = splineInnerBasisFunctionMatrix;

end
