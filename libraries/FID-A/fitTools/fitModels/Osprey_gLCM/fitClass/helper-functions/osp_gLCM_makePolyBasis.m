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
% order       = order of the polynomial baseline
% plotFlag    = Boolean whether spline generation is supposed to be plotted
%               (Default = 0)

function [PolyArray] = osp_gLCM_makePolyBasis(dataToFit, fitRangePPM, order, plotFlag)

% Parse input arguments
if nargin < 4
    plotFlag = 0;
end

% Cut out the frequency range of the spectrum to be fit
[indMin, indMax] = ppmToIndex(dataToFit.ppm, fitRangePPM);
ppm = dataToFit.ppm(indMin:indMax);
% Create vector with the same number of points as the ppm axis
baselineRange   = linspace(1, length(ppm), length(ppm));

PolyBasisFunctionMatrix = zeros(length(baselineRange),order+1);

for rr = 1:size(PolyBasisFunctionMatrix,2)
    p = zeros(rr,1);
    p(1) = 1;
    PolyBasisFunctionMatrix(:,rr) = polyval(p,ppm);
end

% Get imaginary part through Hilbert transform
for rr = 1:size(PolyBasisFunctionMatrix,2)
    PolyBasisFunctionMatrix_Hilb(:,rr) = hilbert(PolyBasisFunctionMatrix(:,rr));
end
PolyBasisFunctionMatrix(:,:,1) = real(PolyBasisFunctionMatrix_Hilb);
PolyBasisFunctionMatrix(:,:,2) = -imag(PolyBasisFunctionMatrix_Hilb);

PolyBasisFunctionMatrix = PolyBasisFunctionMatrix(:,:,1) + 1j.*PolyBasisFunctionMatrix(:,:,2);

% Now cut the functions out only over the length of the frequency axis
eitherSide = round((length(baselineRange) - length(dataToFit.ppm))/2);

% We append with zeros and do a circshift so we can comfortably move the
% spline array along
appendWith = zeros(1000, size(PolyBasisFunctionMatrix,2));
PolyBasisFunctionMatrix = cat(1, appendWith, PolyBasisFunctionMatrix, appendWith);

% find the required shift
dppm = dataToFit.ppm(1) - dataToFit.ppm(2);
neededShift = (mean(dataToFit.ppm) - mean(ppm))/dppm;
PolyBasisFunctionMatrix = circshift(PolyBasisFunctionMatrix, floor(neededShift), 1);

% crop
splineInnerBasisFunctionMatrix = PolyBasisFunctionMatrix(1000+eitherSide+1:1000+eitherSide+length(dataToFit.ppm),:);

% Return variables
PolyArray = splineInnerBasisFunctionMatrix;

end

        function [indMin, indMax] = ppmToIndex(ppmVector, ppmRange)
            
            % Returns indices of the ppm vector corresponding to a given
            % 2-element range
            [~,indMin] = min(abs(ppmVector - ppmRange(1)));
            [~,indMax] = min(abs(ppmVector - ppmRange(2)));
            
        end