function PolyArray = osp_gLCM_makePolyBasis(dataToFit, fitRangePPM, order, plotFlag)
%%  PolyArray = osp_gLCM_makePolyBasis(dataToFit, fitRangePPM, order, plotFlag)
%   This function generates a polynomial baseline basis set
%
%   USAGE:
%       PolyArray = osp_gLCM_makePolyBasis(dataToFit, fitRangePPM, order, plotFlag)
%
%   INPUTS:
%       dataToFit     = FID-A MRS struct (2D is modeled along extra dim)
%       fitRangePPM   = model range
%       order         = order of the polinom
%       plotFlag      = plot flag
%
%   OUTPUTS:
%       PolyArray     = polinomial baseline basis
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%%  Diverge to default options if required 

    if nargin < 4
        plotFlag = 0;                               % Don't plot
    end

%% Create basis
   
    [indMin, indMax] = ppmToIndex(dataToFit.ppm, fitRangePPM);   % Get indices for ppm range
    ppm = dataToFit.ppm(indMin:indMax);                          % Cut out the frequency range of the spectrum to be fit                       
    
    baselineRange   = linspace(1, length(ppm), length(ppm));    % Create vector with the same number of points as the ppm axis    
    PolyBasisFunctionMatrix = zeros(length(baselineRange),order+1); % Create matrix with the same number of points as teh ppm axis and order + 1
    
    for rr = 1:size(PolyBasisFunctionMatrix,2)                  % Loop over polynomial order
        p = zeros(rr,1);
        p(1) = 1;
        PolyBasisFunctionMatrix(:,rr) = polyval(p,ppm);         % Get polynomial value
    end                                                         % End loop over polynomial order
    
    
    for rr = 1:size(PolyBasisFunctionMatrix,2)                  % Get imaginary part through Hilbert transform
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