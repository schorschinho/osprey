function [indMin, indMax] = ppmToIndex(ppmVector, ppmRange)
            
    % Returns indices of the ppm vector corresponding to a given
    % 2-element range
    [~,indMin] = min(abs(ppmVector - ppmRange(1)));
    [~,indMax] = min(abs(ppmVector - ppmRange(2)));
    
end