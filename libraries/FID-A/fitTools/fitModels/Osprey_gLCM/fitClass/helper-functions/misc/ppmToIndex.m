function [indMin, indMax] = ppmToIndex(ppmVector, ppmRange)
%%  [indMin, indMax] = ppmToIndex(ppmVector, ppmRange)
%   Returns indices of the ppm vector corresponding to a given 2-element
%   range
%
%   USAGE:
%       [indMin, indMax] = ppmToIndex(ppmVector, ppmRange)
%
%   INPUTS:
%       ppmVector      = step to plot    
%       ppmRange       = set plot range default is optimFreqFitRange
%
%   OUTPUTS:
%       indMin         = minimum index
%       indMax         = maximum index
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
%%  Get minimum and maximum index        

    [~,indMin] = min(abs(ppmVector - ppmRange(1)));        % Get minimum index
    [~,indMax] = min(abs(ppmVector - ppmRange(2)));        % Get maximum index
    
end