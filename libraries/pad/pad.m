function A = pad(B, newSize, paddedWith)
% Matlab function to pad an array to a desired new size
%
% By:	Christopher C. Wilcox, PhD
%		Naval Research Laboratory
%		Date: Mar 14, 2011
% 
% Usage:	A = pad(B, newSize, [paddedWith]);
%	Input:
%   - B is the input array
%   - newSize is a 2 element vector, [rows, cols], of the desired padded array
%   - paddedWith (Optional) can be used to pad with NaNs or Zeros
%       - 'nan' pads with NaNs (default)
%       - 'zero' pads with Zeros
% 
%   Output:
%   - A is the new array with the input array, B, surrounded by nans or
%   zeros
% 
% Modified by Helge Zoellner (Johns Hopkins University, 2020-05-11) hzoelln2@jhmi.edu  
% to add zeros in to 1-D FIDs.


if nargin > 1
    if nargin == 3
        if strcmp(paddedWith, 'zero')
            A = zeros(newSize,1);
        elseif strcmp(paddedWith, 'nan')
            A = nan(newSize,1);
        else
            error('Unrecognized token for padding value');
        end
    else
        A = nan(newSize);
    end
else
    error('Enter parameters for the new padded array');
end

[m, ~] = size(A);
[q, ~] = size(B);

if m < q
    error('The desired new array must be at least the size of the starting array');
end
    A(1:q, :) = B;
end