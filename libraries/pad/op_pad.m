function A = op_pad(B, newSize, paddedWith)
% Matlab function to pad an array to a desired new size
%
% By:	Christopher C. Wilcox, PhD
%		Naval Research Laboratory
%		Date: Mar 14, 2011
% 
% Usage:	A = op_pad(B, newSize, [paddedWith]);
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
        if ndims(B) == 2 % MRS data
            if strcmp(paddedWith, 'zero')
                A = zeros(newSize,size(B,2));
            elseif strcmp(paddedWith, 'nan')
                A = nan(newSize,size(B,2));
            else
                error('Unrecognized token for padding value');
            end
        else
            mrsi_matrix_sz = size(B);
            mrsi_matrix_sz(1) = newSize;
            if strcmp(paddedWith, 'zero')
                A = zeros(mrsi_matrix_sz);
            elseif strcmp(paddedWith, 'nan')
                A = nan(mrsi_matrix_sz);
            else
                error('Unrecognized token for padding value');
            end
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
    if ndims(B) == 2 % MRS data
        A(1:q, :) = B;
    else if ndims(B) == 3 % 2D MRSI data
        A(1:q, :,:) = B;
    else %3D or 2D multi-slice MRSI data
        A(1:q, :,:,:) = B;
    end
    end
end