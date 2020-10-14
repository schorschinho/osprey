function count = fwriteVAXG(fid, A, precision)
% FWRITEVAXG(FID, A , PRECISION) writes the elements of 'A' to the 
% specified file in VAXG format, translating MATLAB ® values to the specified
% precision. COUNT is the number of elements successfully written.  
% 
% FID is an integer file identifier obtained from FOPEN. FWRITEVAXG requires
% FOPEN to open the output file in IEEE little-endian machineformat.
%
% PRECISION controls the form and size of result. See FREAD for supported
% formats.
%
% Usage:
%   A = rand(3,3); 
% 	fid = fopen('myFile', 'w', 'ieee-le');
%	count = fwriteVAXG(fid, A, 'double');
%   fclose(fid);
%
% The function is intended to be called by the user.
%
% ® 2009 The MathWorks, Inc. MATLAB and Simulink are registered trademarks
% of The MathWorks, Inc. See www.mathworks.com/trademarks for a list of 
% additional trademarks. Other product or brand names may be trademarks or 
% registered trademarks of their respective holders.

% Check for proper number of input arguments
if nargin < 2
    error('Not enough input arguments.')
end

% Check that the file has been open in ieee-le machineformat
[filename, permission, machineformat] = fopen(fid);
if ~strcmp(machineformat, 'ieee-le')
    error('Use FOPEN with ieee-le precision');
end

switch precision

    case {'float32', 'single'}
      rawUINT32 = VAXF_to_uint32le(A);
      count = fwrite(fid, rawUINT32, 'uint32');

    case {'float64', 'double'}
      rawUINT32 = VAXG_to_uint64le(A);
      count = fwrite(fid, rawUINT32, 'uint32');
      count = count/2;%2 UINT32 pieces for each double precision number

    case {'float'}
      if intmax == 2147483647 %32bit OS float is 32 bits
         rawUINT32 = VAXF_to_uint32le(A);
         count = fwrite(fid, rawUINT32, 'uint32');
      else
         rawUINT32 = VAXG_to_uint64le(A);
         count = fwrite(fid, rawUINT32, 'uint32');
         count = count/2;%2 UINT32 pieces for each double precision number
      end

    otherwise

      count = fwrite(fid, A, precision, 'vaxg');

end

end