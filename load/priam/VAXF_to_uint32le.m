
function [ uint32le] = VAXF_to_uint32le(floatVAXF)
%VAXF_TO_UINT32LE Converts from VAXF (single precision) to IEEE-LE (UINT32)


A = 2;      % VAX specific
B = 128;    % VAX specific
C = 0.5;    % VAX specific
D = log(2); % VAX specific

% Determine the sign bit. If -ve transform to positive.
S = zeros(size(floatVAXF));
if any(floatVAXF(:) < 0)
    indices = find(floatVAXF<0);
    floatVAXF(indices) = (-1) .* floatVAXF(indices);
    S = zeros(size(floatVAXF));
    S(indices) = 1;
end

% Decompose the floating point number to SEF (Sign, Exp, Fract)
E = floor((log(floatVAXF)./ D) + 1 + B); 
F = ((floatVAXF ./ A.^(double(E)-B))) - C; 
% Convert floating point fraction to unsigned integer
F = floor(F * 16777216);   %VAX Specific 16777216=2^24

% Shift the bits of S, E and F
S = bitshift(bitshift(uint32(S),0), 31);
E = bitshift(bitshift(uint32(E),0), 24);
F = bitshift(bitshift(uint32(F),0), 9);

% Combine the S, E and F into the unsigned integer value
vaxInt = bitor(bitor(S,bitshift(E, -1)),bitshift(F,-9));

% Swap WORD1 and WORD2
% VAX      <-----WORD1-----><-----WORD2----->
% IEEE-LE  <-----WORD2-----><-----WORD1----->

word1 = bitshift(vaxInt,16);
word2 = bitshift(vaxInt,-16);

uint32le = bitor(word1,word2);

	
end