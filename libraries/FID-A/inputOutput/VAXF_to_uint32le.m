
function [ uint32le] = VAXF_to_uint32le(floatVAXF)
%VAXF_TO_UINT32LE Converts from VAXF (single precision) to IEEE-LE (UINT32)
% This function converts floating point numbers initialized in MATLAB  
% into equivalent raw 32bit unsigned integers (little endian)using the 
% specification for the VAXF floating point number format.
% The VAXF format is the single precision format: 
%  http://www.opengroup.org/onlinepubs/9629399/chap14.htm#tagfcjh_20
%
% The function is intended to be called from FWRITEVAX.
%
% See also VAXG_TO_UINT64LE, VAXD_TO_UINT64LE, FWRITEVAX
%
%  2009 The MathWorks, Inc. MATLAB and Simulink are registered trademarks
% of The MathWorks, Inc. See www.mathworks.com/trademarks for a list of 
% additional trademarks. Other product or brand names may be trademarks or 
% registered trademarks of their respective holders. 

%% Define floating value properties for VAX architecture
% The generic equation for a floating point number is:
% V = (-1)^double(S) * M * A^(double(E)-B);
% Substituting M = C + F 
% V = (-1)^double(S) * (C+F) * A^(double(E)-B);
%
% Performing inverse operations to solve for E and F:
% 0 &lt;= 1 + log(M)/log(2) &lt; 1  (VAX specific)
% E = (floor) (logV / D) + 1 + B   
% F = V / ((A ^ (E-B)) - C)
%
% V = value, S = sign, M = mantissa, A = base, E = exponent, B = exponent 
% bias, C = mantissa constant, F = fraction

A = 2;      % VAX specific
B = 128;    % VAX specific
C = 0.5;    % VAX specific
D = log(2); % VAX specific

%% Determine the sign bit. If -ve transform to positive.
S = zeros(size(floatVAXF));
if any(floatVAXF(:) < 0)
    indices = find(floatVAXF<0);
    floatVAXF(indices) = (-1) .* floatVAXF(indices);
    S = zeros(size(floatVAXF));
    S(indices) = 1;
end

%% Decompose the floating point number to SEF (Sign, Exp, Fract)
E = floor((log(floatVAXF)./ D) + 1 + B); 
F = ((floatVAXF ./ A.^(double(E)-B))) - C; 
% Convert floating point fraction to unsigned integer
F = floor(F * 16777216);   %VAX Specific 16777216=2^24

%% Shift the bits of S, E and F
S = bitshift(bitshift(uint32(S),0), 31);
E = bitshift(bitshift(uint32(E),0), 24);
F = bitshift(bitshift(uint32(F),0), 9);

%% Combine the S, E and F into the unsigned integer value
vaxInt = bitor(bitor(S,bitshift(E, -1)),bitshift(F,-9));

%% Swap WORD1 and WORD2
% VAX      &lt;-----WORD1-----&gt;&lt;-----WORD2-----&gt;
% IEEE-LE  &lt;-----WORD2-----&gt;&lt;-----WORD1-----&gt;

word1 = bitshift(vaxInt,16);
word2 = bitshift(vaxInt,-16);

uint32le = bitor(word1,word2);

%% Appendix 1
% 
% Determining C:
% C is the mantissa constant. VAX and IEEE normalization methods are as
% follows:
% VAX 	 (-1)s x 2e-b x .1f
% IEEE 	(-1)s x 2e-b x 1.f
% 
% By accounting for normalization, the significand or mantissa can be 
% written as:
% M = C+F 
% such that the values of C are: 0.5 (VAX) and 1 (IEEE). These constants
% simply imply the hidden bits:
% 
% References: 
% Floating point formats: http://www.quadibloc.com/comp/cp0201.htm
% Common Floating Point Representations: http://owen.sj.ca.us/~rk/howto/fltpt/index.html

%% Appendix 2
% 
% Range of the significand (M):
%
% For the maximum value of the fraction, all the fraction bits are 1.  
% For VAX, the fraction bits start from 1/4. The maximum value of the fraction equals the 
% sum of the geometric series: 1/4 + ... + 1/2^n, where n = 24
% General equation for the sum of geometric series is: a(1-r^n)/(1-r). For VAX a = 1/4, r = 1/2
% Max value: G = 1/4 + ... + 1/2^n = 1.0 - 1/2^n = 0.5 - 1/2^24
% Min value: 1/2
%
% Hence 1/2&lt;=M&lt;=1-1/2^24 or 1/2&lt;=M&lt;1.  

