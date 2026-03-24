function new = write_xxx(string,old)
% Dynamic changing of regular expressions with regexprep.
% Overwrites strings within a string with the same number of X.
%
% Author: Dr. Georg Oeltzschner
%         Johns Hopkins University, 08/23/2016
%
%   History:
%       2016-08-23: First version of the code.
%       2018-09-13: Updated to generalize string parts to be stripped.

% Strip unneeded " and white space.
old = regexprep(old, '\s*"\s*', '');

% Do actual overwriting.
new = strrep(string,old,repmat('X',[1 length(old)]));
end