%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright:
% Jun Tan
% University of Texas Southwestern Medical Center
% Department of Radiation Oncology
% Last edited: 06/16/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TF = is_handle(h)
% Test for valid handle.
% Difference from MATLAB's built-in ishandle, this function:
% 1) returns false when argument is []. 
% 2) only tests one input value, not array.

% TF = ~isempty(h) && 1 == numel(h) && isa(h, 'double') && ishandle(h);
TF = ~isempty(h) && 1 == numel(h) && ishandle(h);

% 06/16/2015
% Many thanks to Matthew Lewis at UT Southwestern, Radiology.
% He found that function "is_handle" threw an assertion failure in Matlab 2015a.
% By trial and error, he was able to isolate the failed assertion from
% the 4 assertions in this function.
% In older Matlabs, a handle is also of type "double", which does not seem
% to be true since 2015a.
