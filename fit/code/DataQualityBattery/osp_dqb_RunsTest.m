function[Runs_pVal,Stats] = osp_dqb_RunsTest(Residual, Range, Gap)
%% function[Runs_pVal,Stats] = osp_dqb_RunsTest(Residual, Range)
%
% Description: Uses the "runstest" to determine how random the residual is.
% The p-value  refers to the confidence in rejecting the null hypothesis 
% that the values in the data vector x come in random order, against the 
% alternative that they do not.
%
% Input:     Residual = Residual vector (data-model)
% Output:    Runs_pVal,  = p-value of the test
%
% Example usage:
%
% C.W. Davies-Jenkins, Johns Hopkins University 2025
arguments
Residual {mustBeVector} = [];
Range (1,2) {mustBeVector} = 1:length(Residual);
Gap = []
end

% If a Gap is supplied, exclude that range
if ~isempty(Gap)
    Incl = [Range(1):Gap(1), Gap(2):Range(2)];
else
    Incl = Range(1):Range(2);
end
Residual = Residual(Incl);

% "ud" returns a test decision based on the number of runs up or down. Too 
% few runs indicate a trend, while too many runs indicate an oscillation. 
% Values exactly equal to the preceding value are discarded.
[~,Runs_pVal,Stats] = runstest(real(Residual), 'ud');

end
