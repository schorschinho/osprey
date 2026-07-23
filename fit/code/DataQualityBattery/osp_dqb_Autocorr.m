function[SSAC, MaxAC, AC] = osp_dqb_Autocorr(Residual, Range, Gap)
%% function[SSAC, MaxAC, AC] = osp_dqb_Autocorr(Data, Range, Gap)
%
% Description: Calculates the autocorrelation of the supplied data. Useful
% for looking at the degree of structure in the residuals. This uses the
% Matlab xcorr function with max lag length of 75. This function will
% produce a warning if the data vector length<150, filling the ouput 
% vectors with NaNs.
%
% Input:     Residual = Vector to run autocorrelation on. Usually residual.
% Optional:  Range = A pair of indices to truncate the Data vector
%            Gap = A pair of indices to define excludes
% Output:    SSAC = Sum of squared AC values (Lag>0)
%            MaxAC = The absolute max AC (Lag>0)
%            AC = Normalized autocorrelation vector for visualization
%
% C.W. Davies-Jenkins, Johns Hopkins University 2025
arguments
Residual = []
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

Residual = Residual - mean(Residual); % Removes the effect of DC offset in vector.

N = length(Residual);

% If data vector
if N<150
    warning('Data vector too small to determine autocorrelation!')
    SSAC=nan;
    MaxAC=nan;
    AC = nan(size(Residual));
    return
end

% Compute autocorr. Number of lags restricted 50 to avoid noise wobbles
% 'unbiased' option normalizes autocorr to avoid linear trends as a
% function of lag.
[AC,lags] = xcorr(real(Residual), 75, 'unbiased'); 

AC = AC / AC(lags==0); % Normalization

MaxAC = max(abs(AC(~(lags==0))));    % Max AC value (detects large spikes)
SSAC = sum((AC(~(lags==0))).^2);      % Sum of squared AC

end