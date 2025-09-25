function[anyNegative, belowBaseline] = osp_dqb_Verticalshifts(Data,Baseline,Range)
%% function[anyNegative, belowBaseline] = osp_dqb_Verticalshifts(Data,Baseline,Range)
%
% Description: Computes the "vertical shifts" DQ metric proposed by B. 
%     Beroukhim et al (http://doi.org/10.1111/jon.13246). It checks
%     whether the spectrum/model drops below zero/baseline.
%
% Input:     Data = Frequency-domain spectrum (or model) vector
% Optional:  Baseline = Baseline vector from the linear combination model
%            Range = 2-element vector indicating range (indices)
% Output:    AnyNegative = Did data drop below zero? [bool]
%            belowBaseline = Did data drop below baseline? [bool]
%
% Example usage:
%
% C.W. Davies-Jenkins, Johns Hopkins University 2025
arguments
Data {mustBeVector} = []
Baseline {mustBeVector} = []
Range (1,2) {mustBeVector} = []
end

% If a range is supplied, truncate the data (and baseline)
if exist("Range","var") && ~isempty(Range)
    Data = Data(Range(1):Range(2));
    if exist("Baseline","var") && ~isempty(Baseline)
        Baseline = Baseline(Range(1):Range(2));
    end
end

anyNegative = any(Data<0);

if exist("Baseline","var") && ~isempty(Baseline)
    belowBaseline = any(Data<Baseline);
else
    belowBaseline = [];
end

end