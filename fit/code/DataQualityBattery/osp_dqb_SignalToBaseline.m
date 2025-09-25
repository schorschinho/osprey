function[Baseline_2_tCr_Ratio] = osp_dqb_SignalToBaseline(Data, PPM, Baseline, Range)
%% function[Baseline_2_tCr_Ratio] = osp_dqb_SignalToBaseline(Data, PPM, Baseline, Range)
%
% Description: Function that compares the integral of the signal in the
% 4.5–4.85 ppm region to that of tCr (2.93–3.13 ppm region).
%
% Input:     Data = Frequency-domain spectrum (or model)
%            PPM = ppm axis of "Data"
%            Baseline = Modeled baseline
% Optional:  Range = indices of the model range           
% Output:    Wat_2_tCr_Ratio = Ratio of the magnitude integrals of the
%               residual water region to the tCr region.
%
% C.W. Davies-Jenkins, Johns Hopkins University 2025
arguments
Data {mustBeVector} = [];
PPM {mustBeVector} = [];
Baseline = []
Range = []
end

% If a range is supplied, truncate the resiudal (and baseline)
if ~isempty(Range)
    Data = Data(Range(1):Range(2));
    PPM = PPM(Range(1):Range(2));
    Baseline = Baseline(Range(1):Range(2));
end

% Subtract the baseline from the tCr region
range_tcr = [2.93 3.13];
Data_tcr = Data(PPM>range_tcr(1) & PPM<range_tcr(2)) - Baseline(PPM>range_tcr(1) & PPM<range_tcr(2));

% Ratio of the magnitude of the baseline signal to tCr (minus the baseline)
Baseline_2_tCr_Ratio = trapz(abs(Baseline)) ./ trapz(abs(Data_tcr));

end
