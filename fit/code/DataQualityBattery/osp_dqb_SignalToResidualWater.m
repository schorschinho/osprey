function[Wat_2_tCr_Ratio] = osp_dqb_SignalToResidualWater(Data, PPM)
%% function[Wat_2_tCr_Ratio] = osp_dqb_SignalToResidualWater(Data, PPM)
%
% Description: Function that compares the integral of the signal in the
% 4.5–4.85 ppm region to that of tCr (2.93–3.13 ppm region).
%
% Input:     Data = Frequency-domain spectrum
%            PPM = ppm axis of "Data"
% Output:    Wat_2_tCr_Ratio = Ratio of the magnitude integrals of the
%               residual water region to the tCr region.
%
% C.W. Davies-Jenkins, Johns Hopkins University 2025
arguments
Data {mustBeVector} = [];
PPM {mustBeVector} = [];
end

range_tcr = [2.93 3.13];
Data_tcr = Data(PPM>range_tcr(1) & PPM<range_tcr(2));

range_Wat = [4.5, 4.85];
Data_Wat = Data(PPM>range_Wat(1) & PPM<range_Wat(2));

Wat_2_tCr_Ratio = trapz(abs(Data_Wat)) ./ trapz(abs(Data_tcr));

end
