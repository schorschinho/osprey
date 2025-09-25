function[Lip_2_tCr_Ratio] = osp_dqb_SignalToLipid(Data, PPM)
%% function[Lip_2_tCr_Ratio] = osp_dqb_SignalToLipid(Data, PPM)
%
% Description: Function that compares the integral of the signal in the
% 0.5–1.9 ppm region to that of tCr (2.93–3.13 ppm region).
%
% Input:     Data = Frequency-domain spectrum
%            PPM = ppm axis of "Data"
% Output:    Lip_2_tCr_Ratio = Ratio of the magnitude integrals of the
%               lipid region to the tCr region. 
%
% C.W. Davies-Jenkins, Johns Hopkins University 2025
arguments
Data {mustBeVector} = [];
PPM {mustBeVector} = [];
end

range_tcr = [2.93 3.13];
Data_tcr = Data(PPM>range_tcr(1) & PPM<range_tcr(2));

range_Lip = [0.5, 1.9];
Data_Lip = Data(PPM>range_Lip(1) & PPM<range_Lip(2));

Lip_2_tCr_Ratio = trapz(abs(Data_Lip)) ./ trapz(abs(Data_tcr));

end
