function[MeanAbsCurvature] = osp_dqb_BaselineCurvature(Baseline, PPM, Range)
%% function[MeanAbsCurvature] = osp_dqb_BaselineCurvature(Baseline, PPM, Range)
%
% Description: Function that calculates the mean absolute curvature of the
% frequency-domain  (f) baseline (B):
%    MAC = mean(abs( d2B/df2  ./  (1 + dB/df)^2/3 ))
%
% Input:     Baseline = Baseline model
%            PPM = PPM axis
% Optional:  Range = indices defining model range. Curvature vector is
%                    truncated BEFORE calculation
%            
% Output:    MeanAbsCurvature = MAC, as described above
%
% C.W. Davies-Jenkins, Johns Hopkins University 2025
arguments
Baseline = []
PPM = []
Range = []
end

PPMstep = abs(PPM(1)-PPM(2));

% If a range is supplied, truncate the curvature vector before calculating
% the mean
if ~isempty(Range)
    Baseline = Baseline(Range(1):Range(2));
end

dBase  = gradient(Baseline,PPMstep);
d2Base = gradient(dBase,PPMstep);

Curvature = abs(d2Base) ./ (1+dBase.^2).^(3/2);

MeanAbsCurvature = mean(abs(Curvature),'omitnan');


end
