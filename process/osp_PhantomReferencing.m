function [refShift, refFWHM, refSinglets] = osp_PhantomReferencing(dataToFit)
%% [refShift, refFWHM] = osp_PhantomReferencing(dataToFit)
%   Identify appropriate reference signals then derive reference shift
%   accordingly, using the osp_XReferencing implementation.
%
%   Handles three possible phantom scenarios:
%
%   1) BRAINO-like phantom with prominent singlets around 2, 3 and 3.2 ppm
%      (amongst others), in which case default Cho/Cr referencing can be used
%   2) Non-brain-like phantoms with a reference singlet (DSS/TSP) at 0ppm (in
%      which case, use that peak), or
%   3) Non-brain-like phantoms without an identifiable reference peak (in which
%      case, do nothing)
%
%
%   Input:
%       dataToFit   = FID-A data structure
%
%   Output:
%       refShift    = Reference frequency shift (in Hz)
%       refFWHM     = full-width half-maximum of the input spectrum (in ppm)
%       refSinglets = singlets used in osp_XReferencing call
%
%   Author:
%       Alexander R Craven (Haukeland University Hospital, 2024-04-24)
%       alex.craven@uib.no
%
%   History:
%       2024-04-24: First version of the code.
%

possible_references=[0.00,2.01,3.03,3.22];

search_dist=0.5;

% use matlab's findpeaks function to identify major prominences in the spectrum {{{

if dataToFit.dims.averages>0
    dataToFit=op_averaging(dataToFit);
end

ppm=dataToFit.ppm;
spec=abs(dataToFit.specs);

% limit to the region around our expected peaks
mm=ppm<max(possible_references)+search_dist & ppm>min(possible_references)-search_dist;
ppm=ppm(mm);
spec=spec(mm);

spec=spec/max(spec); % normalise
if (ppm(end)<ppm(1)) % findpeaks is particular about order
    ppm=flip(ppm);
    spec=flip(spec);
end

% note that OspreyProcess applies 2Hz linebroadening before calling
% osp_PhantomReferencing; peak selection criteria (esp MinPeakWidth) may need
% to be revised if this function is to be used without such broadening
[pk,loc,w]=findpeaks(spec,ppm,'MinPeakProminence',0.05,'MinPeakWidth',0.02,'MaxPeakWidth',0.2,'MinPeakDistance',0.14,'Annotate','extents');

% }}}

% find distance to nearest match for our possible references {{{
d=repmat(possible_references',1,length(loc))-repmat(loc,length(possible_references),1);
nearest=min(abs(d),[],2);
[best_shift,best_ix]=min(nearest);
% }}}

if (std(nearest(2:4))<0.02)
    % seemingly stable estimate on the "big 3" candidates peaks, NAA, Cr, Cho
    % use default referencing (ref OspreyProcess)
    disp('Phantom: Found standard reference singlets (NAA, Cr, Cho); using default referencing');
    refSinglets=[3.03 3.22];
    [refShift, refFWHM] = osp_XReferencing(dataToFit,[3.03 3.22],[1 1],[1.85 4.2]);
    if abs(refShift) > 10 % This a huge shift. Most likley wrong and we will try it again with tNAA only
        refSinglets=[2.01];
        [refShift, refFWHM] = osp_XReferencing(raw,2.01,1,[1.85 4.2]);% determine frequency shift
    end
elseif (best_ix==1 && best_shift<0.08) || (nearest(1)<0.03)
    % found a close match for a 0ppm reference, use it
    disp('Phantom: Found DSS/TSP reference near 0ppm');
    refSinglets=[0.00];
    [refShift, refFWHM] = osp_XReferencing(dataToFit,[0.00],[1],[-1 1]);
else
    % Otherwise, do nothing.
    % Even if we have a near match to one of the other possible peaks, without
    % context we cannot safely assume that it actually corresponds to the
    % correct metabolite
    disp('Phantom: no reliable reference could be identified!')
    refSinglets=[];
    refShift=0;
    refFWHM=nan;
end
