function [MRSCont] = osp_remove_noise_scan(MRSCont)
%% MRSCont = osp_remove_noise_scan(MRSCont)
%   This function removes the noise reference from the raw data which is
%   added to the raw data acquired with the CMRR sequence (if the flag is
%   set to include a water and noise reference). The function should be
%   called after running OspreyLoad.
%
%   USAGE:
%       MRSCont = osp_remove_noise_scan(MRSCont)
%
%   INPUTS:
%      MRSCont     = Osprey MRS container
%
%   OUTPUTS:
%      MRSCont     = Osprey MRS container
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2022-07-05)
%       hzoelln2@jh.edu
%  
%%
for kk = 1 : MRSCont.nDatasets
    MRSCont.raw_ref{kk}.fids(:,2) = [];
    MRSCont.raw_ref{kk}.specs(:,2) = [];
    MRSCont.raw_ref{kk}.sz(2) = 1;
    MRSCont.raw_ref{kk}.averages = 1;
end

end