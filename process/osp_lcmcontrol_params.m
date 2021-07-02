function LCMparam = osp_lcmcontrol_params(isMEGAdiff)
%%osp_LCMcontrol.m
% This function creates the basic struct for the creation of the LCModel control
% file. You need to fill in all the needed information below to make this
% work. For more information and parameter please refer to the LCModel
% manual. Remember to create the output folder on your Linux machine,
% otherwise the analysis will fail. You need to adapt all path variables.
%   OUTPUT:     Creates control file for LCmodel analysis
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-01-16)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2020-02-11: First version of the code.

%%% 1. DEFINE ARGUMENTS %%%

LCMparam.key = 210387309; %Your LCM key goes here
LCMparam.owner = 'Osprey processed spectra'; %Owner of the LCModel software package goes here
LCMparam.FILBAS = '/storage/myBasisSet_30ms_PRESS.BASIS'; % Location of .BASIS file used in LCModel
LCMparam.DKNTMN = 0.15; % Knot spacing parameter
LCMparam.DOWS = 'T'; % Do water referencing with the provided water reference. It will be turned off if no water reference is provided.
LCMparam.FOLDER = '/storage/LCMoutput'; %Output folder (you need to create this)
LCMparam.ATTH2O = 0.7; % No TE/TR correction. You can do this afterwards. Default is 0.7
LCMparam.ATTMET = 1; % No TE/TR correction.  You can do this afterwards. Default is 1
LCMparam.WCONC = 55510; % Assuming pure water concentration.
LCMparam.NEACH = 99; % Number of metabolites.
LCMparam.WDLINE = 0; % No ugly lines in the output.
LCMparam.CHCOMB={'PCh+GPC','Cr+PCr','NAA+NAAG','Glu+Gln','Glc+Tau'}; % Some combinations (You may need to adapt this).
LCMparam.NOMIT={'Gly','Ser'}; % Do not include on the flight simulations of Gly and Ser.
LCMparam.NAMREL = 'Cr+PCr'; % Calculate tCr ratios.
% LCMparam.CONREL =8; % Normalize this to 8 mM concentration.
LCMparam.DOECC= 'F'; % No eddy current correction in LCModel as this is already performed in Osprey.
if isMEGAdiff
    LCMparam.SPTYPE = 'mega-press-3'; %Set sequence type to MEGA-PRESS in the LCModel control file
end
end