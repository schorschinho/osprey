function LCMparam = osp_editControlParameters(LCMparam, key, value)
%% LCMparam = osp_readlcm_control(controlFile)
%   This function edits a set of control parameters saved as a struct.
%
%   If the assigned value is '', the key is being removed.
%
%   For now, this is relatively unsophisticated, i.e. there are no checks
%   whatsoever whether the key is a valid LCModel control parameter, or
%   whether the assigned value is valid.
%
%   USAGE:
%       LCMparam = osp_readlcm_control(controlFile)
%
%   INPUTS:
%      controlFile  = LCModel .control file
%      key          = string of LCModel control parameter to be changed
%      value        = string of the value to be assigned to 'key'
%
%   OUTPUTS:
%      LCMparam     = struct with LCModel control parameters
%
%   AUTHORS:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2021-08-05)
%       goeltzs1@jhmi.edu
%   

% Throw error if input is incomplete
if nargin < 3
    error('Need to provide three arguments: LCMparam, key, value.');
end

% If the key exists and the value is '', remove the field
if isfield(LCMparam, key)
    if strcmp(value, '')
        LCMparam = rmfield(LCMparam, key);
    else
        LCMparam.(key) = value;
    end
else
    if strcmp(value, '')
    else
        LCMparam.(key) = value;
    end

end        


end
