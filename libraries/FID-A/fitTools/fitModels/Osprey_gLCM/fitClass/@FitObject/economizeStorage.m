function economizeStorage(obj, removeBasis, removeJacobian)
%%  economizeStorage(obj, removeBasis, removeJacobian)
%   This method removes redundant and optional entries from a OspreyFitObj.
%
%   USAGE:
%       obj.economizeStorage(removeBasis, removeJacobian)
%
%   INPUTS:
%       removeBasis     = flag to remove the basis set.  % OPTIONS:    - '1' (default)
%                                                                      - '0' (no)
%       removeJacobian  = flag to remove the basis set.  % OPTIONS:    - '1' (default)
%                                                                      - '0' (no)
%
%   OUTPUTS:
%       obj     = OspreyFitObj.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-03-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%%  Diverge to default options if required 
    if nargin < 3
            removeJacobian = 1;                 % set removeJacobian flag
        if nargin < 2
            removeBasis = 1;                    % set removeBasis flag
        end
    end    
    steps = obj.step;                           % get number of steps
    
%% This part removes the basis set

    if removeBasis
        obj.BaselineBasis   = [];               % remove spline basis
        obj.BasisSets.fids  = [];               % remove basis set
    end
    
%% This part removes the final jacobian from all model steps

    if removeJacobian && isfield(obj.Model{1,1}.info,'jacobian')
        for ss = 1 : steps                      % loop over all model steps
            obj.Model{1,ss}.info.jacobian=[];   % remove jacobian
        end
    end
end