function osp_onOsp( ~, ~)
%% osp_onOsp
%   Callback function on osprey button click. Opens mail contact.
%
%
%   USAGE:
%       osp_onOsp( ~, ~ );
%
%
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-01-16)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-01-16: First version of the code.
%%% 1. Visit websites %%%


email = 'gabamrs@gmail.com'; %Mail
url = ['mailto:',email];
web(url)
end 