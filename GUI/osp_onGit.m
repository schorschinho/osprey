function osp_onGit( ~, ~)
%% osp_onGit
%   Callback function on git button click. Visits GitHub Schorschinho.
%
%
%   USAGE:
%       osp_onGit( ~, ~ );
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

web('https://github.com/schorschinho/osprey', '-browser'); %Github Link
end 