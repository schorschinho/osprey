function osp_onPub( ~, ~)
%% osp_onPub
%   Callback function on publication button. Visits HERCULES paper.
%
%
%   USAGE:
%       osp_onPub( ~, ~ );
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


web('https://www.ncbi.nlm.nih.gov/sites/myncbi/1lo2w1io3qTgqH/collections/59221405/public/', '-browser'); %HERCULES paper will create a pubmed collection for referencing all papers

end 