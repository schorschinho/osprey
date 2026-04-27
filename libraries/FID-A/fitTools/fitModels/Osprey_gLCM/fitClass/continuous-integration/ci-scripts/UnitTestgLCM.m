function [results,rt] = UnitTestgLCM(test,parallel)
%% UnitTestgLCM.m
%   This function performs a quick unit test of different parts in the gLCM
%   algorithm
%
%   USAGE:
%       [results,rt] = UnitTestgLCM
%
%   OUTPUT:
%       results - Test results in MATLAB unittest class format.
%       rt - results in table format
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-07-24)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-06-04: First version of the code.


%% Run different scenarios through the gLCM 
warning('off','all')                            % Turn off warnings

if any(any(contains(struct2cell(ver), 'Parallel Computing Toolbox'))) && parallel
        results{1} = runtests([test '.m'],'UseParallel',true);
        rt{1} = table(results{1});
else
    results{1} = runtests([test '.m']);
    rt{1} = table(results{1});
end

warning('on','all')                             % Turn on warnings
end
