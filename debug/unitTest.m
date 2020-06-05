function [results,rt] = unitTest
%% unitTest.m
%   This function performs a unit test based on the jobTesting.m file in the
%   debug folder. All Osprey modules are tested in dependence of the
%   included dataset. Currently only command line calls of Opsrey are
%   tested. You can change the jobTesting.m file according to your desired
%   needs. The results unittest class includes detailed information about
%   the performance.
%
%   USAGE:
%       [results,rt] = unitTest
%
%   OUTPUT:
%       results - Test results in MATLAB unittest class format.
%       rt - results in table format    
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2019-11-07)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-06-03: First version of the code.

%% Call commandline tests
% You may want to change the input within the jobTesting file in the
% osprey/debug folder.

results = runtests('unitTestOspreyConsole.m');
rt = table(results)
end