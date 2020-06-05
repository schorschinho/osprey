function [results,rt] = quickUnitTest
%% quickUnitTest.m
%   This function performs a quick unit test of the command line calls of
%   Osprey. It is based on two Philips datasets of the Big GABA paper, which
%   include a complete quantification of PRESS, MEGA-PRESS and HERMES
%   spectra. The results are stored in MATLAB unittest class format. The 
%   derivatives (see output folders) should also be checked.
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
%       2020-06-04: First version of the code.

%% 1. Store the jobTesting.m file
    jobTestStr = which('debug/jobTesting.m');
    backUpStr = strrep(jobTestStr,'jobTesting.m','jobTestingBackup.m');
    movefile(jobTestStr, backUpStr);

%% 2. PRESS dataset test
% Rename jobFile
    jobPRESSstr = which('debug/jobPRESSquickDebug.m');
    movefile(jobPRESSstr, jobTestStr);
    
% Run unit test   
    results{1} = runtests('unitTestOspreyConsole.m');
    rt{1} = table(results{1});
    rt{1}
    
% Store derivatives & rename
    derivativeStr = strrep(which('debug/jobTesting.m'),'jobTesting.m','derivatives');
    movefile(derivativeStr, strrep(derivativeStr,'derivatives','derivativesPRESS'));
    movefile(jobTestStr,jobPRESSstr);
%% 2. MEGA-PRESS dataset test
% Rename jobFile
    jobMEGAPRESSstr = which('debug/jobMEGAPRESSquickDebug.m');
    movefile(jobMEGAPRESSstr, jobTestStr);
    
% Run unit test      
    results{2} = runtests('unitTestOspreyConsole.m');
    rt{2} = table(results{2});
    rt{2}    

% Store derivatives & rename
    derivativeStr = strrep(which('debug/jobTesting.m'),'jobTesting.m','derivatives');
    movefile(derivativeStr, strrep(derivativeStr,'derivatives','derivativesMEGAPRESS'));
    movefile(jobTestStr,jobMEGAPRESSstr);
%% 3. HERMES dataset test
% This is currently not working for GSH and Lac edited HERMES data, future
% implemetation will have GSH and GABA HERMES dataset for the quick test.
% % Rename jobFile
%     jobHERMESstr = which('debug/jobHERMESquickDebug.m');
%     movefile(jobHERMESstr, jobTestStr);
%     
% % Run unit test      
%     results{3} = runtests('unitTestOspreyConsole.m');
%     rt{3} = table(results{3});
%     rt{3}    
% 
% % Store derivatives & rename  
%     derivativeStr = strrep(which('debug/jobTesting.m'),'jobTesting.m','derivatives');
%     movefile(derivativeStr, strrep(derivativeStr,'derivatives','derivativesHERMES'));
%     movefile(jobTestStr,jobHERMESstr);
%% 4. restore original jobTesting.m file
    movefile(backUpStr,jobTestStr);
end