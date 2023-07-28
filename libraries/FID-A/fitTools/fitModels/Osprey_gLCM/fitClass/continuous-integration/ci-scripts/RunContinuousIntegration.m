%% RunContinuousIntegration
%   This script triggers the continuous integration tests for the gLCM
%   model inclduing several test scenarios and model settings
%
%   USAGE:
%       results = RunContinuousIntegration;
%
%   INPUTS: None so far
%       
%
%   OUTPUTS:
%       results     = table with results.
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2023-07-24)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)    
%% Run synthetic data generator and create test data
parallel = 0;
[result_CreateSyntheticData,rt_CreateSyntheticData] = UnitTestgLCM('CreateSyntheticData',parallel);

% assert(sum(result_CreateSyntheticData{1}.Failed) == 0 && sum(result_CreateSyntheticData{1}.Incomplete) == 0)
%% Validate 1D baseline models
[result_BaselineModels1D,rt_BaselineModels1D] = UnitTestgLCM('BaselineModels1D',parallel);

%% Validate between step parsing of parameter results
[result_BetweenStepParser,rt_BetweenStepParser] = UnitTestgLCM('BetweenStepParser',parallel);

%% Validate default parameter parsing
[result_DefaultParameterParsing,rt_DefaultParameterParsing] = UnitTestgLCM('DefaultParameterParsing',parallel);

%% Validate external parameter parsing
[result_ExternalParameterParsing,rt_ExternalParameterParsing] = UnitTestgLCM('ExternalParameterParsing',parallel);

%% Validate 2D multi-transient modeling & parametrization of fixed/free parameters
[result_MultiTransient2D,rt_MultiTransient2D] = UnitTestgLCM('MultiTransient2D',parallel);

%% Validate 2D TE-series modeling

%% Valdiate CRLB results 


