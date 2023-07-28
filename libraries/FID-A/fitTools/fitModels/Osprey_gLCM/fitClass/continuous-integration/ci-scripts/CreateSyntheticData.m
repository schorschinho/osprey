%% CreateSyntheticData.m
%   This function performs runs unit tests for the synthetic data generator
%   it is called within the unittest class in MATLAB see
%   RunContinuousIntegration.m
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
%       2020-06-03: First version of the code.
%% 1. Call all modules with a functiontests array %%%
function tests = CreateSyntheticData
% unitTestCreateSyntheticData
    tests = functiontests(localfunctions);
end
%% 2. Subfunctions for different data scenarios
% Please ammend this if you want to include new scenarios into the unit
% testing. The name of the function has to start with "test"

%Test File Availability
function testSyntheticDataFileAvailability(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/press/30/basis_philips_press30.mat')};
end

% Test Scyllo-Inositol peak
function test1DsITE30(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/press/30/basis_philips_press30.mat')};
    parameter.indirDim.flag = 0;                                            % averaged data
    parameter.indirDim.function = '';                                       % averaged data

    alter = [];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(26,1);
    zero.freqShift = zeros(26,1);
    zero.baseAmpl = 1;
    zero.lineShape = 1;


    dataSets = 1;                                                           % Create 1 dataset
    changedComb = 0;
    shareLorentzLB = 1;
    SignalSNR = [3.3,3.5];

    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 5.70 * ones(dataSets,1); %For gaussian LB take 5.7043 (mean of Big PRESS)
    overwrite.lorentzLB = 2.42 * ones(dataSets,26); % For lorentzian LB take 2.4238 (mean of sI Big PRESS)
    overwrite.freqShift = [];
    overwrite.metAmpl = [];
    overwrite.baseAmpl = [];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-06;

    overwrite.metAmpl = zeros(dataSets,26);
    overwrite.metAmpl(:,17) = 0.15; % Write sI amplitude

    zero.metAmpl = ones(34,1);
    indx = [18]; % Keep sI
    zero.metAmpl(indx) = 0;

    outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/1D/sINoBS');
    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite);
end

% Test Scyllo-Inositol peak transients
function test1DsITE30_SeparateTransients(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/press/30/basis_philips_press30.mat')};
    parameter.indirDim.flag = 1;
    parameter.indirDim.function = 'averages';
    parameter.indirDim.length = 16;

    alter = [];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(26,1);
    zero.freqShift = zeros(26,1);
    zero.baseAmpl = 1;
    zero.lineShape = 1;


    dataSets = 1;                                                           % Create 1 dataset
    changedComb = 0;
    shareLorentzLB = 1;
    SignalSNR = [3.3,3.5];

    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 5.70 * ones(dataSets,1); %For gaussian LB take 5.7043 (mean of Big PRESS)
    overwrite.lorentzLB = 2.42 * ones(dataSets,26); % For lorentzian LB take 2.4238 (mean of sI Big PRESS)
    overwrite.freqShift = [];
    overwrite.metAmpl = [];
    overwrite.baseAmpl = [];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-06;

    overwrite.metAmpl = zeros(dataSets,26);
    overwrite.metAmpl(:,17) = 0.15; % Write sI amplitude

    zero.metAmpl = ones(34,1);
    indx = [18]; % Keep sI
    zero.metAmpl(indx) = 0;

    outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/2D/sINoBSSeparateTransients');
    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite);
end



%Test Healthy 1D TE 30 ms no baseline no MM
function test1DNoBaselineNoMMHealthyTE30(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/press/30/basis_philips_press30.mat')};
    parameter.indirDim.flag = 0;                                            % averaged data
    parameter.indirDim.function = '';                                       % averaged data

    % Set parameters for healthy spectrum
    alter = [];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(26,1);
    zero.freqShift = ones(26,1);
    zero.baseAmpl = 1;
    zero.lineShape = 1;  

    dataSets = 1;                                                           % Create 1 dataset
    changedComb = 0;                                                        % Don't alter the combos
    shareLorentzLB = 1;                                                     % Same LB for all metabolites
    SignalSNR = [2.9,3.1];

    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 9.70 * ones(dataSets,1); 
    overwrite.lorentzLB = 3.40 * ones(dataSets,26);
    overwrite.freqShift = [];
    overwrite.metAmpl = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39, ...
                         0, 0, 0, 0, 0, 0, 0, 0];                           % MM amplitudes go here
    overwrite.baseAmpl = [];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-04;
    zero.metAmpl = zeros(34,1);

    outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/1D/NoBSNoMM');

    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite);
end

%Test Healthy 1D TE 30 ms no baseline with MMs
function test1DNoBaselineHealthyTE30(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/press/30/basis_philips_press30.mat')};
    parameter.indirDim.flag = 0;                                            % averaged data
    parameter.indirDim.function = '';                                       % averaged data

    % Set parameters for healthy spectrum
    alter = [];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(26,1);
    zero.freqShift = ones(26,1);
    zero.baseAmpl = 1;
    zero.lineShape = 1;  

    dataSets = 1;                                                           % Create 1 dataset
    changedComb = 0;                                                        % Don't alter the combos
    shareLorentzLB = 1;                                                     % Same LB for all metabolites
    SignalSNR = [2.9,3.1];

    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 9.70 * ones(dataSets,1); 
    overwrite.lorentzLB = 3.40 * ones(dataSets,26);
    overwrite.freqShift = [];
    overwrite.metAmpl = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39, ...
                         1.72, 0.12, 0.61, 0.44, 1.47, 2.16, 2.01, 0.10];
    overwrite.baseAmpl = [];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-04;
    zero.metAmpl = zeros(34,1);

    outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/1D/NoBS');

    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite);
end

%Test Healthy 1D TE 30 ms with Gaussian baseline and MMs
function test1DHealthyTE30(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/press/30/basis_philips_press30.mat')};
    parameter.indirDim.flag = 0;                                            % averaged data
    parameter.indirDim.function = '';                                       % averaged data

    % Set parameters for healthy spectrum
    alter = [];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(26,1);
    zero.freqShift = ones(26,1);
    zero.baseAmpl = 0;
    zero.lineShape = 1;  

    dataSets = 1;                                                           % Create 1 dataset
    changedComb = 0;                                                        % Don't alter the combos
    shareLorentzLB = 1;                                                     % Same LB for all metabolites
    SignalSNR = [2.9,3.1];

    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 9.70 * ones(dataSets,1); 
    overwrite.lorentzLB = 3.40 * ones(dataSets,26);
    overwrite.freqShift = [];
    overwrite.metAmpl = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39, ...
                         1.72, 0.12, 0.61, 0.44, 1.47, 2.16, 2.01, 0.10];
    overwrite.baseAmpl = [0.0161, 0.0024, 0.0136, 0.0873, 0.0884, 0.1412, 0.1335, 0.1338, 0.0144, 0.1055, 0.0932, 0.0463, 0.1648, 0.0923];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-04;
    zero.metAmpl = zeros(34,1);

    outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/1D/BSandMMs');

    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite,0,0,0);
end

%Test Healthy 1D TE 30 ms with exp MMs
function test1DHealthyExpMMTE30(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/press/30/basis_philips_press30.mat')};
    parameter.indirDim.flag = 0;                                            % averaged data
    parameter.indirDim.function = '';                                       % averaged data

    % Set parameters for healthy spectrum
    alter = [];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(19,1);
    zero.freqShift = ones(19,1);
    zero.baseAmpl = 1;
    zero.lineShape = 1;  

    dataSets = 1;                                                           % Create 1 dataset
    changedComb = 0;                                                        % Don't alter the combos
    shareLorentzLB = 1;                                                     % Same LB for all metabolites
    SignalSNR = [2.9,3.1];

    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 9.70 * ones(dataSets,1); 
    overwrite.lorentzLB = 3.40 * ones(dataSets,19);
    overwrite.freqShift = [];
    overwrite.metAmpl = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39, ...
                         7.2];
    overwrite.baseAmpl = [];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-04;
    zero.metAmpl = zeros(34,1);

    outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/1D/ExpMM');

    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite,0,1);
end

%Test Healthy 2D TE 30 ms no baseline no MM separate transients
function test2DNoBaselineNoMMHealthyTE30_SeparateTransients(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/press/30/basis_philips_press30.mat')};
    parameter.indirDim.flag = 1;
    parameter.indirDim.function = 'averages';
    parameter.indirDim.length = 16;

    % Set parameters for healthy spectrum
    alter = [];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(26,1);
    zero.freqShift = ones(26,1);
    zero.baseAmpl = 1;
    zero.lineShape = 1;  

    dataSets = 1;                                                           % Create 1 dataset
    changedComb = 0;                                                        % Don't alter the combos
    shareLorentzLB = 1;                                                     % Same LB for all metabolites
    SignalSNR = [2.9,3.1];

    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 9.70 * ones(dataSets,1); 
    overwrite.lorentzLB = 3.40 * ones(dataSets,26);
    overwrite.freqShift = [];
    overwrite.metAmpl = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39, ...
                         0, 0, 0, 0, 0, 0, 0, 0];                           % MM amplitudes go here
    overwrite.baseAmpl = [];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-04;
    zero.metAmpl = zeros(34,1);

    outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/2D/NoBSNoMMSeparateTransients');

    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite);
end

%Test Healthy 2D TE 30 ms no baseline no MM separate transients
function test2DNoBaselineNoMMHealthyTE30_NAA_Cr_Ins_SeparateTransients(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/press/30/basis_philips_press30.mat')};
    parameter.indirDim.flag = 1;
    parameter.indirDim.function = 'averages';
    parameter.indirDim.length = 16;

    % Set parameters for healthy spectrum
    alter = [];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(26,1);
    zero.freqShift = ones(26,1);
    zero.baseAmpl = 1;
    zero.lineShape = 1;  

    dataSets = 1;                                                           % Create 1 dataset
    changedComb = 0;                                                        % Don't alter the combos
    shareLorentzLB = 1;                                                     % Same LB for all metabolites
    SignalSNR = [2.9,3.1];

    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 9.70 * ones(dataSets,1); 
    overwrite.lorentzLB = 3.40 * ones(dataSets,26);
    overwrite.freqShift = [];
    overwrite.metAmpl = [0, 0, 3.13, 0, 0, 0, 0, 0, 0, 7.06, 0, 8.70, 0, 0, 0, 0, 0, 0, ...
                         0, 0, 0, 0, 0, 0, 0, 0];                           % MM amplitudes go here
    overwrite.baseAmpl = [];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-04;
    zero.metAmpl = zeros(34,1);

    outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/2D/NoBSNoMMNAACrInsSeparateTransients');

    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite);
end

%Test Healthy 2D TE 30 ms no baseline with MMs separate transients
function test2DNoBaselineHealthyTE30_SeparateTransients(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/press/30/basis_philips_press30.mat')};
    parameter.indirDim.flag = 1;
    parameter.indirDim.function = 'averages';
    parameter.indirDim.length = 16;

    % Set parameters for healthy spectrum
    alter = [];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(26,1);
    zero.freqShift = ones(26,1);
    zero.baseAmpl = 1;
    zero.lineShape = 1;  

    dataSets = 1;                                                           % Create 1 dataset
    changedComb = 0;                                                        % Don't alter the combos
    shareLorentzLB = 1;                                                     % Same LB for all metabolites
    SignalSNR = [2.9,3.1];

    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 9.70 * ones(dataSets,1); 
    overwrite.lorentzLB = 3.40 * ones(dataSets,26);
    overwrite.freqShift = [];
    overwrite.metAmpl = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39, ...
                         1.72, 0.12, 0.61, 0.44, 1.47, 2.16, 2.01, 0.10];
    overwrite.baseAmpl = [];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-04;
    zero.metAmpl = zeros(34,1);

    outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/2D/NoBSSeparateTransients');

    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite);
end

%Test Healthy 2D TE 30 ms with Gaussian baseline and MMs separate transients
function test2DHealthyTE30_SeparateTransients(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/press/30/basis_philips_press30.mat')};
    parameter.indirDim.flag = 1;
    parameter.indirDim.function = 'averages';
    parameter.indirDim.length = 16;
   
    % Set parameters for healthy spectrum
    alter = [];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(36,1);
    zero.freqShift = ones(36,1);
    zero.baseAmpl = 0;
    zero.lineShape = 1;  

    dataSets = 1;                                                           % Create 1 dataset
    changedComb = 0;                                                        % Don't alter the combos
    shareLorentzLB = 1;                                                     % Same LB for all metabolites
    SignalSNR = [2.9,3.1];
    
    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 9.70 * ones(dataSets,1); 
    overwrite.lorentzLB = 3.40 * ones(dataSets,26);
    overwrite.freqShift = [];
    overwrite.metAmpl = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39, ...
                         1.72, 0.12, 0.61, 0.44, 1.47, 2.16, 2.01, 0.10];
    overwrite.baseAmpl = [0.0161, 0.0024, 0.0136, 0.0873, 0.0884, 0.1412, 0.1335, 0.1338, 0.0144, 0.1055, 0.0932, 0.0463, 0.1648, 0.0923];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-04;
    zero.metAmpl = zeros(34,1);
    
    outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/2D/BSandMMsSeparateTransients');

    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite,0,0,0);
end

%Test Healthy 2D TE series no baseline and no MMs
function test2DNoBaselineNoMMHealthyTESeries(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/slaser/30/basis_philips_slaser30.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/50/basis_philips_slaser50.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/70/basis_philips_slaser70.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/90/basis_philips_slaser90.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/110/basis_philips_slaser110.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/130/basis_philips_slaser130.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/150/basis_philips_slaser150.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/170/basis_philips_slaser170.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/190/basis_philips_slaser190.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/210/basis_philips_slaser210.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/230/basis_philips_slaser230.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/250/basis_philips_slaser250.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/270/basis_philips_slaser270.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/290/basis_philips_slaser290.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/310/basis_philips_slaser310.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/330/basis_philips_slaser330.mat')};
    parameter.indirDim.flag = 1;
    parameter.indirDim.function = 'T2';
    parameter.indirDim.parameter = 'metAmpl';
    parameter.indirDim.length = 16;
    
    parameter.indirDim.expectation.mean = [105;90;144;105;75;222;82;99;122;229;99;263;90;221;144;86;107;102;26.45;35.45;16.4;16.9;14.25;26.45;25.9;14.25;];
    parameter.indirDim.expectation.SD = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;];
    parameter.indirDim.expectation.meanBS = [22];
    parameter.indirDim.expectation.SDBS = [0];
    
    alter =[];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(26,1);
    zero.freqShift = ones(26,1);
    zero.baseAmpl = 1;
    zero.lineShape = 1;
    
    
    dataSets = 1;
    changedComb = 0;
    shareLorentzLB = 1;
    SignalSNR = [2.9,3.1];
    
    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 9.70 * ones(dataSets,1); 
    overwrite.lorentzLB = 3.40 * ones(dataSets,26);
    overwrite.freqShift = [];
    overwrite.metAmpl = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39, ...
                         0, 0, 0, 0, 0, 0, 0, 0];
    overwrite.baseAmpl = [];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-04;
    zero.metAmpl = zeros(34,1);
    
     outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/2D/TE-SeriesNoBsNoMM');
    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite);
end

%Test Healthy 2D TE series no baseline
function test2DNoBaselineHealthyTESeries(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/slaser/30/basis_philips_slaser30.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/50/basis_philips_slaser50.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/70/basis_philips_slaser70.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/90/basis_philips_slaser90.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/110/basis_philips_slaser110.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/130/basis_philips_slaser130.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/150/basis_philips_slaser150.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/170/basis_philips_slaser170.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/190/basis_philips_slaser190.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/210/basis_philips_slaser210.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/230/basis_philips_slaser230.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/250/basis_philips_slaser250.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/270/basis_philips_slaser270.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/290/basis_philips_slaser290.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/310/basis_philips_slaser310.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/330/basis_philips_slaser330.mat')};
    parameter.indirDim.flag = 1;
    parameter.indirDim.function = 'T2';
    parameter.indirDim.parameter = 'metAmpl';
    parameter.indirDim.length = 16;
    
    parameter.indirDim.expectation.mean = [105;90;144;105;75;222;82;99;122;229;99;263;90;221;144;86;107;102;26.45;35.45;16.4;16.9;14.25;26.45;25.9;14.25;];
    parameter.indirDim.expectation.SD = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;];
    parameter.indirDim.expectation.meanBS = [22];
    parameter.indirDim.expectation.SDBS = [0];
    
    alter =[];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(26,1);
    zero.freqShift = ones(26,1);
    zero.baseAmpl = 1;
    zero.lineShape = 1;
    
    
    dataSets = 1;
    changedComb = 0;
    shareLorentzLB = 1;
    SignalSNR = [2.9,3.1];
    
    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 9.70 * ones(dataSets,1); 
    overwrite.lorentzLB = 3.40 * ones(dataSets,26);
    overwrite.freqShift = [];
    overwrite.metAmpl = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39, ...
                         1.72, 0.12, 0.61, 0.44, 1.47, 2.16, 2.01, 0.10];
    overwrite.baseAmpl = [];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-04;
    zero.metAmpl = zeros(34,1);
    
     outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/2D/TE-SeriesNoBS');
    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite);
end

%Test Healthy 2D TE series
function test2DHealthyTESeries(~)
    % Load parameter struct
    prior_knowledge_folder = fileparts(which(fullfile('utilities','simulation-prior-knowledge','MRS_BigPRESS_Philips_Struct.mat')));
    load(fullfile(prior_knowledge_folder,'MRS_BigPRESS_Philips_Struct.mat'));
    % Add basis set
    parameter.basisSet = {which('/fit/basissets/3T/philips/unedited/slaser/30/basis_philips_slaser30.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/50/basis_philips_slaser50.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/70/basis_philips_slaser70.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/90/basis_philips_slaser90.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/110/basis_philips_slaser110.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/130/basis_philips_slaser130.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/150/basis_philips_slaser150.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/170/basis_philips_slaser170.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/190/basis_philips_slaser190.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/210/basis_philips_slaser210.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/230/basis_philips_slaser230.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/250/basis_philips_slaser250.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/270/basis_philips_slaser270.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/290/basis_philips_slaser290.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/310/basis_philips_slaser310.mat'),...
                          which('/fit/basissets/3T/philips/unedited/slaser/330/basis_philips_slaser330.mat')};
    parameter.indirDim.flag = 1;
    parameter.indirDim.function = 'T2';
    parameter.indirDim.parameter = 'metAmpl';
    parameter.indirDim.length = 16;
    
    parameter.indirDim.expectation.mean = [105;90;144;105;75;222;82;99;122;229;99;263;90;221;144;86;107;102;26.45;35.45;16.4;16.9;14.25;26.45;25.9;14.25;];
    parameter.indirDim.expectation.SD = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;];
    parameter.indirDim.expectation.meanBS = [22];
    parameter.indirDim.expectation.SDBS = [0];
    
    alter =[];
    zero.ph0 = 1;
    zero.ph1 = 1;
    zero.gaussLB = 0;
    zero.lorentzLB = zeros(26,1);
    zero.freqShift = ones(26,1);
    zero.baseAmpl = 0;
    zero.lineShape = 1;
    
    
    dataSets = 1;
    changedComb = 0;
    shareLorentzLB = 1;
    SignalSNR = [2.9,3.1];
    
    overwrite.ph0 = [];
    overwrite.ph1 = [];
    overwrite.gaussLB = 9.70 * ones(dataSets,1); 
    overwrite.lorentzLB = 3.40 * ones(dataSets,26);
    overwrite.freqShift = [];
    overwrite.metAmpl = [0.19, 1.94, 3.13, 2.14, 1.46, 0.54, 1.07, 1.15, 7.23, 7.06, 1.28, 8.70, 1.35, 0.61, 2.63, 1.72, 0.13, 0.39, ...
                         1.72, 0.12, 0.61, 0.44, 1.47, 2.16, 2.01, 0.10];
    overwrite.baseAmpl = [0.0161, 0.0024, 0.0136, 0.0873, 0.0884, 0.1412, 0.1335, 0.1338, 0.0144, 0.1055, 0.0932, 0.0463, 0.1648, 0.0923];
    overwrite.lineShape = [];
    overwrite.noiseAmpl = 2.2788e-04;
    zero.metAmpl = zeros(34,1);
    
     outputFolder = strrep(which('libraries/FID-A/fitTools/fitModels/Osprey_gLCM/fitClass/continuous-integration/ci-scripts/RunContinuousIntegration.m'),'ci-scripts/RunContinuousIntegration.m', 'data/2D/TE-SeriesBSandMM');
    [MRSCont] = OspreyGenerateSpectra(dataSets,outputFolder,1,parameter,alter,changedComb,zero,shareLorentzLB,SignalSNR,overwrite,0,0,1);
end
