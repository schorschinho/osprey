%% unitTest.m
%   This function performs a command line unit test based on the jobTesting.m file in the
%   ci folder. All Osprey modules are tested in dependence of the
%   included dataset. You can change the jobTesting.m file according to your desired
%   needs. The results unittest class includes detailed information about
%   the performance.
%
%   USAGE:
%       tests = unitTestOspreyConsole
%
%   OUTPUT:
%       tests - Test results in MATLAB unittest class format.
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
%% 1. Call all modules with a functiontests array %%%
function tests = Osp_SinglePRESS_CL_w_Seg
% unitTestOspreyCommandLine
    tests = functiontests(localfunctions);
end
%% 2. Subfunctions for all Osprey Modules

%Test OspreyJob
function testOspreyJob(~)

    dir = strrep(which(['ci' filesep 'jobSinglePRESS.m']),'jobSinglePRESS.m','derivatives');
    if ~isempty(dir)
        delete(fullfile(dir,'jobSinglePRESS.mat'));
        delete(fullfile(dir,'LogFile.txt'));
    end
    MRSCont = OspreyJob(which(['ci' filesep 'jobSinglePRESS.m']),0,'11');
    addpath(dir);
end

%Test OspreyLoad
function testOspreyLoad(~)
    load(which(['ci' filesep 'derivatives' filesep 'jobSinglePRESS.mat']));
    MRSCont = OspreyLoad(MRSCont);
end

%Test OspreyProcess
function testOspreyProcess(~)
    load(which(['ci' filesep 'derivatives' filesep 'jobSinglePRESS.mat']));
    MRSCont = OspreyProcess(MRSCont);
end

%Test OspreyFit
function testOspreyFit(~)
    load(which(['ci' filesep 'derivatives' filesep 'jobSinglePRESS.mat']));
    MRSCont = OspreyFit(MRSCont);
end

%Test OspreyCoreg
function testOspreyCoreg(~)
    load(which(['ci' filesep 'derivatives' filesep 'jobSinglePRESS.mat']));
    if ~isempty(MRSCont.files_nii)
        MRSCont = OspreyCoreg(MRSCont);
    end
end

%Test OspreySegment
function testOspreySegment(~)
    load(which(['ci' filesep 'derivatives' filesep 'jobSinglePRESS.mat']));
    if ~isempty(MRSCont.files_nii)
        MRSCont = OspreySeg(MRSCont);
    end
end

%Test OspreyQuantify
function testOspreyQuantify(~)
    load(which(['ci' filesep 'derivatives' filesep 'jobSinglePRESS.mat']));
    MRSCont = OspreyQuantify(MRSCont);
end

%Test OspreyOverview
function testOspreyOverview(~)
    load(which(['ci' filesep 'derivatives' filesep 'jobSinglePRESS.mat']));
    MRSCont = OspreyOverview(MRSCont);
end
