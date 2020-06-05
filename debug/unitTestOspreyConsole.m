%% unitTest.m
%   This function performs a command line unit test based on the jobTesting.m file in the
%   debug folder. All Osprey modules are tested in dependence of the
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
function tests = unitTestOspreyConsole
% unitTestOspreyCommandLine
    tests = functiontests(localfunctions);
end
%% 2. Subfunctions for all Osprey Modules

%Test OspreyJob
function testOspreyJob(~)
    dir = strrep(which('debug/jobTesting.m'),'jobTesting.m','derivatives');
    if ~isempty(dir)
        delete(fullfile(dir,'jobTesting.mat'));
        delete(fullfile(dir,'LogFile.txt'));
    end
    MRSCont = OspreyJob(which('debug/jobTesting.m'),0,'11');
    addpath(dir);
end

%Test OspreyLoad
function testOspreyLoad(~)
    load(which('debug/derivatives/jobTesting.mat'));
    MRSCont = OspreyLoad(MRSCont);
end

%Test OspreyProcess
function testOspreyProcess(~)
    load(which('debug/derivatives/jobTesting.mat'));
    MRSCont = OspreyProcess(MRSCont);
end

%Test OspreyFit
function testOspreyFit(~)
    load(which('debug/derivatives/jobTesting.mat'));
    MRSCont = OspreyFit(MRSCont);
end

%Test OspreyCoreg
function testOspreyCoreg(~)
    load(which('debug/derivatives/jobTesting.mat'));
    if ~isempty(MRSCont.files_nii)
        MRSCont = OspreyCoreg(MRSCont);
    end
end

%Test OspreySegment
function testOspreySegment(~)
    load(which('debug/derivatives/jobTesting.mat'));
    if ~isempty(MRSCont.files_nii)
        MRSCont = OspreySeg(MRSCont);
    end
end

%Test OspreyQuantify
function testOspreyQuantify(~)
    load(which('debug/derivatives/jobTesting.mat'));
    MRSCont = OspreyQuantify(MRSCont);
end

%Test OspreyOverview
function testOspreyOverview(~)
    load(which('debug/derivatives/jobTesting.mat'));
    MRSCont = OspreyOverview(MRSCont);
end

%Reload test reloading the stored MRSContainer
function testOspreyReload(~)
    %Loading the finished MRSContainer as test
    MRSCont = OspreyJob(which('debug/jobTesting.m'),0,'01');
    
    %Test loading and plots
    MRSCont = OspreyLoad(MRSCont);
        osp_plotModule(MRSCont, 'OspreyLoad', 1, 'mets');
        if MRSCont.flags.hasRef
            osp_plotModule(MRSCont, 'OspreyLoad', 1, 'ref');
        end
        if MRSCont.flags.hasWater
            osp_plotModule(MRSCont, 'OspreyLoad', 1, 'w');
        end
        if MRSCont.flags.hasMM
            osp_plotModule(MRSCont, 'OspreyLoad', 1, 'mm');
        end
    
    %Test process and plots    
    MRSCont = OspreyProcess(MRSCont);
        Names = fieldnames(MRSCont.processed);
            for ss = 1 : length(Names)
                osp_plotModule(MRSCont, 'OspreyProcess', 1, Names{ss});
            end
    
    %Test fit and plots
    MRSCont = OspreyFit(MRSCont);
        if strcmp(MRSCont.opts.fit.style, 'Concatenated')
        temp = fieldnames(MRSCont.fit.results);
        if MRSCont.flags.isUnEdited
            Names = fieldnames(MRSCont.fit.results);
        end
        if MRSCont.flags.isMEGA
            Names = {'diff1','sum'};
            if length(temp) == 2
                Names{3} = temp{2};
            else if length(temp) == 3
                Names{3} = temp{2};
                Names{4} = temp{3};
                end
            end
        end
        if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
            Names = {'diff1','diff2','sum'};
            if length(temp) == 2
                Names{4} = temp{2};
            else if length(temp) == 3
                Names{4} = temp{2};
                Names{5} = temp{3};
                end
            end
        end
        else
            Names = fieldnames(MRSCont.fit.results);  
        end
        for ss = 1 : length(Names)
            osp_plotModule(MRSCont, 'OspreyFit', 1, Names{ss});
        end
    
    %Test coreg and plots
    if ~isempty(MRSCont.files_nii)
        MRSCont = OspreyCoreg(MRSCont);
            osp_plotModule(MRSCont, 'OspreyCoreg', 1);
    end
    
    %Test segment and plots   
    if ~isempty(MRSCont.files_nii)
        MRSCont = OspreySeg(MRSCont);
            osp_plotModule(MRSCont, 'OspreySeg', 1);
    end
    
    %Test quantify and plots
    MRSCont = OspreyQuantify(MRSCont);
    
    %Test overview and plots
    MRSCont = OspreyOverview(MRSCont);
    Names = fieldnames(MRSCont.processed);
    for ss = 1 : length(Names)
        osp_plotModule(MRSCont, 'OspreySpecOverview', 1, Names{ss});
        osp_plotModule(MRSCont, 'OspreyMeanOverview', 1, Names{ss});
    end
    
    if MRSCont.flags.isUnEdited
        osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, 'off-tCr', 'tNAA');
        osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, 'off-tCr', 'tCho');
        osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, 'off-tCr', 'Ins');
        osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, 'off-tCr', 'Glx');
    
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'off-tCr', 'tNAA', 'SNR');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'off-tCr', 'tCho', 'SNR');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'off-tCr', 'Ins', 'SNR');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'off-tCr', 'Glx', 'SNR');

        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'off-tCr', 'tNAA', 'FWHM');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'off-tCr', 'tCho', 'FWHM');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'off-tCr', 'Ins', 'FWHM');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'off-tCr', 'Glx', 'FWHM');
    end
    if MRSCont.flags.isMEGA
        if ~strcmp(MRSCont.opts.fit.style, 'Concatenated')
            osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, 'diff1-tCr', MRSCont.opts.editTarget{1});
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'diff1-tCr', MRSCont.opts.editTarget{1}, 'SNR');
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'diff1-tCr', MRSCont.opts.editTarget{1}, 'FWHM');
        else
            osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, 'conc-tCr', MRSCont.opts.editTarget{1});
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'conc-tCr', MRSCont.opts.editTarget{1}, 'SNR');
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'conc-tCr', MRSCont.opts.editTarget{1}, 'FWHM');    
        end
    end
    if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
        if ~strcmp(MRSCont.opts.fit.style, 'Concatenated')
            osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, 'diff1-tCr', MRSCont.opts.editTarget{1});
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'diff1-tCr', MRSCont.opts.editTarget{1}, 'SNR');
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'diff1-tCr', MRSCont.opts.editTarget{1}, 'FWHM');
            osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, 'diff1-tCr', MRSCont.opts.editTarget{2});
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'diff1-tCr', MRSCont.opts.editTarget{2}, 'SNR');
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'diff1-tCr', MRSCont.opts.editTarget{2}, 'FWHM');
        else
            osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, 'conc-tCr', MRSCont.opts.editTarget{1});
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'conc-tCr', MRSCont.opts.editTarget{1}, 'SNR');
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'conc-tCr', MRSCont.opts.editTarget{1}, 'FWHM');
            osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, 'conc-tCr', MRSCont.opts.editTarget{2});
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'conc-tCr', MRSCont.opts.editTarget{2}, 'SNR');
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, 'conc-tCr', MRSCont.opts.editTarget{2}, 'FWHM');
        end
    end
end