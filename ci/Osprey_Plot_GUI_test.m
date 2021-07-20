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
function tests = Osprey_Plot_GUI_test
% unitTestOspreyCommandLine
    tests = functiontests(localfunctions);
end
%% 2. Subfunctions to call GUI and plots


%Test OspreyGUI
function testOspreyGUI(~)
    folder_path = strrep(which(['ci' filesep 'UnitTest.m']),'UnitTest.m','derivatives');
    filestruct =  dir([folder_path filesep '*.mat']);
    load([filestruct.folder filesep filestruct.name]);
    gui = OspreyGUI(MRSCont);
    
    delete( gui.figure );
    
    %Test loading and plots
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
    Names = fieldnames(MRSCont.processed);
    for ss = 1 : length(Names)
        osp_plotModule(MRSCont, 'OspreyProcess', 1, Names{ss});
    end

        
    %Test fit and plots 
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
        

    osp_plotModule(MRSCont, 'OspreyCoreg', 1);



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