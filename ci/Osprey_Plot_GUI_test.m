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
        osp_plotModule(MRSCont, 'OspreyLoad', 1, [1 1], 'metabolites');
        if MRSCont.flags.hasRef
            osp_plotModule(MRSCont, 'OspreyLoad',1, [1 1], 'ref');
        end
        if MRSCont.flags.hasWater
            osp_plotModule(MRSCont, 'OspreyLoad', 1, [1 1], 'w');
        end
        if MRSCont.flags.hasMM
            osp_plotModule(MRSCont, 'OspreyLoad', 1, [1 1], 'MM');
        end

    %Test process and plots    
    Names = fieldnames(MRSCont.processed);
    for mm = 1 : length(Names)
        for ss = 1 : length(MRSCont.processed.(Names{mm}){1}.names)
            if (~contains(MRSCont.processed.(Names{mm}){1}.names{ss},'spline')) && (~contains(MRSCont.processed.(Names{mm}){1}.names{ss},'clean'))
                osp_plotModule(MRSCont, 'OspreyProcess', 1,[1 ss], Names{mm});
            end
        end
    end

        
    %Test fit and plots 
    Names = {'metab'}; 
    if MRSCont.flags.hasMM
        Names{end+1} = 'mm';
    end
    if MRSCont.flags.hasRef
        Names{end+1} = 'ref';
    end
    if MRSCont.flags.hasWater
        Names{end+1} = 'w';
    end

    for mm = 1 : length(Names)
        if isfield(MRSCont.fit.results,Names{mm})
            for bb = 1 : size(MRSCont.fit.results.(Names{mm}).fitParams,1)
                for ss = 1 : size(MRSCont.fit.results.(Names{mm}).fitParams,3)
                    osp_plotModule(MRSCont, 'OspreyFit', 1,[bb ss], Names{mm});
                end
            end
        end
    end
        

    osp_plotModule(MRSCont, 'OspreyCoreg', 1);

    SubNames = fieldnames(MRSCont.overview.SubSpecNamesStruct);
       k=1;
       if ~isempty(SubNames) 
           for i = 1 : length(SubNames)
               for j = 1 :size(MRSCont.overview.SubSpecNamesStruct.(SubNames{i}),2)
                    tempSubNames{k} = [SubNames{i}, ' ', MRSCont.overview.SubSpecNamesStruct.(SubNames{i}){1,j}]; 
                    k=k+1;
               end
           end
       end

       FitNames = fieldnames(MRSCont.overview.FitSpecNamesStruct);
       k=1;
       if ~isempty(FitNames) 
           for i = 1 : length(FitNames)
               for j = 1 :size(MRSCont.overview.FitSpecNamesStruct.(FitNames{i}),2)
                    tempFitNames{k} = ['Model ', FitNames{i}, ' ', MRSCont.overview.FitSpecNamesStruct.(FitNames{i}){1,j}]; 
                    k=k+1;
               end
           end
       end
       Names = [tempSubNames';tempFitNames'];
    for ss = 1 : length(Names)
        osp_plotModule(MRSCont, 'OspreySpecOverview', 1,1, Names{ss});
    end

    SubNames = fieldnames(MRSCont.overview.SubSpecNamesStruct);
       k=1;
       if ~isempty(SubNames) 
           for i = 1 : length(SubNames)
               for j = 1 :size(MRSCont.overview.SubSpecNamesStruct.(SubNames{i}),2)
                   if ~isempty(find(strcmp(FitNames,SubNames{i})))
                        if (~isempty(find(strcmp(MRSCont.overview.FitSpecNamesStruct.(SubNames{i}),MRSCont.overview.SubSpecNamesStruct.(SubNames{i}){1,j}))))
                            Names{k} = ['Model ' , SubNames{i}, ' ', MRSCont.overview.SubSpecNamesStruct.(SubNames{i}){1,j}];
                        else
                            Names{k} = [SubNames{i}, ' ', MRSCont.overview.SubSpecNamesStruct.(SubNames{i}){1,j}]; 
                        end
                   else
                       Names{k} = [SubNames{i}, ' ', MRSCont.overview.SubSpecNamesStruct.(SubNames{i}){1,j}]; 
                   end
                    k=k+1;
               end
           end
       end
    for ss = 1 : length(Names)                
        osp_plotModule(MRSCont, 'OspreyMeanOverview', 1,1, Names{ss});
    end

    if MRSCont.flags.isUnEdited
        osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, [1 1], 'metab-A-tCr', 'tNAA');
        osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, [1 1], 'metab-A-tCr', 'tCho');
        osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, [1 1], 'metab-A-tCr', 'mI');
        osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, [1 1], 'metab-A-tCr', 'Glx');

        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'tNAA', 'SNR');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'tCho', 'SNR');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'mI', 'SNR');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'Glx', 'SNR');

        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'tNAA', 'FWHM');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'tCho', 'FWHM');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'mI', 'FWHM');
        osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'Glx', 'FWHM');
    end
    if MRSCont.flags.isMEGA
            osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, [1 1], 'metab-diff1-tCr', MRSCont.opts.editTarget{1});
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-diff1-tCr', MRSCont.opts.editTarget{1}, 'SNR');
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-diff1-tCr', MRSCont.opts.editTarget{1}, 'FWHM');
    end
    if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
            osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, [1 1], 'metab-diff1-tCr', MRSCont.opts.editTarget{1});
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-diff1-tCr', MRSCont.opts.editTarget{1}, 'SNR');
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-diff1-tCr', MRSCont.opts.editTarget{1}, 'FWHM');
            osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, [1 1], 'metab-diff1-tCr', MRSCont.opts.editTarget{2});
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-diff1-tCr', MRSCont.opts.editTarget{2}, 'SNR');
            osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-diff1-tCr', MRSCont.opts.editTarget{2}, 'FWHM');
    end    
end