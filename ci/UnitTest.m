function [results,rt] = UnitTest(sequence,segmentation,GUItest)
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


%% Run jobFile with all Osprey modules
warning('off','all')
if segmentation == 0
    % Run fast unit test without segmentation
    results{1} = runtests(['Osp_' sequence '_CL_no_Seg.m']);
    rt{1} = table(results{1});

    if GUItest
        results{2} = runtests('Osprey_Plot_GUI_test.m');
        rt{2} = table(results{2});
    end

    if ~strcmp(sequence,'SinglePRESS')
        dir = strrep(which(['ci' filesep 'job' sequence '.m']),['job' sequence '.m'],'derivatives');
        rmdir(dir,'s')
    end

else if segmentation == 1
    % Run unit test with segmentation
    results{1} = runtests(['Osp_' sequence '_CL_w_Seg.m']);
    rt{1} = table(results{1});


    if GUItest
        results{2} = runtests('Osprey_Plot_GUI_test.m');
        rt{2} = table(results{2});
    end

    dir = strrep(which(['ci' filesep 'job' sequence '.m']),['job' sequence '.m'],'derivatives');
    rmdir(dir,'s')

    else
        % Run unit test with segmentation
    if   strcmp(sequence,'SinglePRESS')
        results{1} = runtests(['Osp_' sequence '_CL_Seg_downstream.m']);
        rt{1} = table(results{1});

        results{2} = runtests('Osprey_Plot_GUI_test.m');
        rt{2} = table(results{2});

        dir = strrep(which(['ci' filesep 'job' sequence '.m']),['job' sequence '.m'],'derivatives');
        rmdir(dir,'s')
    end

end



warning('on','all')
end
