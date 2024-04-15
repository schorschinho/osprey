function osp_onSave( ~, ~,gui)
%% osp_onSave
%   Callback function on save button click. Saves MRSCont.
%
%
%   USAGE:
%       osp_onSave( ~, ~ ,gui);
%
%   INPUT:      gui      = gui class containing all handles and the MRSCont 
%
%   OUTPUT:     Changes in gui parameters and MRSCont are written into the
%               gui class
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2020-01-16)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2020-01-22: First version of the code.
%%% 1. SAVE %%%
    % Save the output structure to the output folder
    % Determine output folder
    refProcessTime = tic;
    fprintf('Saving MRS container...\n');
    MRSCont = getappdata(gui.figure,'MRSCont'); % Get MRSCont from hidden container in gui class  
    set(gui.figure,'HandleVisibility','off');
    outputFolder    = MRSCont.outputFolder;
    outputFile      = MRSCont.outputFile;
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end
    MRSCont.flags.isGUI = 0;
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    MRSCont.flags.isGUI = 1;
    fprintf('... done.\n');
    setappdata(gui.figure,'MRSCont',MRSCont); % Write MRSCont into hidden container in gui class
    set(gui.figure,'HandleVisibility','on');
    toc(refProcessTime);

    if gui.process.ManualManipulation
        fprintf('Rewriting MRS data after manual manipulation...\n');
        % Optional: write edited files to LCModel .RAW files
        if MRSCont.opts.saveLCM && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
            [MRSCont] = osp_saveLCM(MRSCont);
        end
        
        % Optional: write edited files to jMRUI .txt files
        if MRSCont.opts.savejMRUI && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
            [MRSCont] = osp_saveJMRUI(MRSCont);
        end
        
        % Optional: write edited files to vendor specific format files readable to
        % LCModel and jMRUI
        % SPAR/SDAT if Philips
        % RDA if Siemens
        if MRSCont.opts.saveVendor && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
            [MRSCont] = osp_saveVendor(MRSCont);
        end
        
        % Optional: write edited files to NIfTI-MRS format
        if MRSCont.opts.saveNII && ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
            [MRSCont] = osp_saveNII(MRSCont);
        end
        fprintf('... done.\n');
    end

end % onExit