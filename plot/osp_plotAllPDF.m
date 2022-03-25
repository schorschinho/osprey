function osp_plotAllPDF(MRSCont, Module)
%% osp_plotModule
%   Callback function on print figure button click.
%
%
%   USAGE:
%       osp_plotModule(MRSCont, Module)
%
%   INPUT:     MRSCont  = Osprey data container.
%              Module       = String for the Module     
%              OPTIONS:    - 'OspreyLoad' (default)
%                          - 'OspreyProcess'
%                          - 'OspreyFit'
%                          - 'OspreyCoreg'
%                          - 'OspreySeg'
%                          - 'OspreyOverview'
%   OUTPUT:     all figures
%
%
%   AUTHORS:
%       Dr. Helge Zoellner (Johns Hopkins University, 2021-02-12)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions from the FID-A toolbox by
%       Dr. Jamie Near (McGill University)
%       https://github.com/CIC-methods/FID-A
%       Simpson et al., Magn Reson Med 77:23-33 (2017)
%
%   HISTORY:
%       2021-02-12: First version of the code.
%%% 1. Create plots accroding to the module%%%
if  ~MRSCont.flags.isPRIAM && ~MRSCont.flags.isMRSI
    switch Module
        case 'OspreyLoad'
            for kk = 1 : MRSCont.nDatasets(1)
                osp_plotModule(MRSCont, 'OspreyLoad', kk,[1 1], 'metabolites');
                if MRSCont.flags.hasRef
                    osp_plotModule(MRSCont, 'OspreyLoad', kk,[1 1], 'ref');
                end
                if MRSCont.flags.hasWater
                    osp_plotModule(MRSCont, 'OspreyLoad', kk,[1 1], 'w');
                end
                if MRSCont.flags.hasMM
                    osp_plotModule(MRSCont, 'OspreyLoad', kk,[1 1], 'MM');
                end
            end
        case 'OspreyProcess'
            Names = fieldnames(MRSCont.processed);
            for kk = 1 : MRSCont.nDatasets(1)
                for mm = 1 : length(Names)
                    for ss = 1 : length(MRSCont.processed.(Names{mm}){kk}.names)
                        if (~contains(MRSCont.processed.(Names{mm}){kk}.names{ss},'spline')) && (~contains(MRSCont.processed.(Names{mm}){kk}.names{ss},'clean'))
                            osp_plotModule(MRSCont, 'OspreyProcess', kk,[1 ss], Names{mm});
                        end
                    end
                end
            end
        case 'OspreyFit'
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
                
            for kk = 1 : MRSCont.nDatasets(1)
                for mm = 1 : length(Names)
                    if isfield(MRSCont.fit.results,Names{mm})
                        for bb = 1 : size(MRSCont.fit.results.(Names{mm}).fitParams,1)
                            for ss = 1 : size(MRSCont.fit.results.(Names{mm}).fitParams,3)
                                osp_plotModule(MRSCont, 'OspreyFit', kk,[bb ss], Names{mm});
                            end
                        end
                    end
                end
            end
        case 'OspreyCoreg'
            for kk = 1 : MRSCont.nDatasets(1)
              osp_plotModule(MRSCont, 'OspreyCoreg', kk);
            end
        case 'OspreySeg'   
             for kk = 1 : MRSCont.nDatasets(1)
                 osp_plotModule(MRSCont, 'OspreySeg', kk);
             end 
        case 'OspreyOverview'
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
                osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, [1 1], 'metab-A-tCr', 'Ins');
                osp_plotModule(MRSCont, 'OspreyRaincloudOverview', 1, [1 1], 'metab-A-tCr', 'Glx');

                osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'tNAA', 'SNR');
                osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'tCho', 'SNR');
                osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'Ins', 'SNR');
                osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'Glx', 'SNR');

                osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'tNAA', 'FWHM');
                osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'tCho', 'FWHM');
                osp_plotModule(MRSCont, 'OspreyScatterOverview', 1, [1 1], 'metab-A-tCr', 'Ins', 'FWHM');
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
end