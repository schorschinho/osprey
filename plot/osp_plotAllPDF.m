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
            for kk = 1 : MRSCont.nDatasets
                osp_plotModule(MRSCont, 'OspreyLoad', kk, 'mets');
                if MRSCont.flags.hasRef
                    osp_plotModule(MRSCont, 'OspreyLoad', kk, 'ref');
                end
                if MRSCont.flags.hasWater
                    osp_plotModule(MRSCont, 'OspreyLoad', kk, 'w');
                end
                if MRSCont.flags.hasMM
                    osp_plotModule(MRSCont, 'OspreyLoad', kk, 'mm');
                end
            end
        case 'OspreyProcess'
            Names = fieldnames(MRSCont.processed);
            for kk = 1 : MRSCont.nDatasets
                for ss = 1 : length(Names)
                    osp_plotModule(MRSCont, 'OspreyProcess', kk, Names{ss});
                end
            end
        case 'OspreyFit'
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
            for kk = 1 : MRSCont.nDatasets
                for ss = 1 : length(Names)
                    osp_plotModule(MRSCont, 'OspreyFit', kk, Names{ss});
                end
            end
        case 'OspreyCoreg'
            for kk = 1 : MRSCont.nDatasets
              osp_plotModule(MRSCont, 'OspreyCoreg', kk);
            end
        case 'OspreySeg'   
             for kk = 1 : MRSCont.nDatasets
                 osp_plotModule(MRSCont, 'OspreySeg', kk);
             end 
        case 'OspreyOverview'
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
end