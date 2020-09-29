function osp_updatecorrOvWindow(gui) 
%% osp_updatecorrOvWindow
%   This function updates the correlation overview tab.
%
%
%   USAGE:
%       osp_updatecorrOvWindow(gui);
%
%   INPUT:  
%           gui      = gui class containing all handles and the MRSCont             
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
%       2020-01-16: First version of the code.
%%% 1. INITIALIZE %%%
            MRSCont = getappdata(gui.figure,'MRSCont');  % Get MRSCont from hidden container in gui class
            rmpath(genpath([gui.folder.spmversion filesep]));
            delete(gui.Plot.corrOv.Children(3).Children)
            Selection = gui.quant.popMenuNames{gui.quant.Selected.Quant};
            split_Selection = strsplit(Selection,'-');  
%%% 2. VISUALIZATION PART OF THIS TAB %%%
           if strcmp(split_Selection{2},'AlphaCorrWaterScaled') || strcmp(split_Selection{2},'AlphaCorrWaterScaledGroupNormed')
                metab = 'GABA';
            else
                metab = MRSCont.quantify.metabs.(split_Selection{1}){gui.overview.Selected.Metab};
           end 
            if gui.overview.Selected.CorrChoice == 1
                temp = osp_plotScatter(MRSCont,split_Selection{1},split_Selection{2},metab,gui.overview.CorrMeas{gui.overview.Selected.Corr},gui.overview.Names.Corr{gui.overview.Selected.Corr});
            else if gui.overview.Selected.CorrChoice == 2
                temp = osp_plotScatter(MRSCont,split_Selection{1},split_Selection{2},metab,MRSCont.quantify.metabs.(split_Selection{1}){gui.overview.Selected.Corr},MRSCont.quantify.metabs.(split_Selection{1}){gui.overview.Selected.Corr});
                else
                    switch gui.overview.Selected.Corr
                        case 1
                        temp = osp_plotScatter(MRSCont,split_Selection{1},split_Selection{2},metab,MRSCont.QM.SNR.A',gui.overview.Names.QM{gui.overview.Selected.Corr});
                        case 2
                        temp = osp_plotScatter(MRSCont,split_Selection{1},split_Selection{2},metab,MRSCont.QM.FWHM.A',gui.overview.Names.QM{gui.overview.Selected.Corr});
                    end
                end
            end
            set(gui.Plot.corrOv.Children(3), 'XLim', temp.Children(2).XLim);
            set(gui.Plot.corrOv.Children(3), 'YLim', temp.Children(2).YLim);
            set(temp.Children(2).Children, 'Parent', gui.Plot.corrOv.Children(3) );
            set(gui.Plot.corrOv.Children(3), 'XLabel', temp.Children(2).XLabel);
            set(gui.Plot.corrOv.Children(3), 'YLabel', temp.Children(2).YLabel);
            close(temp);
            set(gui.Plot.corrOv,'Heights', [-0.07 -0.90 -0.03]);
            gui.Plot.corrOv.Children(3).Legend.Location = 'North'; %Update legend
            if gui.overview.Selected.CorrChoice == 1
                set(gui.Plot.corrOv.Children(3).Title, 'String', [gui.overview.Names.Corr{gui.overview.Selected.Corr} ' vs ' metab]) %Update title
            else if gui.overview.Selected.CorrChoice == 2
                set(gui.Plot.corrOv.Children(3).Title, 'String', [MRSCont.quantify.metabs.(split_Selection{1}){gui.overview.Selected.Corr} ' vs ' metab]) %Update title
                else
                    switch gui.overview.Selected.Corr
                        case 1
                        set(gui.Plot.corrOv.Children(3).Title, 'String', [metab ' vs SNR']) %Update title
                        case 2
                        set(gui.Plot.corrOv.Children(3).Title, 'String', [metab ' vs FHWM (ppm)']) %Update title
                    end
                end
            end
        setappdata(gui.figure,'MRSCont',MRSCont); %Write  MRSCont into hidden container in gui class            
end