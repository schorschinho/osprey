function osp_updatequantOvWindow(gui)
%% osp_updatequantOvWindow
%   This function updates the quantify overview tab.
%
%
%   USAGE:
%       osp_updatequantOvWindow(gui);
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
        Selection = gui.quant.popMenuNames{gui.quant.Selected.Quant};
        split_Selection = strsplit(Selection,'-');
%This function updates the quantification table overview tab
        if strcmp(split_Selection{2},'AlphaCorrWaterScaled') || strcmp(split_Selection{2},'AlphaCorrWaterScaledGroupNormed')
                QuantTextOv = cell(MRSCont.nDatasets+1,1);
                QuantTextOv(1,:) = {'GABA'};
                QuantTextOv(2:end,:) = table2cell(MRSCont.quantify.tables.(split_Selection{1}).(split_Selection{2})(:,:));
        else
            QuantTextOv = cell(MRSCont.nDatasets+1,gui.quant.Number.Metabs);
            QuantTextOv(1,:) = MRSCont.quantify.metabs;
            QuantTextOv(2:end,:) = table2cell(MRSCont.quantify.tables.(split_Selection{1}).(split_Selection{2})(:,:));
        end
        temp=uimulticollist ( 'units', 'normalized', 'position', [0 0 1 1], 'string', QuantTextOv);
        set( temp, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        set( temp, 'Parent', gui.Results.quantOv );
        set(gui.Results.quantOv, 'Title', ['Results: ' (gui.quant.Names.Quants{gui.quant.Selected.Quant})]);
end
