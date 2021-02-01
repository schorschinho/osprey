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
        if ~strcmp(Selection,'Quality') 
            split_Selection = strsplit(Selection,'-');
    %This function updates the quantification table overview tab
            if strcmp(split_Selection{2},'AlphaCorrWaterScaled') || strcmp(split_Selection{2},'AlphaCorrWaterScaledGroupNormed')
                if isfield(MRSCont,'exclude') && ~isempty(MRSCont.exclude)
                    exclude = length(MRSCont.exclude);
                else
                    exclude = 0;
                end
                if ~strcmp(MRSCont.opts.fit.coMM3, 'none')
                    QuantTextOv = cell(MRSCont.nDatasets+1-exclude,2);
                    QuantTextOv(1,:) = {'GABA','GABA+'};
                else
                   QuantTextOv = cell(MRSCont.nDatasets+1-exclude,1);
                   QuantTextOv(1,:) = {'GABA'}; 
                end
                QuantTextOv(2:end,:) = table2cell(MRSCont.quantify.tables.(split_Selection{1}).(split_Selection{2}).Voxel_1(:,:));
            else
                if isfield(MRSCont,'exclude') && ~isempty(MRSCont.exclude)
                    exclude = length(MRSCont.exclude);
                else
                    exclude = 0;
                end            
                QuantTextOv = cell(MRSCont.nDatasets+1-exclude,length(MRSCont.quantify.metabs.(split_Selection{1})));
                QuantTextOv(1,:) = MRSCont.quantify.metabs.(split_Selection{1});
                QuantTextOv(2:end,:) = table2cell(MRSCont.quantify.tables.(split_Selection{1}).(split_Selection{2}).Voxel_1(:,:));
            end
        else
            if isfield(MRSCont,'exclude') && ~isempty(MRSCont.exclude)
                exclude = length(MRSCont.exclude);
            else
                exclude = 0;
            end  
            QuantTextOv = cell(MRSCont.nDatasets+1-exclude,length(MRSCont.QM.tables.Properties.VariableNames));
            QuantTextOv(1,:) = MRSCont.QM.tables.Properties.VariableNames;
            QuantTextOv(2:end,:) = table2cell(MRSCont.QM.tables(:,:));
        end
        tempF = figure( 'Visible', 'off' );
        temp=uimulticollist ( 'Parent',tempF,'units', 'normalized', 'position', [0 0 1 1], 'string', QuantTextOv,'Visible','off');
        set( temp, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
        if ~(isfield(MRSCont.flags,'isPRIAM') || isfield(MRSCont.flags,'isMRSI')) || ~(MRSCont.flags.isPRIAM || MRSCont.flags.isMRSI)
            set( temp, 'Parent', gui.Results.quantOv );
            close(tempF);
            if ~strcmp(Selection,'Quality') 
                set(gui.Results.quantOv, 'Title', ['Results: ' Selection]);
            else
                set(gui.Results.quantOv, 'Title', 'Results: Quality' );
            end
        else
            set( temp, 'Parent', gui.Results.quantOv1 );
            close(tempF);
            if ~strcmp(Selection,'Quality') 
                set(gui.Results.quantOv1, 'Title', ['Results Voxel 1: ' Selection]);
            else
                set(gui.Results.quantOv1, 'Title', 'Results Voxel 1: Quality' );
            end  
            
            % Set up a second table
                if ~strcmp(Selection,'Quality') 
                    split_Selection = strsplit(Selection,'-');
            %This function updates the quantification table overview tab
                    if strcmp(split_Selection{2},'AlphaCorrWaterScaled') || strcmp(split_Selection{2},'AlphaCorrWaterScaledGroupNormed')
                        if isfield(MRSCont,'exclude') && ~isempty(MRSCont.exclude)
                            exclude = length(MRSCont.exclude);
                        else
                            exclude = 0;
                        end
                        if ~strcmp(MRSCont.opts.fit.coMM3, 'none')
                            QuantTextOv = cell(MRSCont.nDatasets+1-exclude,2);
                            QuantTextOv(1,:) = {'GABA','GABA+'};
                        else
                           QuantTextOv = cell(MRSCont.nDatasets+1-exclude,1);
                           QuantTextOv(1,:) = {'GABA'}; 
                        end
                        QuantTextOv(2:end,:) = table2cell(MRSCont.quantify.tables.(split_Selection{1}).(split_Selection{2}).Voxel_2(:,:));
                    else
                        if isfield(MRSCont,'exclude') && ~isempty(MRSCont.exclude)
                            exclude = length(MRSCont.exclude);
                        else
                            exclude = 0;
                        end            
                        QuantTextOv = cell(MRSCont.nDatasets+1-exclude,length(MRSCont.quantify.metabs.(split_Selection{1})));
                        QuantTextOv(1,:) = MRSCont.quantify.metabs.(split_Selection{1});
                        QuantTextOv(2:end,:) = table2cell(MRSCont.quantify.tables.(split_Selection{1}).(split_Selection{2}).Voxel_2(:,:));
                    end
                else
                    if isfield(MRSCont,'exclude') && ~isempty(MRSCont.exclude)
                        exclude = length(MRSCont.exclude);
                    else
                        exclude = 0;
                    end  
                    QuantTextOv = cell(MRSCont.nDatasets+1-exclude,length(MRSCont.QM.tables.Properties.VariableNames));
                    QuantTextOv(1,:) = MRSCont.QM.tables.Properties.VariableNames;
                    QuantTextOv(2:end,:) = table2cell(MRSCont.QM.tables(:,:));
                end
                tempF = figure( 'Visible', 'off' );
                temp=uimulticollist ( 'Parent',tempF,'units', 'normalized', 'position', [0 0 1 1], 'string', QuantTextOv,'Visible','off');
                set( temp, 'BackgroundColor',gui.colormap.Background,'ForegroundColor', gui.colormap.Foreground);
                set( temp, 'Parent', gui.Results.quantOv2 );
                close(tempF);
                if ~strcmp(Selection,'Quality') 
                    set(gui.Results.quantOv2, 'Title', ['Results Voxel 2: ' Selection]);
                else
                    set(gui.Results.quantOv2, 'Title', 'Results Voxel 2: Quality' );
                end  
        end
            
end
