function out = osp_plotInteractiveAtlasAnalysis(MRSCont,kk,quantification,metabolites,target_image,CurrLabel,AtlasThreshold,Interactive,convention)
%% out = osp_plotInteractiveAtlasAnalysis(MRSCont,kk,quantification,metabolites,target_image,CurrLabel,AtlasThreshold,Interactive,convention)
%   Creates a interactive figure to look at the atlas analysis results.
%
%   USAGE:
%       out = osp_plotInteractiveAtlasAnalysis(MRSCont,kk,quantification,metabolites,target_image,CurrLabel,AtlasThreshold,Interactive,convention)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       quantification = Which quantification to plot use 'TissCorrWaterScaled'
%       metabolites   = Target metabolites as cell {'tNAA','tCR'}
%       target_image = Target image to overlay on
%       CurrLabel = atlas region to start with in the plot
%       AtlasThreshold = partial fraction of MRSI volume threshold
%       Interactive     = flag for interactive mode
%       convention      = 'radiological' (L on right) or 'neurological' (L on left)
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.

% Fall back to defaults if not provided
if nargin < 9
    convention = 'neurological';
    if nargin < 8
        Interactive = 1;
        if nargin < 7
            AtlasThreshold = MRSCont.opts.MRSI.atlas.AtlasThreshold;
            if nargin < 6
                switch MRSCont.opts.MRSI.atlas.name
                    case 'AAL'
                        CurrLabel = 'Pallidum_L';
                    case 'neuromorphometrics'
                        CurrLabel = 'Left Pallidum';
                end        
                if nargin < 5
                    if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
                        target_image = 'MRSIloc_rMRSI'; 
                    else
                        target_image = 'T1w_rMRSI';
                    end
                    if nargin < 4
                        metabolites = MRSCont.opts.MRSI.atlas.metabolites;
                        if nargin < 3
                           quantification = 'amplitudes'; 
                           if nargin<2
                                kk = 1;
                                if nargin<1
                                    error('ERROR: no input Osprey container specified.  Aborting!!');
                                end
                           end
                        end
                    end
                end
            end
        end
    end
end

% Check that OspreyCoreg has been run before
if ~MRSCont.flags.didLoadData
    error('Trying to plot MRSI atlas analysis, but data has not been loaded yet. Run OspreyLoad first.')
end

if ~MRSCont.flags.didCoreg
    error('Trying to plot atlas analysis, but coregistration has not been performed yet. Run OspreyCoreg first.')
end
if ~MRSCont.flags.didFit
    error('Trying to plot MRSI atlas analysis, but data has not been modelled yet. Run OspreyFit first.')
end



%%% 1. PARSE INPUT ARGUMENTS %%%

[~,~,T1ext]   = fileparts(MRSCont.files_nii{1});
if strcmp(T1ext,'.gz')
    gunzip(MRSCont.files_nii{1});
    MRSCont.files_nii{1} = strrep(MRSCont.files_nii{1},'.gz','');
end

if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
    [~,~,MRSIlocext]   = fileparts(MRSCont.files_nii_MRSIloc{1});
    if strcmp(MRSIlocext,'.gz')
        gunzip(MRSCont.files_nii_MRSIloc{1});
        MRSCont.files_nii_MRSIloc{1} = strrep(MRSCont.files_nii_MRSIloc{1},'.gz','');
    end
end


% Get gifti vertices
g=gifti(MRSCont.gii_filename_VoxelGrid{1});
vertices = g.vertices;

NumberOfMetabolites = length(metabolites);

%% Load and prepare images based on target_image type
switch target_image
    case 'MRSIloc'
        [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_MRSIloc(MRSCont, vertices);
        
    case 'T1w_rMRSIloc'
        [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_T1w_rMRSIloc(MRSCont, vertices);
        
    case 'T1w_rMRSI'
        [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_T1w_rMRSI(MRSCont, vertices);
        
    case 'MRSIloc_rMRSI'
        [Coreg_img, AffineMat, vertices_voxel, MRSI_vol] = osp_load_MRSIloc_rMRSI(MRSCont, vertices);
        
    otherwise
        error('Unknown target_image type: %s', target_image);
end

%% Analyze orientation
[orientation_info] = osp_analyze_orientation(AffineMat);

%% Prepare image for display (with permutation and flips)
[Coreg_img, display_info] = osp_prepare_image_for_display(Coreg_img, orientation_info, convention);

%% Transform vertices to match displayed image
[vertices_display] = osp_transform_vertices_for_display(vertices_voxel, Coreg_img, orientation_info, display_info, MRSI_vol);
%%

max_val = prctile(Coreg_img(:), 99.5);

[VoxColors]=cbrewer('qual', 'Set1', 9);

% Get labels
switch MRSCont.opts.MRSI.atlas.name
    case 'AAL'
        atlas_path =  which('/libraries/AAL/AAL3v1_1mm.nii'); 
        load(which('/libraries/AAL/modROI_MNI_V7_1mm_List.mat'));
        AtlasLabels = {ROI.Nom_L};
        ROIid = [ROI.ID];
    case 'neuromorphometrics'
        atlas_path =  which('/libraries/neuromorphometrics/neuromorphometrics.nii');
        AtlasTable = readtable(which('/libraries/neuromorphometrics/neuromorphometrics.csv')); 
        ROIid = AtlasTable{:,1};
        AtlasLabels = AtlasTable{:,2};
end

CurrLabelIndex = find(strcmp(AtlasLabels,CurrLabel));
CurrQuantIndex = find(strcmp(fieldnames(MRSCont.quantify),quantification));

MRSCont.opts.MRSI.atlas.labels = AtlasLabels;
MRSCont.opts.MRSI.atlas.CurrLabelIndex = CurrLabelIndex;
MRSCont.opts.MRSI.atlas.CoregImage = Coreg_img;
MRSCont.opts.MRSI.atlas.kk = kk;
MRSCont.opts.MRSI.atlas.AtlasThreshold = AtlasThreshold;
MRSCont.opts.MRSI.atlas.vertices_voxel = vertices_voxel;
MRSCont.opts.MRSI.atlas.metabolites = metabolites;
MRSCont.opts.MRSI.atlas.quantification = quantification;

out = figure;  
set(out, 'units','normalized','outerposition',[0 0 1 1]);
set(out, 'units','pixel');
setappdata(out,'MRSCont',MRSCont);
tiledlayout(2, 2 + round(NumberOfMetabolites/2) ,'TileSpacing','compact')

if Interactive
 % Create ui menus
FracThreshText = uicontrol(  'Style', 'Text', ...
              'String', 'Atlas Voxel Fraction', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/5, out.InnerPosition(4)*0.03, 150, 30]);

SNRThreshText = uicontrol(  'Style', 'Text', ...
              'String', 'SNR > x', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/5+150, out.InnerPosition(4)*0.03, 150, 30]);

FWHMThreshText = uicontrol(  'Style', 'Text', ...
              'String', 'FWHM < x Hz', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/5+300, out.InnerPosition(4)*0.03, 150, 30]);

CRLBThreshText = uicontrol(  'Style', 'Text', ...
              'String', 'CRLB < x %', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/5+450, out.InnerPosition(4)*0.03, 150, 30]);

SDThreshText = uicontrol(  'Style', 'Text', ...
              'String', 'value < x*SD', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/5+600, out.InnerPosition(4)*0.03, 150, 30]);

FracThresh = uicontrol(  'Style', 'Edit', ...
              'String', num2str(AtlasThreshold), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/5, out.InnerPosition(4)*0.01, 150, 30]);

SNRThresh = uicontrol(  'Style', 'Edit', ...
              'String', num2str(MRSCont.opts.MRSI.atlas.SNRThreshold), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/5+150, out.InnerPosition(4)*0.01, 150, 30]);

FWHMThresh = uicontrol(  'Style', 'Edit', ...
              'String', num2str(MRSCont.opts.MRSI.atlas.FWHMThreshold), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/5+300, out.InnerPosition(4)*0.01, 150, 30]);

CRLBThresh = uicontrol(  'Style', 'Edit', ...
              'String', num2str(MRSCont.opts.MRSI.atlas.CRLBThreshold(1)), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/5+450, out.InnerPosition(4)*0.01, 150, 30]);

SDThresh = uicontrol(  'Style', 'Edit', ...
              'String', num2str(MRSCont.opts.MRSI.atlas.SDThreshold), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/5+600, out.InnerPosition(4)*0.01, 150, 30]);


UpdateButton = uicontrol(  'Style', 'pushbutton', ...
              'String', 'Update', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/5+750, out.InnerPosition(4)*0.01, 150, 30], ...
              'Callback', {@updatePlot,out,FracThresh,SNRThresh,FWHMThresh,CRLBThresh,SDThresh});

dropdownROI = uicontrol('Style', 'popupmenu', ...
                    'ForegroundColor',[11/255 71/255 111/255],...
                    'BackgroundColor',[1 1 1],...
                    'String', AtlasLabels, ...
                    'Position', [20, 20, 150, 30], ...
                    'Value',CurrLabelIndex,...
                    'Callback', {@updatePlot,out,FracThresh,SNRThresh,FWHMThresh,CRLBThresh,SDThresh});

dropdownQuant = uicontrol('Style', 'popupmenu', ...
                    'ForegroundColor',[11/255 71/255 111/255],...
                    'BackgroundColor',[1 1 1],...
                    'String', fieldnames(MRSCont.quantify), ...
                    'Position', [170, 20, 150, 30], ...
                    'Value', CurrQuantIndex,...
                    'Callback', {@updatePlot,out,FracThresh,SNRThresh,FWHMThresh,CRLBThresh,SDThresh});


markerCheckLR = uicontrol('Style', 'checkbox', ...
                            'ForegroundColor',[11/255 71/255 111/255],...
                            'BackgroundColor',[1 1 1],...
                           'String', 'Combine Left/Right', ...
                           'Position', [320, 20, 150, 30], ...
                           'Value', 0, ...
                           'Callback', {@updatePlot,out,FracThresh,SNRThresh,FWHMThresh,CRLBThresh,SDThresh});

end

AtlasMask = squeeze(MRSCont.atlas{kk}.fAtlas(CurrLabelIndex,:,:,:));
[~,AtlasMaskMaxIndex] = max(AtlasMask,[],'all');
[~, ~, zIndex] = ind2sub(size(AtlasMask), AtlasMaskMaxIndex);
AtlasMask(AtlasMask > AtlasThreshold) = 1;
AtlasMask(AtlasMask < AtlasThreshold) = 0;

AtlasMaskIndices = find(AtlasMask > 0);
[xAtlasMaskIndices, yAtlasMaskIndices, zAtlasMaskIndices] = ind2sub(size(AtlasMask), AtlasMaskIndices);
xAtlasMaskIndices(zAtlasMaskIndices~=zIndex)=[];
yAtlasMaskIndices(zAtlasMaskIndices~=zIndex)=[];
zAtlasMaskIndices(zAtlasMaskIndices~=zIndex)=[];
xAtlasMaskIndices = MRSCont.processed.A{kk}.sz(2)-xAtlasMaskIndices;
yAtlasMaskIndices = MRSCont.processed.A{kk}.sz(3)-yAtlasMaskIndices-1;



% Starting with the anatomical image and the atlas mask
nexttile(1)

imagesc(squeeze(Coreg_img(:,:,zIndex)),[0 max_val])
colormap gray;
axis image
hold on
ax = gca;
ax.XTickLabel = {};
ax.YTickLabel = {};
title(['Overlay ' AtlasLabels{CurrLabelIndex} ' fractional content > ' num2str(AtlasThreshold * 100) '%'],'Color',[11/255 71/255 111/255], 'Interpreter', 'none')

for vox = 1 : length(xAtlasMaskIndices)
    ind = (yAtlasMaskIndices(vox)*MRSCont.processed.A{kk}.sz(2)+xAtlasMaskIndices(vox)) * 8;
    tri = delaunay(double(vertices_voxel((1:4)+ind, 1)), double(vertices_voxel((1:4)+ind, 2)));
        trisurf(tri, vertices_voxel((1:4)+ind, 1), vertices_voxel((1:4)+ind, 2), vertices_voxel((1:4)+ind, 3),...
               'LineWidth',1,'FaceColor',VoxColors(1,:),'EdgeColor', 'none', 'FaceAlpha',0.60);
end

% All spectra from the voxels are plotted next
nexttile
hold on
tempSpec = zeros(MRSCont.processed.A{kk}.sz(1),length(xAtlasMaskIndices));
for vox = 1 : length(xAtlasMaskIndices)
    if ~MRSCont.flags.isMEGA
        tempSpec(:,vox) = squeeze(MRSCont.processed.A{kk}.specs(:,xAtlasMaskIndices(vox),yAtlasMaskIndices(vox),zAtlasMaskIndices(vox)));
    else
        tempSpec(:,vox) = squeeze(MRSCont.processed.diff1{kk}.specs(:,xAtlasMaskIndices(vox),yAtlasMaskIndices(vox),zAtlasMaskIndices(vox)));
    end
    plot(MRSCont.processed.A{kk}.ppm,real(tempSpec(:,vox)),'Color', [11/255 71/255 111/255]);
end
set(gca,'XDir','reverse','XLim',[0,4],'TickDir','out','YTickLabel',{},'YTick',{},'XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255])
xlabel('chemical shift (ppm)');
title(['Spectra from ' AtlasLabels{CurrLabelIndex}],'Color',[11/255 71/255 111/255], 'Interpreter', 'none')

nexttile(2 + round(NumberOfMetabolites/2) + 2)
hold on
plot(MRSCont.processed.A{kk}.ppm,real(mean(tempSpec,2)),'Color', [11/255 71/255 111/255]);
set(gca,'XDir','reverse','XLim',[0,4],'TickDir','out','YTickLabel',{},'YTick',{},'XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255])
xlabel('chemical shift (ppm)');
title(['Mean spectrum from ' AtlasLabels{CurrLabelIndex}],'Color',[11/255 71/255 111/255], 'Interpreter', 'none')


nexttile(2 + round(NumberOfMetabolites/2) + 1)
imagesc(flip((rot90(squeeze(AtlasMask(:,:,zIndex)))),2),[0 1])
colormap gray;
axis image
hold on
ax = gca;
ax.XTickLabel = {};
ax.YTickLabel = {};
title([AtlasLabels{CurrLabelIndex}  ' (' num2str(length(xAtlasMaskIndices)) ' total MRSI voxels)'],'Color',[11/255 71/255 111/255], 'Interpreter', 'none')


for mets = 1 : length(metabolites)
    if ~strcmp(metabolites{mets},'water') && ~strcmp(metabolites{mets},'CRLBs')
        plotMap = MRSCont.quantify.(quantification).(metabolites{mets});
    else
        plotMap = MRSCont.quantify.(metabolites{mets});
    end
    
    values = plotMap(AtlasMask == 1);
    mean_values = nanmean(values);
    std_values = nanstd(values);

    nexttile
    histogram(values,'FaceColor',VoxColors(mets,:));
    xline(mean_values,'--',{[num2str(round(mean_values,2)) '+-' num2str(round(std_values,2))]},'LineWidth',2,'Color',[11/255 71/255 111/255])
    title(metabolites{mets},'Color',[11/255 71/255 111/255], 'Interpreter', 'none')
    switch quantification
        case {'amplitudes','water'}
            xlabel('raw amplitudes')
        case {'CRLBs'}
            xlabel('relative CRLB (%)')
        case {'rawWaterScaled','CSFWaterScaled','TissCorrWaterScaled'}
            xlabel([quantification ' concentration (i.u.)'])
    end
    ylabel('Number of voxels')
    box off
    set(gca,'TickDir','out','XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255])
end
set(out, 'Color', [1 1 1]);

if strcmp(T1ext,'.gz')
    delete(MRSCont.files_nii{1});
    MRSCont.files_nii{1} = strrep(MRSCont.files_nii{1},'.nii','.nii.gz');
end

if isfield(MRSCont, 'MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
    if strcmp(MRSIlocext,'.gz')
        delete(MRSCont.files_nii_MRSIloc{1});
        MRSCont.files_nii_MRSIloc{1} = strrep(MRSCont.files_nii_MRSIloc{1},'.nii','.nii.gz');
     end
end


end
function updatePlot(source, ~,out,FracThresh,SNRThresh,FWHMThresh,CRLBThresh,SDThresh)
    MRSCont = getappdata(out,'MRSCont');  % Get MRSCont from hidden container in gui class  
    Coreg_img = MRSCont.opts.MRSI.atlas.CoregImage; 
    vertices_voxel = MRSCont.opts.MRSI.atlas.vertices_voxel;
    kk = MRSCont.opts.MRSI.atlas.kk;
    metabolites = MRSCont.opts.MRSI.atlas.metabolites;

    CombineLR = MRSCont.opts.MRSI.atlas.CombineLR;
    AtlasThreshold = str2num(FracThresh.String);
    SNRThreshold = str2num(SNRThresh.String);
    FWHMThreshold = str2num(FWHMThresh.String);
    CRLBThreshold = str2num(CRLBThresh.String);
    SDThreshold = str2num(SDThresh.String);

    if length(CRLBThreshold) > length(metabolites)
        CRLBThreshold = CRLBThreshold(1:length(metabolites));
    end
    if length(CRLBThreshold) < length(metabolites)
        if length(CRLBThreshold) == 1
            CRLBThreshold = repmat(CRLBThreshold,1,length(metabolites));
        else
            CRLBThreshold = cat(1,CRLBThreshold,repmat(CRLBThreshold(1),1,length(metabolites)-length(CRLBThreshold)));
        end
    end

    style = get(source,'Style');
    selection = source.Value;
    switch style
        case 'popupmenu'
            if length(source.String) > length(fieldnames(MRSCont.quantify))
                AtlasLabel = MRSCont.opts.MRSI.atlas.labels{selection};
                quantification = MRSCont.opts.MRSI.atlas.quantification;
                selectionROI = selection;
                MRSCont.opts.MRSI.atlas.CurrLabelIndex =selection;
            else
                quantification = fieldnames(MRSCont.quantify);
                quantification = quantification{selection};
                MRSCont.opts.MRSI.atlas.quantification = quantification;
                AtlasLabel = MRSCont.opts.MRSI.atlas.labels{MRSCont.opts.MRSI.atlas.CurrLabelIndex};
                selectionROI = MRSCont.opts.MRSI.atlas.CurrLabelIndex;
            end
        case 'checkbox'
            String = get(source,'String');
            if strcmp(String,'Combine Left/Right')
                CombineLR = selection;
                MRSCont.opts.MRSI.atlas.CombineLR = selection;
            end
            quantification = MRSCont.opts.MRSI.atlas.quantification;
            AtlasLabel = MRSCont.opts.MRSI.atlas.labels{MRSCont.opts.MRSI.atlas.CurrLabelIndex};
            selectionROI = MRSCont.opts.MRSI.atlas.CurrLabelIndex;
        otherwise
            quantification = MRSCont.opts.MRSI.atlas.quantification;
            AtlasLabel = MRSCont.opts.MRSI.atlas.labels{MRSCont.opts.MRSI.atlas.CurrLabelIndex};
            selectionROI = MRSCont.opts.MRSI.atlas.CurrLabelIndex;

    end
    
   
    [VoxColors]=cbrewer('qual', 'Set1', 9);

    

    if CombineLR      
        AtlasLabel = AtlasLabel(1:end-2);
        selectionROI = find(contains(MRSCont.opts.MRSI.atlas.labels,AtlasLabel));
        AtlasMask = squeeze(sum(MRSCont.atlas{kk}.fAtlas(selectionROI,:,:,:),1));
    else
        AtlasMask = squeeze(MRSCont.atlas{kk}.fAtlas(selectionROI,:,:,:));    
    end

    % Now lets apply other thresholds
    AtlasMask(MRSCont.quickMaps.A.FWHM > FWHMThreshold) = 0;
    AtlasMask(MRSCont.quickMaps.A.SNR < SNRThreshold) = 0;

    % For CRLBs we need it per metabolite
    for m = 1 : length(metabolites)
        tempAtlasMask = AtlasMask;
        tempAtlasMask(MRSCont.quantify.CRLBs.(metabolites{m}) > CRLBThreshold(m)) = 0;
        AtlasMaskMetabs{m} = tempAtlasMask;
    end

    
    [~,AtlasMaskMaxIndex] = max(AtlasMask,[],'all');
    [~, ~, zIndex] = ind2sub(size(AtlasMask), AtlasMaskMaxIndex);
    AtlasMask(AtlasMask > AtlasThreshold) = 1;
    AtlasMask(AtlasMask < AtlasThreshold) = 0;

    for m = 1 : length(metabolites)
        tempAtlasMask = AtlasMaskMetabs{m};
        tempAtlasMask(tempAtlasMask > AtlasThreshold) = 1;
        tempAtlasMask(tempAtlasMask < AtlasThreshold) = 0;
        AtlasMaskMetabs{m} = tempAtlasMask;
    end

    AtlasMaskIndices = find(AtlasMask > 0);
    [xAtlasMaskIndices, yAtlasMaskIndices, zAtlasMaskIndices] = ind2sub(size(AtlasMask), AtlasMaskIndices);
    xAtlasMaskIndices(zAtlasMaskIndices~=zIndex)=[];
    yAtlasMaskIndices(zAtlasMaskIndices~=zIndex)=[];
    zAtlasMaskIndices(zAtlasMaskIndices~=zIndex)=[];
    xAtlasMaskIndices = MRSCont.processed.A{kk}.sz(2)-xAtlasMaskIndices;
    yAtlasMaskIndices = MRSCont.processed.A{kk}.sz(3)-yAtlasMaskIndices-1;


    AllAxes = out.Children(end).Children;
    max_val = prctile(Coreg_img(:), 99.5);

    MaxIndexAx = length(out.Children(end).Children);
    MetabIndex = 1;
    

    for aa = length(out.Children(end).Children) : -1 : 1
        cla(AllAxes(aa));
        set(out,'CurrentAxes',AllAxes(aa))
        switch aa
            case MaxIndexAx
        
                imagesc(squeeze(Coreg_img(:,:,zIndex)),[0 max_val]);
                   
                if ~isempty(AtlasMaskIndices)
                    for vox = 1 : length(xAtlasMaskIndices)
                        ind = (yAtlasMaskIndices(vox)*MRSCont.processed.A{kk}.sz(2)+xAtlasMaskIndices(vox)) * 8;
                        tri = delaunay(double(vertices_voxel((1:4)+ind, 1)), double(vertices_voxel((1:4)+ind, 2)));
                        trisurf(tri, vertices_voxel((1:4)+ind, 1), vertices_voxel((1:4)+ind, 2), vertices_voxel((1:4)+ind, 3),...
                               'LineWidth',1,'FaceColor',VoxColors(1,:),'EdgeColor', 'none', 'FaceAlpha',0.60);
                    end
                end
                title(['Overlay ' AtlasLabel ' fractional content > ' num2str(AtlasThreshold * 100) '%'],'Color',[11/255 71/255 111/255], 'Interpreter', 'none')
            case MaxIndexAx - 1
                if ~isempty(AtlasMaskIndices)
                    tempSpec = zeros(MRSCont.processed.A{kk}.sz(1),length(xAtlasMaskIndices));
                    for vox = 1 : length(xAtlasMaskIndices)
                        if ~MRSCont.flags.isMEGA
                            tempSpec(:,vox) = squeeze(MRSCont.processed.A{kk}.specs(:,xAtlasMaskIndices(vox),yAtlasMaskIndices(vox),zAtlasMaskIndices(vox)));
                        else
                            tempSpec(:,vox) = squeeze(MRSCont.processed.diff1{kk}.specs(:,xAtlasMaskIndices(vox),yAtlasMaskIndices(vox),zAtlasMaskIndices(vox)));
                        end
                        plot(MRSCont.processed.A{kk}.ppm,real(tempSpec(:,vox)),'Color', [11/255 71/255 111/255]);
                    end
                else
                    text(3.9,0.1,'No MRSI voxel found for this region and threshold','Color',[11/255 71/255 111/255])
                end
                title(['Spectra from ' AtlasLabel],'Color',[11/255 71/255 111/255], 'Interpreter', 'none')
            case MaxIndexAx - 2
                if ~isempty(AtlasMaskIndices)
                    tempSpec = zeros(MRSCont.processed.A{kk}.sz(1),length(xAtlasMaskIndices));
                    for vox = 1 : length(xAtlasMaskIndices)
                        if ~MRSCont.flags.isMEGA
                            tempSpec(:,vox) = squeeze(MRSCont.processed.A{kk}.specs(:,xAtlasMaskIndices(vox),yAtlasMaskIndices(vox),zAtlasMaskIndices(vox)));
                        else
                            tempSpec(:,vox) = squeeze(MRSCont.processed.diff1{kk}.specs(:,xAtlasMaskIndices(vox),yAtlasMaskIndices(vox),zAtlasMaskIndices(vox)));
                        end
                    end
                    plot(MRSCont.processed.A{kk}.ppm,real(mean(tempSpec,2)),'Color', [11/255 71/255 111/255]);
                else
                    text(3.9,0.1,'No MRSI voxel found for this region and threshold','Color',[11/255 71/255 111/255])
                end
                title(['Mean spectrum from ' AtlasLabel],'Color',[11/255 71/255 111/255], 'Interpreter', 'none')
            case MaxIndexAx - 3
                imagesc(flip(rot90(squeeze(AtlasMask(:,:,zIndex))),2),[0 1])
                title([AtlasLabel ' (' num2str(length(xAtlasMaskIndices)) ' total MRSI voxels)'],'Color',[11/255 71/255 111/255], 'Interpreter', 'none')
            otherwise
                if ~isempty(AtlasMaskIndices)
                    if ~strcmp(metabolites{MetabIndex},'water') && ~strcmp(metabolites{MetabIndex},'CRLBs')
                        plotMap = MRSCont.quantify.(quantification).(metabolites{MetabIndex});
                    else
                        plotMap = MRSCont.quantify.(metabolites{MetabIndex});
                    end
                    
                    values = plotMap(AtlasMaskMetabs{MetabIndex} == 1);
                    mean_values = nanmean(values);
                    std_values = nanstd(values);

                    values(values > (mean_values + (SDThreshold*std_values))) = [];
                    values(values < (mean_values - (SDThreshold*std_values))) = [];

                    mean_values = nanmean(values);
                    std_values = nanstd(values);
    
                    histogram(values,'FaceColor',VoxColors(MetabIndex,:));
                    xline(mean_values,'--',{[num2str(round(mean_values,2)) '+-' num2str(round(std_values,2))]},'LineWidth',2,'Color',[11/255 71/255 111/255])
                    title(metabolites{MetabIndex},'Color',[11/255 71/255 111/255], 'Interpreter', 'none')
                    set(gca,'TickDir','out','XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255])
                    box off
                    switch quantification
                        case {'amplitudes','water'}
                            xlabel('raw amplitudes')
                        case {'CRLBs'}
                            xlabel('relative CRLB (%)')
                        case {'rawWaterScaled','CSFWaterScaled','TissCorrWaterScaled'}
                            xlabel([quantification ' concentration (i.u.)'])
                    end
                else
                    text(1.55,0.1,'No MRSI voxel found for this region and threshold','Color',[11/255 71/255 111/255])
                    set(gca,'XLim',[1.5 2])
                end
                MetabIndex = MetabIndex + 1;

        end

    end
    setappdata(out,'MRSCont',MRSCont);
  
end
   