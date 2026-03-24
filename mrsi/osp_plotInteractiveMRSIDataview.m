function out = osp_plotInteractiveMRSIDataview(MRSCont, target_image, target_module,FitTarget,FitPlotFunction,convention)
%% out = osp_plotInteractiveMRSIDataview(MRSCont, target, kk)
%   Creates a interactive figure to define a maks to remove voxels that are
%   not interesting from the fitting as well as the segmentation.
%
%   USAGE:
%       out = osp_plotInteractiveMRSIDataview(MRSCont, kk)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.

% Fall back to defaults if not provided
if nargin < 6
convention = 'neurological';
    if nargin < 5
        FitPlotFunction = 'Fit1DStack';
        if nargin < 4
            switch target_module
                case 'OspreyLoad'
                    FitTarget = 'raw';
                case 'OspreyProcess'
                    FitTarget = 'A';
                case 'OspreyFit'
                    FitTarget = 'metab';
            end            
            if nargin < 3
                target_module = 'OspreyLoad';
                if nargin < 2
                   if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
                        target_image = 'MRSIloc_rMRSI'; 
                    else
                        target_image = 'T1w_rMRSI';
                    end 
                    if nargin<1
                        error('ERROR: no input Osprey container specified.  Aborting!!');
                    end
                end
            end
        end
    end
end


%% Validate prerequisites
if ~MRSCont.flags.didLoadData
    error('Trying to visualize MRSI data, but data has not been loaded yet. Run OspreyLoad first.')
end

if ~MRSCont.flags.didCoreg
    error('Trying to visualize MRSI data, but MRSI grid mask has not been created yet. Run OspreyCoreg first.')
end

if strcmp(target_module,'OspreyProcess') && ~MRSCont.flags.didProcess
    error('Trying to visualize processed MRSI data, but data has not been loaded yet. Run OspreyProcess first.')
end

if strcmp(target_module,'OspreyFit') && ~MRSCont.flags.didFit
    error('Trying to visualize MRSI fits, but data has not been modeled yet. Run OspreyFit first.')
end


%% Handle compressed files
[~, ~, T1ext] = fileparts(MRSCont.files_nii{1});
if strcmp(T1ext, '.gz')
    gunzip(MRSCont.files_nii{1});
    MRSCont.files_nii{1} = strrep(MRSCont.files_nii{1}, '.gz', '');
end

MRSIlocext = '';
if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
    [~, ~, MRSIlocext] = fileparts(MRSCont.files_nii_MRSIloc{1});
    if strcmp(MRSIlocext, '.gz')
        gunzip(MRSCont.files_nii_MRSIloc{1});
        MRSCont.files_nii_MRSIloc{1} = strrep(MRSCont.files_nii_MRSIloc{1}, '.gz', '');
    end
end

%% Setup colors
[VoxColors] = cbrewer('qual', 'Set1', 9);

%% Get gifti vertices
g = gifti(MRSCont.gii_filename_VoxelGrid{1});
vertices = g.vertices;

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
[Coreg_img_display, display_info] = osp_prepare_image_for_display(Coreg_img, orientation_info, convention);

%% Transform vertices to match displayed image
[vertices_display] = osp_transform_vertices_for_display(vertices_voxel, Coreg_img, orientation_info, display_info, MRSI_vol);

%% Setup figure
out = figure;  
tiledlayout(1, 2 ,'TileSpacing','compact')
max_val = prctile(Coreg_img_display(:), 99.5);


slice = round(MRSCont.raw{1}.nZvoxels/2);

set(out, 'units','normalized','outerposition',[0 0 1 1]);
set(out, 'units','pixel');
set(out, 'Color', [1 1 1]);

nexttile

% Prepare image 
Coreg_img_3D = Coreg_img_display;
Coreg_img = squeeze(Coreg_img_3D(:,:,slice));

imagesc(Coreg_img,[0 max_val])

%% Add orientation labels
tile_height = size(Coreg_img_display, 1);
tile_width = size(Coreg_img_display, 2);

osp_add_orientation_labels_to_montage(1, 1, 1, ...
    tile_width, tile_height, orientation_info, display_info);
colormap gray;
axis image
axis off
box off



hold on;


MRSCont.opts.MRSI.InteractiveDataView.CoregImage = Coreg_img_3D;
MRSCont.opts.MRSI.InteractiveDataView.vertices_voxel = vertices_display;
MRSCont.opts.MRSI.InteractiveDataView.kk = 1;
MRSCont.opts.MRSI.InteractiveDataView.currVoxel = [round(MRSCont.raw{1}.nXvoxels/2), round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2)];
MRSCont.opts.MRSI.InteractiveDataView.target_module = target_module;
MRSCont.opts.MRSI.InteractiveDataView.FitPlotFunction = FitPlotFunction;
MRSCont.opts.MRSI.InteractiveDataView.target_module_labels = {'OspreyLoad','OspreyProcess','OspreyFit'};
MRSCont.opts.MRSI.InteractiveDataView.FitTarget = FitTarget;
MRSCont.opts.MRSI.InteractiveDataView.orientation_info = orientation_info;
MRSCont.opts.MRSI.InteractiveDataView.display_info = display_info;




setappdata(out,'MRSCont',MRSCont);


VoxelIndices = [round(MRSCont.raw{1}.nXvoxels/2), round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2)];
VoxelIndices(:,1:2) = VoxelIndices(:,1:2)-1;

for vox = 1 : size(VoxelIndices,1)
    ind = (VoxelIndices(vox,2)*MRSCont.raw{1}.nXvoxels+VoxelIndices(vox,1)) * 8;
    tri = delaunay(double(vertices_voxel((1:4)+ind, 2)), double(vertices_voxel((1:4)+ind, 1)));
    trisurf(tri, vertices_voxel((1:4)+ind, 2), vertices_voxel((1:4)+ind, 1), vertices_voxel((1:4)+ind, 3),...
           'LineWidth',1,'FaceColor',[255/255 140/255 0/255],'EdgeColor', 'none', 'FaceAlpha',0.85);
end


% addGrid    
% Find vertices near the current slice (in voxel space, z-axis)
valid_vertices = abs(vertices_display(:, 3) - slice) < 1; 

% Transform voxel coordinates to montage coordinates        
plot(vertices_display(valid_vertices, 2) , ...
     vertices_display(valid_vertices, 1),'.','Color',[254/255 186/255 47/255], 'MarkerSize', 6);

axis image
title(['MRSI voxel (' num2str(round(MRSCont.raw{1}.nXvoxels/2)) ',' num2str(round(MRSCont.raw{1}.nYvoxels/2)) ',' num2str(round(MRSCont.raw{1}.nZvoxels/2)) ')'],...
    'Color',[11/255 71/255 111/255], 'Interpreter', 'none');


nexttile

switch target_module
    case 'OspreyLoad'
        if ~MRSCont.flags.isMEGA 
            plot(MRSCont.raw{1}.ppm,...
                squeeze(real(MRSCont.raw{1}.specs(:,round(MRSCont.raw{1}.nXvoxels/2), round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2)))), 'Color', [11/255 71/255 111/255],'LineWidth',1.5);
        else
            plot(MRSCont.raw{1}.ppm,...
                squeeze(real(MRSCont.raw{1}.specs(:,1,round(MRSCont.raw{1}.nXvoxels/2), round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2)))), 'Color', [11/255 71/255 111/255],'LineWidth',1.5);
        end
        set(gca,'XLim',[0.5,4.2],'XDir','reverse','TickDir','out','XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255]);
        xlabel('chemical shift (ppm)');
        box off
    case 'OspreyProcess'
        plot(MRSCont.processed.(FitTarget){1}.ppm,...
            squeeze(real(MRSCont.processed.(FitTarget){1}.specs(:,round(MRSCont.raw{1}.nXvoxels/2), round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2)))), 'Color', [11/255 71/255 111/255],'LineWidth',1.5);
        set(gca,'XLim',[0.5,4.2],'XDir','reverse','TickDir','out','XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255]);
        xlabel('chemical shift (ppm)');
        box off
    case 'OspreyFit'
        ModelMatrix = MRSCont.fit.(FitTarget);
        switch FitPlotFunction
            case 'Fit1DStack'
                ModelMatrix{round(MRSCont.raw{1}.nXvoxels/2),round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2)}.plotFit1DStack(0)
            case 'Fit3D'
                ModelMatrix{round(MRSCont.raw{1}.nXvoxels/2),round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2)}.plotFit3D(0)
        end

end

uiStruct.LeftButton = uicontrol(  'Style', 'pushbutton', ...
              'String', '<', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)/4 - 70, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.zVoxelInd = uicontrol(  'Style', 'Edit', ...
              'String', round(MRSCont.raw{1}.nZvoxels/2), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)/4, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.RightButton = uicontrol(  'Style', 'pushbutton', ...
              'String', '>', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)/4 + 70, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.target_module = uicontrol('Style', 'popupmenu', ...
                    'String', {'OspreyLoad','OspreyProcess','OspreyFit'}, ...
                    'Position', [out.InnerPosition(3)/4 + 140, out.InnerPosition(4)*0.03, 150, 30]);


uiStruct.Textppm1 = uicontrol(  'Style', 'Text', ...
              'String', 'ppm(1)', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/4-60, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.Textppm2 = uicontrol(  'Style', 'Text', ...
              'String', 'ppm(2)', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/4 + 10, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.ppm1 = uicontrol(  'Style', 'Edit', ...
              'String', '0.5', ...
               'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/4-60, out.InnerPosition(4)*0.01, 60, 30]);

uiStruct.ppm2 = uicontrol(  'Style', 'Edit', ...
              'String', '4.2', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/4+10, out.InnerPosition(4)*0.01, 60, 30]);


uiStruct.TextIntensity1 = uicontrol(  'Style', 'Text', ...
              'String', 'Intensity(1)', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/4 + 70, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.TextIntensity2 = uicontrol(  'Style', 'Text', ...
              'String', 'Intensity(2)', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/4 + 130, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.Intensity1 = uicontrol(  'Style', 'Edit', ...
              'String', 'NaN', ...
               'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/4+70, out.InnerPosition(4)*0.01, 60, 30]);

uiStruct.Intensity2 = uicontrol(  'Style', 'Edit', ...
              'String', 'NaN', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/4+130, out.InnerPosition(4)*0.01, 60, 30]);

uiStruct.TextZerofill = uicontrol(  'Style', 'Text', ...
              'String', 'Zero-filling factor', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/4 + 200, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.TextLB = uicontrol(  'Style', 'Text', ...
              'String', 'Gaussian LB', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/4 + 260, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.Zerofill = uicontrol(  'Style', 'Edit', ...
              'String', '0', ...
               'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/4+200, out.InnerPosition(4)*0.01, 60, 30]);

uiStruct.LB = uicontrol(  'Style', 'Edit', ...
              'String', '0', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/4+260, out.InnerPosition(4)*0.01, 60, 30]);

uiStruct.TextSubSpec = uicontrol(  'Style', 'Text', ...
              'String', 'Subspec #', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*2/4 + 320, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.SubSpec = uicontrol(  'Style', 'Edit', ...
              'String', '1', ...
               'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*2/4+320, out.InnerPosition(4)*0.01, 60, 30]);

uiStruct.Magnitude = uicontrol( 'Style', 'checkbox', ...
              'String', 'Magnitude Spectrum', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Value', 0,...
              'Position', [out.InnerPosition(3)*2/4+380, out.InnerPosition(4)*0.01, 150, 30]);

set(uiStruct.LeftButton,'Callback', {@updatePlot,out,uiStruct})
set(uiStruct.RightButton,'Callback', {@updatePlot,out,uiStruct})
set(uiStruct.target_module,'Callback', {@updatePlot,out,uiStruct})

% Set up the click callback function
set(gcf, 'WindowButtonDownFcn',{@mouseClick,out,uiStruct});

    [img, ~, ~] = imread('osprey.png', 'BackgroundColor', [1,1,1]);
    [img2] = imresize(img, 0.08);
    ax2 = axes( ...
    'Units','normalized', ...
    'Position',[0.9 0.85 0.05 0.05], ...
    'YDir','reverse', ...
    'Visible','off');
    image(img2)
    box off
    axis image
    axis off

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

% Callback function for mouse clicks
function mouseClick(~, ~,out,uiStruct)
    MRSCont = getappdata(out,'MRSCont');  % Get MRSCont from hidden container in gui class 
    Coreg_img_3D = MRSCont.opts.MRSI.InteractiveDataView.CoregImage;
    Magnitude = uiStruct.Magnitude.Value;
    ppm1 = str2num(uiStruct.ppm1.String);
    ppm2 = str2num(uiStruct.ppm2.String);
    Intensity1 = str2num(uiStruct.Intensity1.String);
    Intensity2 = str2num(uiStruct.Intensity2.String);
    zerofill = str2num(uiStruct.Zerofill.String);
    GaussianLB = str2num(uiStruct.LB.String);
    Subspectrum = str2num(uiStruct.SubSpec.String);
    target_module = MRSCont.opts.MRSI.InteractiveDataView.target_module;
    FitPlotFunction = MRSCont.opts.MRSI.InteractiveDataView.FitPlotFunction;
    FitTarget= MRSCont.opts.MRSI.InteractiveDataView.FitTarget;
    orientation_info = MRSCont.opts.MRSI.InteractiveDataView.orientation_info;
    display_info = MRSCont.opts.MRSI.InteractiveDataView.display_info;

    
    
    % Get the current point
    currentPoint = get(gca, 'CurrentPoint');
    x = round(currentPoint(1, 1));
    y = round(currentPoint(1, 2));
    slice = str2num(uiStruct.zVoxelInd.String);

    vertices_voxel = MRSCont.opts.MRSI.InteractiveDataView.vertices_voxel;
    valid_vertices = abs(vertices_voxel(:, 3) - slice) < 1; 
    vertices_voxel_click = vertices_voxel;
    vertices_voxel_click(~valid_vertices,3)=MRSCont.raw{1}.nZvoxels*10;
    distances = sqrt((vertices_voxel_click(:,2) - x).^2 + (vertices_voxel_click(:,1) - y).^2 + (vertices_voxel_click(:,3) - slice).^2 );
    distances = sum(reshape(distances,8,[]),1);
    [~, row_index] = min(distances);
    [MRSI_x,MRSI_y,~]=ind2sub([MRSCont.raw{1}.nXvoxels, MRSCont.raw{1}.nYvoxels,MRSCont.raw{1}.nZvoxels],row_index);
    MRSI_z = slice;
    
    % Display click coordinates for debugging
    % fprintf('Clicked at: (%.2f, %.2f)\n', x, y);
    % fprintf('MRSI voxel: (%.2f, %.2f, %.2f)\n', MRSI_x, MRSI_y, MRSI_z);
    
    AllAxes = out.Children(end).Children;
    cla(AllAxes(2));
    set(out,'CurrentAxes',AllAxes(2))
    max_val = prctile(Coreg_img_3D(:), 99.5);

    % Prepare image 
    Coreg_img = squeeze(Coreg_img_3D(:,:,slice));
    

    imagesc(Coreg_img,[0 max_val])
    colormap gray;
    axis image
    
    hold on;
    
    MRSCont.opts.MRSI.InteractiveDataView.currVoxel = [MRSI_x, MRSI_y,MRSI_z];
    VoxelIndices = [MRSI_x, MRSI_y,MRSI_z];
    VoxelIndices(:,1:2) = VoxelIndices(:,1:2)-1;
    
    for vox = 1 : size(VoxelIndices,1)
        ind = (VoxelIndices(vox,2)*MRSCont.raw{1}.nXvoxels+VoxelIndices(vox,1)) * 8;
        tri = delaunay(double(vertices_voxel((1:4)+ind, 2)), double(vertices_voxel((1:4)+ind, 1)));
        trisurf(tri, vertices_voxel((1:4)+ind, 2), vertices_voxel((1:4)+ind, 1), vertices_voxel((1:4)+ind, 3),...
               'LineWidth',1,'FaceColor',[255/255 140/255 0/255],'EdgeColor', 'none', 'FaceAlpha',0.85);
    end
    
    
    
    % addGrid    
    % Find vertices near the current slice (in voxel space, z-axis)
    valid_vertices = abs(vertices_voxel(:, 3) - slice) < 1; 
    
    % Transform voxel coordinates to montage coordinates        
    plot(vertices_voxel(valid_vertices, 2) , ...
         vertices_voxel(valid_vertices, 1),'.','Color',[254/255 186/255 47/255], 'MarkerSize', 6);
    
    axis image

    % % There's a flip in the location...
    MRSI_x = MRSCont.raw{1}.nXvoxels - MRSI_x + 1;
    % MRSI_y = MRSCont.raw{1}.nYvoxels - MRSI_y + 1;

    title(['MRSI voxel (' num2str(MRSI_x) ',' num2str(MRSI_y) ',' num2str(MRSI_z) ')'],...
        'Color',[11/255 71/255 111/255], 'Interpreter', 'none');

    tile_height = size(Coreg_img, 1);
    tile_width = size(Coreg_img, 2);
    
    osp_add_orientation_labels_to_montage(1, 1, 1, ...
        tile_width, tile_height, orientation_info, display_info);

    setappdata(out,'MRSCont',MRSCont);

    cla(AllAxes(1));
    set(out,'CurrentAxes',AllAxes(1))

    switch target_module
        case 'OspreyLoad'
            tempSpec = op_takeVoxel(MRSCont.raw{1},[MRSI_x MRSI_y MRSI_z]);
            try
                tempSpec = op_takesubspec(tempSpec,Subspectrum);
            catch
            end
            if zerofill > 1
                tempSpec = op_zeropad(tempSpec,zerofill);
            end
            if GaussianLB >= 1
                tempSpec = op_filter(tempSpec,GaussianLB);
            end
            if Magnitude
                    plot(tempSpec.ppm,...
                        squeeze(abs(tempSpec.specs)), 'Color', [11/255 71/255 111/255],'LineWidth',1.5);
            else
                    plot(tempSpec.ppm,...
                        squeeze(real(tempSpec.specs)), 'Color', [11/255 71/255 111/255],'LineWidth',1.5);
            end
            set(gca,'XLim',[ppm1,ppm2],'XDir','reverse','TickDir','out','XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255]);
            xlabel('chemical shift (ppm)');
            if ~isnan(Intensity1) && ~isnan(Intensity2)
                set(gca,'YLim',[Intensity1,Intensity2])
            end
            box off
        case 'OspreyProcess'
            specNames = fieldnames(MRSCont.processed);
            tempSpec = op_takeVoxel(MRSCont.processed.((FitTarget)){1},[MRSI_x MRSI_y MRSI_z]);
            if zerofill > 1
                tempSpec = op_zeropad(tempSpec,zerofill);
            end
            if GaussianLB >= 1
                tempSpec = op_filter(tempSpec,GaussianLB);
            end
            if Magnitude
                    plot(tempSpec.ppm,...
                        squeeze(abs(tempSpec.specs)), 'Color', [11/255 71/255 111/255],'LineWidth',1.5);
            else
                    plot(tempSpec.ppm,...
                        squeeze(real(tempSpec.specs)), 'Color', [11/255 71/255 111/255],'LineWidth',1.5);
            end
            set(gca,'XLim',[ppm1,ppm2],'XDir','reverse','TickDir','out','XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255]);
            xlabel('chemical shift (ppm)');
            if ~isnan(Intensity1) && ~isnan(Intensity2)
                set(gca,'YLim',[Intensity1,Intensity2])
            end
            box off
       case 'OspreyFit'
           if isfield(MRSCont.fit,FitTarget)
                ModelMatrix = MRSCont.fit.(FitTarget);
           else
               ModelMatrix = MRSCont.fit.metab;
           end
           if ~isempty(ModelMatrix{MRSI_x, MRSI_y,MRSI_z})
                switch FitPlotFunction
                    case 'Fit1DStack'
                        ModelMatrix{MRSI_x, MRSI_y,MRSI_z}.plotFit1DStack(0)
                    case 'Fit3D'
                        ModelMatrix{MRSI_x, MRSI_y,MRSI_z}.plotFit3D(0)
                end   
           end
    end

end

 function updatePlot(source, ~,out,uiStruct)
    style = get(source,'Style');
    MRSCont = getappdata(out,'MRSCont');  % Get MRSCont from hidden container in gui class  
    Coreg_img_3D = MRSCont.opts.MRSI.InteractiveDataView.CoregImage;
    vertices_voxel = MRSCont.opts.MRSI.InteractiveDataView.vertices_voxel;
    Magnitude = uiStruct.Magnitude.Value;
    ppm1 = str2num(uiStruct.ppm1.String);
    ppm2 = str2num(uiStruct.ppm2.String);
    Intensity1 = str2num(uiStruct.Intensity1.String);
    Intensity2 = str2num(uiStruct.Intensity2.String);
    zerofill = str2num(uiStruct.Zerofill.String);
    GaussianLB = str2num(uiStruct.LB.String);
    Subspectrum = str2num(uiStruct.SubSpec.String);
    target_module = MRSCont.opts.MRSI.InteractiveDataView.target_module;    
    FitTarget= MRSCont.opts.MRSI.InteractiveDataView.FitTarget;
    FitPlotFunction = MRSCont.opts.MRSI.InteractiveDataView.FitPlotFunction;
    kk = MRSCont.opts.MRSI.InteractiveDataView.kk;
    orientation_info = MRSCont.opts.MRSI.InteractiveDataView.orientation_info;
    display_info = MRSCont.opts.MRSI.InteractiveDataView.display_info;


    switch style
        case 'popupmenu'
            target_module = MRSCont.opts.MRSI.InteractiveDataView.target_module_labels{source.Value};
            MRSCont.opts.MRSI.InteractiveDataView.target_module = target_module;
        otherwise
            if ~strcmp(source.String,'Update')
                if strcmp(source.String,'>')
                    uiStruct.zVoxelInd.String = num2str(str2num(uiStruct.zVoxelInd.String) + 1);
                    if  str2num(uiStruct.zVoxelInd.String) > MRSCont.raw{1}.nZvoxels
                        uiStruct.zVoxelInd.String = num2str(MRSCont.raw{1}.nZvoxels);
                    end
                end
                if strcmp(source.String,'<')
                    uiStruct.zVoxelInd.String = num2str(str2num(uiStruct.zVoxelInd.String) - 1);
                    if  str2num(uiStruct.zVoxelInd.String) < 1
                        uiStruct.zVoxelInd.String = '1';
                    end
                end
            end           
    end

    % Now we need to read in all the inputs
    slice = str2num(uiStruct.zVoxelInd.String);
    AllAxes = out.Children(end).Children;
    cla(AllAxes(2));
    set(out,'CurrentAxes',AllAxes(2))
    max_val = prctile(Coreg_img_3D(:), 99.5);

    % Prepare image 
    Coreg_img = squeeze(Coreg_img_3D(:,:,slice));
    

    imagesc(Coreg_img,[0 max_val])
    colormap gray;
    axis image
    
    hold on;
    
    MRSI_x =  MRSCont.opts.MRSI.InteractiveDataView.currVoxel(1);
    MRSI_y =  MRSCont.opts.MRSI.InteractiveDataView.currVoxel(2);
    MRSI_z =  slice;
    MRSCont.opts.MRSI.InteractiveDataView.currVoxel(3) = MRSI_z;
    VoxelIndices = MRSCont.opts.MRSI.InteractiveDataView.currVoxel;
    VoxelIndices(:,1:2) = VoxelIndices(:,1:2)-1;
    
    for vox = 1 : size(VoxelIndices,1)
        ind = (VoxelIndices(vox,2)*MRSCont.raw{1}.nXvoxels+VoxelIndices(vox,1)) * 8;
        tri = delaunay(double(vertices_voxel((1:4)+ind, 1)), double(vertices_voxel((1:4)+ind, 2)));
        trisurf(tri, vertices_voxel((1:4)+ind, 1), vertices_voxel((1:4)+ind, 2), vertices_voxel((1:4)+ind, 3),...
               'LineWidth',1,'FaceColor',[255/255 140/255 0/255],'EdgeColor', 'none', 'FaceAlpha',0.85);
    end
        
    
    
    % addGrid    
    % Find vertices near the current slice (in voxel space, z-axis)
    valid_vertices = abs(vertices_voxel(:, 3) - slice) < 1; 
    
    % Transform voxel coordinates to montage coordinates        
    plot(vertices_voxel(valid_vertices, 1) , ...
         vertices_voxel(valid_vertices, 2),'.','Color',[254/255 186/255 47/255], 'MarkerSize', 4);
    
    % % There's a flip in the location...
    MRSI_x = MRSCont.raw{1}.nXvoxels - MRSI_x + 1;
    % MRSI_y = MRSCont.raw{1}.nYvoxels - MRSI_y + 1;

    title(['MRSI voxel (' num2str(MRSI_x) ',' num2str(MRSI_y) ',' num2str(MRSI_z) ')'],...
          'Color',[11/255 71/255 111/255], 'Interpreter', 'none');

      tile_height = size(Coreg_img, 1);
    tile_width = size(Coreg_img, 2);
    
    osp_add_orientation_labels_to_montage(1, 1, 1, ...
        tile_width, tile_height, orientation_info, display_info);

    axis image
    setappdata(out,'MRSCont',MRSCont);

    cla(AllAxes(1));
    set(out,'CurrentAxes',AllAxes(1))
    switch target_module
        case 'OspreyLoad'
            tempSpec = op_takeVoxel(MRSCont.raw{1},[MRSI_x MRSI_y MRSI_z]);
            try
                tempSpec = op_takesubspec(tempSpec,Subspectrum);
            catch
            end
            if zerofill > 1
                tempSpec = op_zeropad(tempSpec,zerofill);
            end
            if GaussianLB >= 1
                tempSpec = op_filter(tempSpec,GaussianLB);
            end
            if Magnitude
                    plot(tempSpec.ppm,...
                        squeeze(abs(tempSpec.specs)), 'Color', [11/255 71/255 111/255],'LineWidth',1.5);
            else
                    plot(tempSpec.ppm,...
                        squeeze(real(tempSpec.specs)), 'Color', [11/255 71/255 111/255],'LineWidth',1.5);
            end
            set(gca,'XLim',[ppm1,ppm2],'XDir','reverse','TickDir','out','XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255]);
            xlabel('chemical shift (ppm)');
            if ~isnan(Intensity1) && ~isnan(Intensity2)
                set(gca,'YLim',[Intensity1,Intensity2])
            end
            box off
        case 'OspreyProcess'
            specNames = fieldnames(MRSCont.processed);
            tempSpec = op_takeVoxel(MRSCont.processed.((FitTarget)){1},[MRSI_x MRSI_y MRSI_z]);
            if zerofill > 1
                tempSpec = op_zeropad(tempSpec,zerofill);
            end
            if GaussianLB >= 1
                tempSpec = op_filter(tempSpec,GaussianLB);
            end
            if Magnitude
                    plot(tempSpec.ppm,...
                        squeeze(abs(tempSpec.specs)), 'Color', [11/255 71/255 111/255],'LineWidth',1.5);
            else
                    plot(tempSpec.ppm,...
                        squeeze(real(tempSpec.specs)), 'Color', [11/255 71/255 111/255],'LineWidth',1.5);
            end
            set(gca,'XLim',[ppm1,ppm2],'XDir','reverse','TickDir','out','XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255]);
            xlabel('chemical shift (ppm)');
            if ~isnan(Intensity1) && ~isnan(Intensity2)
                set(gca,'YLim',[Intensity1,Intensity2])
            end
            box off
        case 'OspreyFit'
           if isfield(MRSCont.fit,FitTarget)
                ModelMatrix = MRSCont.fit.(FitTarget);
           else
               ModelMatrix = MRSCont.fit.metab;
           end
            switch FitPlotFunction
                case 'Fit1DStack'
                    ModelMatrix{MRSI_x, MRSI_y,MRSI_z}.plotFit1DStack(0)
                case 'Fit3D'
                    ModelMatrix{MRSI_x, MRSI_y,MRSI_z}.plotFit3D(0)
            end    
    end
    
 end