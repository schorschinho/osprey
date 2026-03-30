function out = osp_plotInteractiveMRSIMask(MRSCont, target_image, convention)
%% out = osp_plotInteractiveMRSIMask(MRSCont, target_image, convention)
%   Creates a figure showing the quick integration maps of an MRSI scan
%
%   USAGE:
%       out = osp_plotInteractiveMRSIMask(MRSCont, target_image, convention)
%
%   INPUTS:
%       MRSCont  = Osprey data container.
%       target_image = Target image to overlay on
%       convention   = 'radiological' (L on right) or 'neurological' (L on left)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.
%% Fall back to defaults if not provided
if nargin < 3
    convention = 'neurological';
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


%% Validate prerequisites
if ~MRSCont.flags.didLoadData
    error('Trying to plot coregistration, but data has not been loaded yet. Run OspreyLoad first.')
end

if ~MRSCont.flags.didCoreg
    error('Trying to plot coregistration, but coregistration has not been performed yet. Run OspreyCoreg first.')
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

% Also transform for MRSI voxel space (needed for slice selection)
vertices_homogeneous = [vertices, ones(size(vertices, 1), 1)];
vertices_MRSI_voxel = (MRSI_vol.mat \ vertices_homogeneous')';
vertices_MRSI_voxel = vertices_MRSI_voxel(:, 1:3);

%% Setup figure
out = figure;  

set(out, 'units','normalized','outerposition',[0 0 0.75 1]);
set(out, 'units','pixel');

set(out, 'Color', [1 1 1]);


max_val = prctile(Coreg_img_display(:), 99.5);
slice = round(MRSCont.raw{1}.nZvoxels/2);
imagesc(flip(rot90(squeeze(Coreg_img_display(:,:,slice))),2),[0 max_val])
colormap gray;
axis image
axis off
box off
set(gca,'color', [0 0 0]);
currentX = xlim;
currentY = ylim;
buffer = 10;
xlim([currentX(1) - buffer, currentX(2) + buffer]);
ylim([currentY(1) - buffer, currentY(2) + buffer]);


hold on;


MRSCont.opts.MRSI.outerMask.mask = zeros(MRSCont.raw{1}.nXvoxels,MRSCont.raw{1}.nYvoxels,MRSCont.raw{1}.nZvoxels);
MRSCont.opts.MRSI.outerMask.mask(MRSCont.opts.MRSI.outerMask.x(1):MRSCont.opts.MRSI.outerMask.x(2),...
MRSCont.opts.MRSI.outerMask.y(1):MRSCont.opts.MRSI.outerMask.y(2),...
MRSCont.opts.MRSI.outerMask.z(1):MRSCont.opts.MRSI.outerMask.z(2)) = 1;
MRSCont.opts.MRSI.outerMask.CoregImage = Coreg_img_display;
MRSCont.opts.MRSI.outerMask.vertices_voxel = vertices_voxel;
MRSCont.opts.MRSI.outerMask.kk = 1;
MRSCont.opts.MRSI.outerMask.saved = 0;

setappdata(out,'MRSCont',MRSCont);

% addGrid    
% Find vertices near the current slice (in voxel space, z-axis)
valid_vertices = abs(vertices_display(:, 3) - slice) < 1; 

% Transform voxel coordinates to montage coordinates        
plot(vertices_display(valid_vertices, 2) , ...
     vertices_display(valid_vertices, 1),'.','Color',[254/255 186/255 47/255], 'MarkerSize', 6);

% osp_plot_grid_overlay(vertices_display, vertices_MRSI_voxel, ...
%         slice, slice, 1, tile_width, tile_height, ...
%         orientation_info, display_info);


axis image


uiStruct.LeftButton = uicontrol(  'Style', 'pushbutton', ...
              'String', '<', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)/2 - 70, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.zVoxelInd = uicontrol(  'Style', 'Edit', ...
              'String', round(MRSCont.raw{1}.nZvoxels/2), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)/2, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.RightButton = uicontrol(  'Style', 'pushbutton', ...
              'String', '>', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)/2 + 70, out.InnerPosition(4)*0.03, 60, 30]);

uiStruct.Textouter = uicontrol(  'Style', 'Text', ...
              'String', 'MRSI slice index', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)/2, out.InnerPosition(4)*0.07, 150, 20]);

uiStruct.Textouter = uicontrol(  'Style', 'Text', ...
              'String', 'Outer MRSI Mask Indices', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*0.85, out.InnerPosition(4)*0.85, 150, 20]);

uiStruct.Textouter = uicontrol(  'Style', 'Text', ...
              'String', 'Index 1', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*0.85, out.InnerPosition(4)*0.825, 100, 20]);

uiStruct.Textouter = uicontrol(  'Style', 'Text', ...
              'String', 'Index 2', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*0.9, out.InnerPosition(4)*0.825, 100, 20]);

uiStruct.outerMaskX1 = uicontrol(  'Style', 'Edit', ...
              'String', MRSCont.opts.MRSI.outerMask.x(1), ...
               'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*0.85, out.InnerPosition(4)*0.8, 40, 20]);

uiStruct.outerMaskX2 = uicontrol(  'Style', 'Edit', ...
              'String', MRSCont.opts.MRSI.outerMask.x(2), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*0.9, out.InnerPosition(4)*0.8, 40, 20]);

uiStruct.TextX = uicontrol(  'Style', 'Text', ...
              'String', ['X max(' num2str(MRSCont.raw{1}.nXvoxels) ')'], ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*0.95, out.InnerPosition(4)*0.8, 100, 20]);

uiStruct.outerMaskY1 = uicontrol(  'Style', 'Edit', ...
              'String', MRSCont.opts.MRSI.outerMask.y(1), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*0.85, out.InnerPosition(4)*0.75, 40, 20]);

uiStruct.outerMaskY2 = uicontrol(  'Style', 'Edit', ...
              'String', MRSCont.opts.MRSI.outerMask.y(2), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*0.9, out.InnerPosition(4)*0.75, 40, 20]);

uiStruct.TextY = uicontrol(  'Style', 'Text', ...
              'String', ['Y max(' num2str(MRSCont.raw{1}.nYvoxels) ')'], ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*0.95, out.InnerPosition(4)*0.75, 100, 20]);


uiStruct.outerMaskZ1 = uicontrol(  'Style', 'Edit', ...
              'String', MRSCont.opts.MRSI.outerMask.z(1), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*0.85, out.InnerPosition(4)*0.7, 40, 20]);

uiStruct.outerMaskZ2 = uicontrol(  'Style', 'Edit', ...
              'String', MRSCont.opts.MRSI.outerMask.z(2), ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*0.9, out.InnerPosition(4)*0.7, 40, 20]);

uiStruct.TextZ = uicontrol(  'Style', 'Text', ...
              'String', ['Z max(' num2str(MRSCont.raw{1}.nZvoxels) ')'], ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*0.95, out.InnerPosition(4)*0.7, 100, 20]);


uiStruct.UpdateButton = uicontrol(  'Style', 'pushbutton', ...
              'String', 'Update', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*0.85, out.InnerPosition(4)*0.65, 60, 30]);

uiStruct.SaveButton = uicontrol(  'Style', 'pushbutton', ...
              'String', 'Save', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*0.90, out.InnerPosition(4)*0.65, 60, 30]);

uiStruct.TextNeigh = uicontrol(  'Style', 'Text', ...
              'String', 'Neighborhood Size', ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'HorizontalAlignment','left',...
              'Position', [out.InnerPosition(3)*0.85, out.InnerPosition(4)*0.60, 150, 20]);

uiStruct.Neigh = uicontrol(  'Style', 'Edit', ...
              'String', 1, ...
              'ForegroundColor',[11/255 71/255 111/255],...
              'BackgroundColor',[1 1 1],...
              'Position', [out.InnerPosition(3)*0.85, out.InnerPosition(4)*0.58, 40, 20]);

uiStruct.AddRem = uicontrol('Style', 'checkbox', ...
                           'String', 'Add/Remove Mask', ...
                           'ForegroundColor',[11/255 71/255 111/255],...
                           'BackgroundColor',[1 1 1],...
                           'Position', [out.InnerPosition(3)*0.85, out.InnerPosition(4)*0.55, 150, 30], ...
                           'Value', 1);

set(uiStruct.LeftButton,'Callback', {@updatePlot,out,uiStruct})
set(uiStruct.RightButton,'Callback', {@updatePlot,out,uiStruct})
set(uiStruct.UpdateButton,'Callback', {@updatePlot,out,uiStruct})
set(uiStruct.SaveButton,'Callback', {@SaveMRSCont,out})

% Set up the click callback function
set(gcf, 'WindowButtonDownFcn',{@mouseClick,out,uiStruct});

% Set up the close callback function
set(gcf, 'CloseRequestFcn', @CloseFunction);


%% Cleanup compressed files
if strcmp(T1ext, '.gz')
    delete(MRSCont.files_nii{1});
    MRSCont.files_nii{1} = strrep(MRSCont.files_nii{1}, '.nii', '.nii.gz');
end

if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
    if strcmp(MRSIlocext, '.gz')
        delete(MRSCont.files_nii_MRSIloc{1});
        MRSCont.files_nii_MRSIloc{1} = strrep(MRSCont.files_nii_MRSIloc{1}, '.nii', '.nii.gz');
    end
end
end

% Callback function for mouse clicks
function mouseClick(~, ~,out,uiStruct)
    MRSCont = getappdata(out,'MRSCont');  % Get MRSCont from hidden container in gui class 
    Coreg_img = MRSCont.opts.MRSI.outerMask.CoregImage;
    AddRem = uiStruct.AddRem.Value;
    Neigh = str2num(uiStruct.Neigh.String);

    % Get the current point
    currentPoint = get(gca, 'CurrentPoint');
    x = round(currentPoint(1, 1));
    y = round(currentPoint(1, 2));
    slice_index = str2num(uiStruct.zVoxelInd.String);

    vertices_voxel = MRSCont.opts.MRSI.outerMask.vertices_voxel;
    valid_vertices = abs(vertices_voxel(:, 3) - slice_index) < 1; 
    vertices_voxel_click = vertices_voxel;
    vertices_voxel_click(~valid_vertices,3)=MRSCont.raw{1}.nZvoxels*10;
    distances = sqrt((vertices_voxel_click(:,1) - x).^2 + (vertices_voxel_click(:,2) - y).^2 + (vertices_voxel_click(:,3) - slice_index).^2 );
    distances = sum(reshape(distances,8,[]),1);
    [~, row_index] = min(distances);
    [MRSI_x,MRSI_y,~]=ind2sub(size(MRSCont.opts.MRSI.outerMask.mask),row_index);
    MRSI_z = slice_index;
    
    % Display click coordinates
    % fprintf('Clicked at: (%.2f, %.2f)\n', x, y);
    % fprintf('MRSI voxel: (%.2f, %.2f, %.2f)\n', MRSI_x, MRSI_y, MRSI_z);

    
    if AddRem
        MRSCont.opts.MRSI.outerMask.mask(MRSI_x-Neigh+1:MRSI_x,MRSI_y-Neigh+1:MRSI_y,MRSI_z) = 0;
    else
        MRSCont.opts.MRSI.outerMask.mask(MRSI_x-Neigh+1:MRSI_x,MRSI_y-Neigh+1:MRSI_y,MRSI_z) = 1;
    end

    cla(gca);
    max_val = prctile(Coreg_img(:), 99.5);

    imagesc(flip(rot90(squeeze(Coreg_img(:,:,slice_index))),2),[0 max_val])
    colormap gray;
    axis image
    
    hold on;
    
    [VoxColors]=cbrewer('qual', 'Set1', 9);
    
    
    linear_indices = find(MRSCont.opts.MRSI.outerMask.mask == 0);
    [xVoxelIndices, yVoxelIndices, zVoxelIndices] = ind2sub(size(MRSCont.opts.MRSI.outerMask.mask), linear_indices);
    VoxelIndices = [xVoxelIndices, yVoxelIndices,zVoxelIndices];
    VoxelIndices(:,1:2) = VoxelIndices(:,1:2)-1;
    
    VoxelIndices(zVoxelIndices ~= slice_index,:) = [];
    
    for vox = 1 : size(VoxelIndices,1)
        ind = (VoxelIndices(vox,2)*MRSCont.raw{1, 1}.nXvoxels+VoxelIndices(vox,1)) * 8;
        tri = delaunay(double(vertices_voxel((1:4)+ind, 1)), double(vertices_voxel((1:4)+ind, 2)));
        trisurf(tri, vertices_voxel((1:4)+ind, 1), vertices_voxel((1:4)+ind, 2), vertices_voxel((1:4)+ind, 3),...
               'LineWidth',1,'FaceColor',VoxColors(1,:),'EdgeColor', 'none', 'FaceAlpha',0.75);
    end
    
    
    
    % addGrid    
    % Find vertices near the current slice (in voxel space, z-axis)
    valid_vertices = abs(vertices_voxel(:, 3) - slice_index) < 1; 
    
    % Transform voxel coordinates to montage coordinates        
    plot(vertices_voxel(valid_vertices, 1) , ...
         vertices_voxel(valid_vertices, 2),'.','Color',[254/255 186/255 47/255], 'MarkerSize', 4);
    
    axis image

    % Make sure that update triggers the not saved flag
    MRSCont.opts.MRSI.outerMask.saved = 0;

    setappdata(out,'MRSCont',MRSCont);

end

 function SaveMRSCont(~,~,out)
    MRSCont = getappdata(out,'MRSCont');  % Get MRSCont from hidden container in gui class 
    outputFolder    = MRSCont.outputFolder;
    outputFile      = MRSCont.outputFile;
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end
    
    save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
    fprintf('Saved Mask to MRSCont: %s\n', fullfile(outputFolder, outputFile));

    % Mask has been saved so set flag
    MRSCont.opts.MRSI.outerMask.saved = 1;
 end

 function CloseFunction(out,~)
    MRSCont = getappdata(out,'MRSCont');  % Get MRSCont from hidden container in gui class 
    if ~MRSCont.opts.MRSI.outerMask.saved
        outputFolder    = MRSCont.outputFolder;
        outputFile      = MRSCont.outputFile;
        if ~exist(outputFolder,'dir')
            mkdir(outputFolder);
        end
        
        save(fullfile(outputFolder, outputFile), 'MRSCont','-v7.3');
        fprintf('Saved Mask to MRSCont: %s\n', fullfile(outputFolder, outputFile));
    
        % Mask has been saved so set flag
        MRSCont.opts.MRSI.outerMask.saved = 1;
    end
    delete(out);
 end

 function updatePlot(source, ~,out,uiStruct)
    MRSCont = getappdata(out,'MRSCont');  % Get MRSCont from hidden container in gui class  
    Coreg_img = MRSCont.opts.MRSI.outerMask.CoregImage; 
    vertices_voxel = MRSCont.opts.MRSI.outerMask.vertices_voxel;
    kk = MRSCont.opts.MRSI.outerMask.kk;
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

    % Now we need to read in all the inputs
    slice = str2num(uiStruct.zVoxelInd.String);
    MRSCont.opts.MRSI.outerMask.x(1) = str2num(uiStruct.outerMaskX1.String);
    MRSCont.opts.MRSI.outerMask.x(2) = str2num(uiStruct.outerMaskX2.String);
    MRSCont.opts.MRSI.outerMask.y(1) = str2num(uiStruct.outerMaskY1.String);
    MRSCont.opts.MRSI.outerMask.y(2) = str2num(uiStruct.outerMaskY2.String);
    MRSCont.opts.MRSI.outerMask.z(1) = str2num(uiStruct.outerMaskZ1.String);
    MRSCont.opts.MRSI.outerMask.z(2) = str2num(uiStruct.outerMaskZ2.String);


    % Update mask
    if strcmp(source.String,'Update')
        MRSCont.opts.MRSI.outerMask.mask = zeros(MRSCont.raw{1}.nXvoxels,MRSCont.raw{1}.nYvoxels,MRSCont.raw{1}.nZvoxels);
        MRSCont.opts.MRSI.outerMask.mask(MRSCont.opts.MRSI.outerMask.x(1):MRSCont.opts.MRSI.outerMask.x(2),...
        MRSCont.opts.MRSI.outerMask.y(1):MRSCont.opts.MRSI.outerMask.y(2),...
        MRSCont.opts.MRSI.outerMask.z(1):MRSCont.opts.MRSI.outerMask.z(2)) = 1;
    end


    cla(gca);
    max_val = prctile(Coreg_img(:), 99.5);

    imagesc(flip(rot90(squeeze(Coreg_img(:,:,slice))),2),[0 max_val])
    colormap gray;
    axis image
    
    hold on;
    
    [VoxColors]=cbrewer('qual', 'Set1', 9);
    
    
    linear_indices = find(MRSCont.opts.MRSI.outerMask.mask == 0);
    [xVoxelIndices, yVoxelIndices, zVoxelIndices] = ind2sub(size(MRSCont.opts.MRSI.outerMask.mask), linear_indices);
    VoxelIndices = [xVoxelIndices, yVoxelIndices,zVoxelIndices];
    VoxelIndices(:,1:2) = VoxelIndices(:,1:2)-1;
    
    VoxelIndices(zVoxelIndices ~= slice,:) = [];
    
    for vox = 1 : size(VoxelIndices,1)
        ind = (VoxelIndices(vox,2)*MRSCont.raw{1, 1}.nXvoxels+VoxelIndices(vox,1)) * 8;
        tri = delaunay(double(vertices_voxel((1:4)+ind, 1)), double(vertices_voxel((1:4)+ind, 2)));
        trisurf(tri, vertices_voxel((1:4)+ind, 1), vertices_voxel((1:4)+ind, 2), vertices_voxel((1:4)+ind, 3),...
               'LineWidth',1,'FaceColor',VoxColors(1,:),'EdgeColor', 'none', 'FaceAlpha',0.75);
    end
    
    
    
    % addGrid    
    % Find vertices near the current slice (in voxel space, z-axis)
    valid_vertices = abs(vertices_voxel(:, 3) - slice) < 1; 
    
    % Transform voxel coordinates to montage coordinates        
    plot(vertices_voxel(valid_vertices, 1) , ...
         vertices_voxel(valid_vertices, 2),'.','Color',[254/255 186/255 47/255], 'MarkerSize', 4);
    
    axis image

    % Make sure that update triggers the not saved flag
    MRSCont.opts.MRSI.outerMask.saved = 0;
    setappdata(out,'MRSCont',MRSCont);
 end