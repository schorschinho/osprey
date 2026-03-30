function out = osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices, target_image,target_module,FitTarget,FitPlotFunction,addGrid,zerofill,GaussianLB,Magnitude,fillOut,convention)
%% out = osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices, target_image,target_module,FitTarget,FitPlotFunction,addGrid,zerofill,GaussianLB,Magnitude,fillOut,convention)
%   Creates a figure showing the MRSI spectra and voxel locations
%
%   USAGE:
%       out = osp_plotSpecAndLocMRSI(MRSCont,VoxelIndices, target_image,target_module,FitTarget,FitPlotFunction,addGrid,zerofill,GaussianLB,Magnitude,fillOut,convention)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   INPUTS:
%       MRSCont  = Osprey data container.
%       VoxelIndices = Vector with indices; three per column
%       target_image = Target image to overlay on
%       target_module = Target module to plot (Load, Process, Fit)
%       FitTarget = Target to plot for fits or an index for subspectra
%       FitPlotFunction = Plot function to use for fits (Fit1DStack, Fit3D)
%       addGrid  = add Grid dots on image
%       zerofill = zero-fill factor to add
%       GaussianLB = Gaussian linebroadening to add
%       Magnitude = plot magnitude spectra
%       fillOut = fill the MRSI voxel locations
%       convention   = 'radiological' (L on right) or 'neurological' (L on left)
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-08-04)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-08-04: First version of the code.

% Fall back to defaults if not provided
if nargin < 12
    convention = 'neurological';
    if nargin < 11
        fillOut = 1;
        if nargin < 10
            Magnitude = 0;
            if nargin < 9
                GaussianLB = 1;
                if nargin < 8
                    zerofill = 2;
                    if nargin<7
                        addGrid = 1;
                        if nargin < 6
                            FitPlotFunction = 'Fit1DStack';
                            if nargin < 5
                                FitTarget = 'metab';
                                if nargin < 4
                                    target_module = 'OspreyLoad';
                                    if nargin < 3
                                        if isfield(MRSCont, 'files_nii_MRSIloc') && ~isempty(MRSCont.files_nii_MRSIloc{1})
                                            target_image = 'T1w_rMRSIloc'; 
                                        else
                                            target_image = 'T1w_rMRSI';
                                        end
                                        if nargin<2
                                            VoxelIndices = [round(MRSCont.raw{1}.nXvoxels/2),round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2);
                                                            round(MRSCont.raw{1}.nXvoxels/2)+1,round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2);
                                                            round(MRSCont.raw{1}.nXvoxels/2)+2,round(MRSCont.raw{1}.nYvoxels/2),round(MRSCont.raw{1}.nZvoxels/2);];
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

switch  target_module
    case 'OspreyProcess'
        if ~MRSCont.flags.didProcess
            error('Trying to plot processed spectra, but processing has not been performed yet. Run OspreyProcess first.')
        end
    case 'OspreyFit'
        if ~MRSCont.flags.didFit
            error('Trying to plot model results, but fitting has not been performed yet. Run OspreyFit first.')
        end
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
[VoxColors]=cbrewer('qual', 'Set1', 9);

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

max_val = prctile(Coreg_img_display(:), 99.5);


VoxelIndicesOverlay = VoxelIndices;
VoxelIndicesOverlay(:,1:2) = VoxelIndicesOverlay(:,1:2)-1;


NumberOfSpecs = size(VoxelIndices,1); 

out = figure;   
tiledlayout(1,NumberOfSpecs + 1,'TileSpacing','compact')
nexttile
imagesc(squeeze(Coreg_img_display(:,:,VoxelIndices(1,3))),[0 max_val])
colormap gray;
axis image
hold on
box off
set(gca,'color', [0 0 0]);
currentX = xlim;
currentY = ylim;
buffer = 10;
xlim([currentX(1) - buffer, currentX(2) + buffer]);
ylim([currentY(1) - buffer, currentY(2) + buffer]);

for vox = 1 : NumberOfSpecs
    ind = (VoxelIndicesOverlay(vox,2)*MRSCont.raw{1}.nXvoxels+VoxelIndicesOverlay(vox,1)) * 8;
    if fillOut
        tri = delaunay(double(vertices_voxel((1:4)+ind, 1)), double(vertices_voxel((1:4)+ind, 2)));
        trisurf(tri, vertices_voxel((1:4)+ind, 1), vertices_voxel((1:4)+ind, 2), vertices_voxel((1:4)+ind, 3),...
               'LineWidth',1,'FaceColor',VoxColors(vox,:),'EdgeColor', 'none', 'FaceAlpha',0.95);
    else
        plot(vertices_voxel((1:4)+ind, 1), vertices_voxel((1:4)+ind, 2),'.','Color',VoxColors(vox,:), 'MarkerSize', 12)
    end
end

if addGrid
    valid_vertices = abs(vertices_voxel(:, 3) - VoxelIndices(1,3)) < 1; 

      plot(vertices_voxel(valid_vertices, 1) , ...
     vertices_voxel(valid_vertices, 2),'.','Color',[254/255 186/255 47/255], 'MarkerSize', 4);    
end

ax = gca;
ax.XTickLabel = {};
ax.YTickLabel = {};

if strcmp(target_module,'OspreyFit')
    if (MRSCont.raw{1}.nZvoxels > 1)
        ModelMatrix = flip(MRSCont.fit.(FitTarget),length(size(MRSCont.fit.(FitTarget))));
    else
        ModelMatrix = MRSCont.fit.(FitTarget);
    end
    ModelMatrix = MRSCont.fit.(FitTarget);
end

for vox = 1 : NumberOfSpecs
    nexttile
    switch  target_module
        case 'OspreyLoad'
            Subspectrum = FitTarget;
            tempSpec = op_takeVoxel(MRSCont.raw{1},[VoxelIndices(vox,1) VoxelIndices(vox,2) VoxelIndices(vox,3)]);
            if tempSpec.dims.subSpecs > 0
                tempSpec = op_takesubspec(tempSpec,Subspectrum);
            end
            if zerofill > 1
                tempSpec = op_zeropad(tempSpec,zerofill);
            end
            if GaussianLB >= 1
                tempSpec = op_filter(tempSpec,GaussianLB);
            end
            if Magnitude
                    plot(tempSpec.ppm,...
                        squeeze(abs(tempSpec.specs)), 'Color', VoxColors(vox,:),'LineWidth',1.5);
            else
                    plot(tempSpec.ppm,...
                        squeeze(real(tempSpec.specs)), 'Color', VoxColors(vox,:),'LineWidth',1.5);
            end
            set(gca,'XLim',[0.5,6],'XDir','reverse','TickDir','out','XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255]);
            xlabel('chemical shift (ppm)');
            box off
        case 'OspreyProcess'
            Subspectrum = FitTarget;
            specNames = fieldnames(MRSCont.processed);
            tempSpec = op_takeVoxel(MRSCont.processed.(Subspectrum){1},[VoxelIndices(vox,1) VoxelIndices(vox,2) VoxelIndices(vox,3)]);
            if zerofill > 1
                tempSpec = op_zeropad(tempSpec,zerofill);
            end
            if GaussianLB >= 1
                tempSpec = op_filter(tempSpec,GaussianLB);
            end
            % tempSpec = op_autophase(tempSpec,1.9,2.1);
            if Magnitude
                    plot(tempSpec.ppm,...
                        squeeze(abs(tempSpec.specs)), 'Color', VoxColors(vox,:),'LineWidth',1.5);
            else
                    plot(tempSpec.ppm,...
                        squeeze(real(tempSpec.specs)), 'Color', VoxColors(vox,:),'LineWidth',1.5);
            end
            if ~strcmp(FitTarget,'w' )
                set(gca,'XLim',[0.5,4],'XDir','reverse','TickDir','out','XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255]);
            else
                set(gca,'XLim',[4.65-3,4.65+3],'XDir','reverse','TickDir','out','XColor', [11/255 71/255 111/255],'YColor', [11/255 71/255 111/255]);
            end
            xlabel('chemical shift (ppm)');
            box off
        case 'OspreyFit'
            switch FitPlotFunction
                case 'Fit1DStack'
                    ModelMatrix{VoxelIndices(vox,1),VoxelIndices(vox,2),VoxelIndices(vox,3)}.plotFit1DStack(0)
                case 'Fit3D'
                    ModelMatrix{VoxelIndices(vox,1),VoxelIndices(vox,2),VoxelIndices(vox,3)}.plotFit3D(0)
            end
            children = allchild(gca);
            for ch = 1 : 4
                children(ch).Color = VoxColors(vox,:);
            end
            for ch = 6 :2 :  length(children)
                children(ch).Color = VoxColors(vox,:);
            end
            if vox ~= NumberOfSpecs
                for ch = 5 :2 :  length(children)
                    delete(children(ch));
                end
            end
    end
    title(['(' num2str(VoxelIndices(vox,1)) ', ' num2str(VoxelIndices(vox,2)) ', ' num2str(VoxelIndices(vox,3)) ')' ],'Color',[110/255 136/255 164/255])
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

   