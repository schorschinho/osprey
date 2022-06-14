function out = osp_plotSegmentPET(MRSCont, kk)
%% out = osp_plotSegmentPET(MRSCont, kk)
%   Creates a figure showing segmentation of the PET image and the 
%   MRS voxel stored in an Osprey data container
%
%   USAGE:
%       out = osp_plotSegmentPET(MRSCont, kk)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2019-11-26)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2019-11-26: First version of the code.

% Check that OspreyCoreg has been run before
if ~MRSCont.flags.didCoreg
    error('Trying to plot segmentation data, but data has not been processed yet. Run OspreyCoreg first.')
end

%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
if nargin < 2
    kk = 1;
    if nargin<1
        error('ERROR: no input Osprey container specified.  Aborting!!');
    end
end

%%% 2. LOAD DATA TO PLOT %%%
% Load T1 image, mask volume, T1 max value, and voxel center
[~,filename_voxel,fileext_voxel]   = fileparts(MRSCont.files{kk});
[~,filename_image,fileext_image]   = fileparts(MRSCont.coreg.vol_image_pet{kk}.fname);

Vimage=spm_vol(MRSCont.coreg.vol_image_pet{kk}.fname);
Vmask=spm_vol(MRSCont.coreg.vol_mask_2nd{kk}.fname);
voxel_ctr = MRSCont.coreg.voxel_ctr{kk};

%%% 3. SET UP THREE PLANE IMAGE %%%
% Generate three plane image for the output
% Transform structural image and co-registered voxel mask from voxel to
% world space for output (MM: 180221)
[img_t,img_c,img_s] = voxel2world_space(Vimage,voxel_ctr);
[mask_t,mask_c,mask_s] = voxel2world_space(Vmask,voxel_ctr);

img_t = flipud(img_t/MRSCont.coreg.PETmax{kk});
img_c = flipud(img_c/MRSCont.coreg.PETmax{kk});
img_s = flipud(img_s/MRSCont.coreg.PETmax{kk});

img_t = img_t + 0.225*flipud(mask_t);
img_c = img_c + 0.225*flipud(mask_c);
img_s = img_s + 0.225*flipud(mask_s);

size_max = max([max(size(img_t)) max(size(img_c)) max(size(img_s))]);
three_plane_img = zeros([size_max 3*size_max]);
three_plane_img(:,1:size_max)              = image_center(img_t, size_max);
three_plane_img(:,size_max+(1:size_max))   = image_center(img_s, size_max);
three_plane_img(:,size_max*2+(1:size_max)) = image_center(img_c, size_max);

%%% 4. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
if ~MRSCont.flags.isGUI
    out = figure;
    title(['Coregistration with PET: ' filename_voxel fileext_voxel ' & '  filename_image fileext_image], 'Interpreter', 'none','FontSize', 16);
    set(gcf, 'Color', 'w');
else
    out = figure('Visible','off');
    title(['Coregistration with PET: ' filename_voxel fileext_voxel ' & '  filename_image fileext_image], 'Interpreter', 'none','FontSize', 16,'Color', MRSCont.colormap.Foreground);
end

axesHandles.coregImage = subplot(2, 3, [1,2,3]);
imagesc(axesHandles.coregImage, three_plane_img);
colormap('gray');
caxis([0 1])
axis equal;
axis tight;
axis off;
if ~MRSCont.flags.isGUI
    title(['Coregistration with PET: ' filename_voxel fileext_voxel ' & '  filename_image fileext_image], 'Interpreter', 'none','FontSize', 16);
else
    title(['Coregistration with PET: ' filename_voxel fileext_voxel ' & '  filename_image fileext_image], 'Interpreter', 'none','FontSize', 16,'Color', MRSCont.colormap.Foreground);
end

%%% 5. ADD PET IMAGE INTENSITY DISTRIBUTION PLOT %%%

%First, plot overall distribution
axesHandles.IntensityDistributionPlot{1} = subplot(2, 3, 4);
xIntensity              = MRSCont.coreg.pet.histogram.xIntensity{kk};
yIntensity              = MRSCont.coreg.pet.histogram.yIntensity{kk};
GaussModelParams        = MRSCont.coreg.pet.histogram.fitParams{kk};
yFit = GaussModel(GaussModelParams,xIntensity);
plot(xIntensity, yIntensity, 'LineWidth', 1, 'Color', MRSCont.colormap.Foreground, 'Parent', axesHandles.IntensityDistributionPlot{1});
hold on;
plot(xIntensity, yFit, 'LineWidth', 1.6, 'Color', MRSCont.colormap.Accent, 'Parent', axesHandles.IntensityDistributionPlot{1});
hold off;
if ~MRSCont.flags.isGUI
    title('Full voxel PET intensities and Gaussian fit', 'Interpreter', 'none','FontSize', 16);
else
    title('Full voxel PET intensities and Gaussian fit', 'Interpreter', 'none','FontSize', 16, 'Color', MRSCont.colormap.Foreground);
end
set(axesHandles.IntensityDistributionPlot{1}, 'LineWidth', 1, 'TickDir', 'out');
set(axesHandles.IntensityDistributionPlot{1}, 'FontSize', 16);
set(axesHandles.IntensityDistributionPlot{1}, 'Units', 'normalized');
box off;
xlab = 'Full voxel PET Image Intensity (a.u.)';
ylab = 'Number of pixels';
%legend(axesHandles.IntensityDistributionPlot, 'Image data', 'Gaussian fit');
xlabel(axesHandles.IntensityDistributionPlot{1}, xlab, 'FontSize', 16);
ylabel(axesHandles.IntensityDistributionPlot{1}, ylab, 'FontSize', 16);

% Next plot GM/WM distributions
rr=1; %CWDJ Assume 1st voxel
if MRSCont.flags.isPRIAM
    warning('Assuming 1st voxel')
end
tissueTypes = {'GM', 'WM'};
for tt = 1:length(tissueTypes)
    % Navigate to the appropriate subplot tile
    axesHandles.IntensityDistributionPlot{tt+1} = subplot(2, 3, 4+tt);
    % Extract the parameters from the MRSCont
    rawPETIntensitySum      = MRSCont.seg.pet.rawPETIntensitySum.(tissueTypes{tt})(kk);
    xIntensity              = MRSCont.seg.pet.histogram.xIntensity.(tissueTypes{tt}){kk}{rr};
    yIntensity              = MRSCont.seg.pet.histogram.yIntensity.(tissueTypes{tt}){kk}{rr};
    GaussModelParams        = MRSCont.seg.pet.histogram.fitParams.(tissueTypes{tt}){kk}{rr};
    mostFrequentIntensity   = MRSCont.seg.pet.histogram.mostFrequentIntensity.(tissueTypes{tt})(kk);
    distFWHM                = MRSCont.seg.pet.histogram.distFWHM.(tissueTypes{tt})(kk);
    % Evaluate model
    yFit = GaussModel(GaussModelParams,xIntensity);
    plot(xIntensity, yIntensity, 'LineWidth', 1, 'Color', MRSCont.colormap.Foreground, 'Parent', axesHandles.IntensityDistributionPlot{tt+1});
    hold on;
    plot(xIntensity, yFit, 'LineWidth', 1.6, 'Color', MRSCont.colormap.Accent, 'Parent', axesHandles.IntensityDistributionPlot{tt+1});
    hold off;
    if ~MRSCont.flags.isGUI
        title([tissueTypes{tt} ' PET intensities and Gaussian fit'], 'Interpreter', 'none','FontSize', 16);
    else
        title([tissueTypes{tt} ' PET intensities and Gaussian fit'], 'Interpreter', 'none','FontSize', 16, 'Color', MRSCont.colormap.Foreground);
    end
    set(axesHandles.IntensityDistributionPlot{tt+1}, 'LineWidth', 1, 'TickDir', 'out');
    set(axesHandles.IntensityDistributionPlot{tt+1}, 'FontSize', 16);
    set(axesHandles.IntensityDistributionPlot{tt+1}, 'Units', 'normalized');
    box off;
    xlab = [tissueTypes{tt} ' PET Image Intensity (a.u.)'];
    ylab = 'Number of pixels';
    %legend(axesHandles.IntensityDistributionPlot, 'Image data', 'Gaussian fit');
    xlabel(axesHandles.IntensityDistributionPlot{tt+1}, xlab, 'FontSize', 16);
    ylabel(axesHandles.IntensityDistributionPlot{tt+1}, ylab, 'FontSize', 16);
end

%%% 6. ADD OSPREY LOGO %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end
end

   