function out = osp_plotCoregPET(MRSCont, kk)
%% out = osp_plotCoregPET(MRSCont, kk)
%   Creates a figure showing coregistration between the PET image and the 
%   MRS voxel stored in an Osprey data container
%
%   USAGE:
%       out = osp_plotCoregPET(MRSCont, kk, GUI)
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
    error('Trying to plot coregistration data, but data has not been processed yet. Run OspreyCoreg first.')
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

axesHandles.coregImage = subplot(2, 1, 1);
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
% Extract the parameters from the MRSCont
axesHandles.IntensityDistributionPlot = subplot(2, 1, 2);
rawPETIntensitySum      = MRSCont.coreg.pet.metrics.rawPETIntensitySum(kk);
xIntensity              = MRSCont.coreg.pet.histogram.xIntensity{kk};
yIntensity              = MRSCont.coreg.pet.histogram.yIntensity{kk};
GaussModelParams        = MRSCont.coreg.pet.histogram.fitParams{kk};
mostFrequentIntensity   = MRSCont.coreg.pet.histogram.mostFrequentIntensity(kk);
distFWHM                = MRSCont.coreg.pet.histogram.distFWHM(kk);
% Evaluate model
yFit = GaussModel(GaussModelParams,xIntensity);
plot(xIntensity, yIntensity, 'LineWidth', 1, 'Color', MRSCont.colormap.Foreground, 'Parent', axesHandles.IntensityDistributionPlot);
hold on;
plot(xIntensity, yFit, 'LineWidth', 1.6, 'Color', MRSCont.colormap.Accent, 'Parent', axesHandles.IntensityDistributionPlot);
hold off;
if ~MRSCont.flags.isGUI
    title('PET image intensity distribution and Gaussian fit', 'Interpreter', 'none','FontSize', 16);
else
    title('PET image intensity distribution and Gaussian fit', 'Interpreter', 'none','FontSize', 16, 'Color', MRSCont.colormap.Foreground);
end
set(axesHandles.IntensityDistributionPlot, 'LineWidth', 1, 'TickDir', 'out');
set(axesHandles.IntensityDistributionPlot, 'FontSize', 16);
set(axesHandles.IntensityDistributionPlot, 'Units', 'normalized');
% Determine the plot limits
lowestXValue = min(cell2mat(MRSCont.coreg.pet.histogram.xIntensity));
highestXValue = max(cell2mat(MRSCont.coreg.pet.histogram.xIntensity));
set(axesHandles.IntensityDistributionPlot, 'xLim', [lowestXValue highestXValue]);
box off;
xlab = 'PET Image Intensity (a.u.)';
ylab = 'Number of pixels';
%legend(axesHandles.IntensityDistributionPlot, 'Image data', 'Gaussian fit');
xlabel(axesHandles.IntensityDistributionPlot, xlab, 'FontSize', 16);
ylabel(axesHandles.IntensityDistributionPlot, ylab, 'FontSize', 16);

%%% 6. ADD OSPREY LOGO %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end
end

   