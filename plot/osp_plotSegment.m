function out = osp_plotSegment(MRSCont, kk, GUI)
%% out = osp_plotSegment(MRSCont, kk, GUI)
%   Creates a figure showing output from the segemntation routine
%   stored in an Osprey data container
%
%   USAGE:
%       out = osp_plotSegment(MRSCont, kk, GUI)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       GUI      = flag to decide whether plot is used in GUI
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2019-11-26)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2019-11-26: First version of the code.

% Check that OspreySeg has been run before
if ~MRSCont.flags.didSeg
    error('Trying to plot segmentation data, but data has not been processed yet. Run OspreySeg first.')
end

%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
if nargin < 3
    GUI = 0;
    if nargin < 2
        kk = 1;
        if nargin<1
            error('ERROR: no input Osprey container specified.  Aborting!!');
        end
    end
end

%%% 2. LOAD DATA TO PLOT %%%
% Load T1 image, mask volume, T1 max value, and voxel center
% Get the input file name
[path_voxel,filename_voxel,fileext_voxel]   = fileparts(MRSCont.files{kk});
[~,filename_image,fileext_image]   = fileparts(MRSCont.coreg.vol_image{kk}.fname);
% For batch analysis, get the last two sub-folders (e.g. site and
% subject)
path_split          = regexp(path_voxel,filesep,'split');
if length(path_split) > 2
    saveName = [path_split{end-1} '_' path_split{end} '_' filename_voxel];
end

segDestination = fullfile(MRSCont.outputFolder, 'SegMaps');
GM  = fullfile(segDestination, [saveName '_GM.nii']);
WM  = fullfile(segDestination, [saveName '_WM.nii']);
CSF = fullfile(segDestination, [saveName '_CSF.nii']);

vol_GM_mask  = spm_vol(GM);
vol_WM_mask  = spm_vol(WM);
vol_CSF_mask = spm_vol(CSF);

Vimage=spm_vol(MRSCont.coreg.vol_image{kk}.fname);

Vmask=spm_vol(MRSCont.coreg.vol_mask{kk}.fname);

voxel_ctr = MRSCont.coreg.voxel_ctr{kk};
%%% 3. SET UP THREE PLANE IMAGE %%%
% Generate three plane image for the output
% Transform structural image and co-registered voxel mask from voxel to
% world space for output (MM: 180221)
img_t     = flipud(voxel2world_space(Vimage, voxel_ctr));
vox_t     = flipud(voxel2world_space(Vmask, voxel_ctr));
vox_t_GM  = flipud(voxel2world_space(vol_GM_mask, voxel_ctr));
vox_t_WM  = flipud(voxel2world_space(vol_WM_mask, voxel_ctr));
vox_t_CSF = flipud(voxel2world_space(vol_CSF_mask, voxel_ctr));
img_t = img_t/MRSCont.coreg.T1_max{kk};
img_montage = [img_t+0.225*vox_t, img_t+0.3*vox_t_GM, img_t+0.225*vox_t_WM, img_t+0.4*vox_t_CSF];

%%% 4. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
if ~GUI
    out = figure;
    set(gcf, 'Color', 'w');
else
    out = figure('Visible','off');
end
imagesc(img_montage);
text(floor(size(vox_t,2)/2), 15, 'Voxel', 'Color', [1 1 1], 'FontSize', 14, 'HorizontalAlignment', 'center');
text(floor(size(vox_t,2)) + floor(size(vox_t,2)/2), 15, 'GM', 'Color', [1 1 1], 'FontSize', 14, 'HorizontalAlignment', 'center');
text(2*floor(size(vox_t,2)) + floor(size(vox_t,2)/2), 15, 'WM', 'Color', [1 1 1], 'FontSize', 14, 'HorizontalAlignment', 'center');
text(3*floor(size(vox_t,2)) + floor(size(vox_t,2)/2), 15, 'CSF', 'Color', [1 1 1], 'FontSize', 14, 'HorizontalAlignment', 'center');
colormap('gray');
caxis([0 1])
axis equal;
axis tight;
axis off;
if ~GUI
    title(['Segmentation: ' filename_voxel fileext_voxel ' & '  filename_image fileext_image], 'Interpreter', 'none','FontSize', 16);
else
    title(['Segmentation: ' filename_voxel fileext_voxel ' & '  filename_image fileext_image], 'Interpreter', 'none','FontSize', 16,'Color', MRSCont.colormap.Foreground);
end


%%% 5. ADD TISSUE COMPOSITION %%%
text(floor(size(vox_t,2)/2), size(img_montage,1)-15, 'voxel fraction', 'Color', [1 1 1], 'FontSize', 14, 'HorizontalAlignment', 'center');
tmp1 = sprintf('%.2f', MRSCont.seg.tissue.fGM(kk));
text(floor(size(vox_t,2)) + floor(size(vox_t,2)/2), size(img_montage,1)-15, tmp1, 'Color', [1 1 1], 'FontSize', 14, 'HorizontalAlignment', 'center');
tmp1 = sprintf('%.2f', MRSCont.seg.tissue.fWM(kk));
text(2*floor(size(vox_t,2)) + floor(size(vox_t,2)/2), size(img_montage,1)-15, tmp1, 'Color', [1 1 1], 'FontSize', 14, 'HorizontalAlignment', 'center');
tmp1 = sprintf('%.2f', MRSCont.seg.tissue.fCSF(kk));
text(3*floor(size(vox_t,2)) + floor(size(vox_t,2)/2), size(img_montage,1)-15, tmp1, 'Color', [1 1 1], 'FontSize', 14, 'HorizontalAlignment', 'center');

%%% 6. ADD OSPREY LOGO %%%
if ~GUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end
end

   