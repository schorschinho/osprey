function out = osp_plotVoxelOverlap(MRSCont,voxelCenter)
%% out = osp_plotVoxelOverlap(MRSCont)
%   Creates a figure showing overlap of all MRS voxels in spm152 space
%
%   USAGE:
%       out = osp_plotVoxelOverlap(MRSCont)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   INPUTS:
%       MRSCont  = Osprey data container.
%       voxelCenter = center coordinates to be used for the three plane
%       image (optinal). If no coordinates are supplied the center of mass
%       of all voxels is calculated to identify the voxel center.
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2022-07-15)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2022-07-15: First version of the code.
%
% For the SPM152 template we are using the MRIcroGL spm152 template.

% Check that OspreySeg has been run before
if ~MRSCont.flags.didSeg
    error('Trying to plot voxel overlap, but data has not been processed yet. Run OspreySeg first.')
end

%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided

if nargin<1
    error('ERROR: no input Osprey container specified.  Aborting!!');
end
if nargin<2
    voxelCenter =[];
end

%%% 2. LOAD DATA TO PLOT %%%
if ~(isfield(MRSCont.flags,'addImages') && (MRSCont.flags.addImages == 1))
    % Load spm152 template and mask overlap
     template = which('osprey/libraries/MRIcroGL/templates/spm152.nii');
    [~,filename_image,fileext_image]   = fileparts(template);
    [~,filename_mask,~]   = fileparts(MRSCont.seg.overlapfile);

    if ~exist(MRSCont.seg.overlapfile,'file')
        gunzip([MRSCont.seg.overlapfile, '.gz']);
    end

    %%% 3. SET UP THREE PLANE IMAGE %%%
    [three_plane_img,three_plane_overlay,vox_t_size] = osp_extract_three_plane_image_overlay(template, MRSCont.seg.overlapfile,voxelCenter,90);

    %%% 4. SET UP FIGURE LAYOUT %%%
    % Generate a new figure and keep the handle memorized

    if ~MRSCont.flags.didSeg
        if exist([MRSCont.coreg.vol_mask{kk}.fname, '.gz'],'file')
            delete(MRSCont.coreg.vol_mask{kk}.fname);
        end
        if exist([MRSCont.coreg.vol_image{kk}.fname, '.gz'],'file')
            delete(MRSCont.coreg.vol_image{kk}.fname);
        end
    end
else
    [~,filename_voxel,fileext_voxel]   = fileparts(MRSCont.files{kk});
    [~,filename_image,fileext_image]   = fileparts(MRSCont.coreg.vol_image{kk}.fname);
    three_plane_img = MRSCont.coreg.three_plane_img{kk};
end


if ~MRSCont.flags.isGUI
    out = figure('Visible','on');
    set(gcf, 'Color', 'w');
else
    out = figure('Visible','off');
end

bar = repmat(linspace(0,1,size(three_plane_overlay,2)),[15 1]);
three_plane_overlay = vertcat(bar,three_plane_overlay);
three_plane_img = vertcat(0.7*ones(size(bar)),three_plane_img);
rgb = mat2im(three_plane_overlay,jet(1000));
for col = 1 : 3
    temp = rgb(:,:,col);
    temp(three_plane_overlay==0)=0;
    rgb(:,:,col) = temp;
end
ha= imshow(rgb);
axis equal;
axis tight;
axis off;
hold on;
hb = imshow(mat2gray(three_plane_img));
% caxis([0 mean(three_plane_img(three_plane_img > 0.01)) + 3*std(three_plane_img(three_plane_img > 0.01))]);
% colormap('gray');
AlphaChannel = three_plane_overlay;

set(hb,'AlphaData',1-AlphaChannel);

text(floor(vox_t_size/4), size(three_plane_img,1)-15, 'L', 'Color', MRSCont.colormap.Background, 'FontSize', 14, 'HorizontalAlignment', 'center');
text(floor(vox_t_size), size(three_plane_img,1)-15, 'R', 'Color', MRSCont.colormap.Background, 'FontSize', 14, 'HorizontalAlignment', 'center');
text(1, 25+15, 'Overlap [%]', 'Color', MRSCont.colormap.Background, 'FontSize', 16, 'HorizontalAlignment', 'left');
text(1, 8+15, '0', 'Color', MRSCont.colormap.Background, 'FontSize', 16, 'HorizontalAlignment', 'left');
text(size(three_plane_img,2)*1/4, 8+15, '25', 'Color', MRSCont.colormap.Background, 'FontSize', 16, 'HorizontalAlignment', 'center');
text(size(three_plane_img,2)*2/4, 8+15, '50', 'Color', MRSCont.colormap.Background, 'FontSize', 16, 'HorizontalAlignment', 'center');
text(size(three_plane_img,2)*3/4, 8+15, '75', 'Color', MRSCont.colormap.Background, 'FontSize', 16, 'HorizontalAlignment', 'center');
text(size(three_plane_img,2)-25, 8+15, '100', 'Color', MRSCont.colormap.Background, 'FontSize', 16, 'HorizontalAlignment', 'left');


if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
    titleStr = sprintf(['Voxel Consistency:\n ' filename_mask ' & SPM 152 template']);
else
    titleStr = sprintf(['Coregistration:\n ' filename_voxel fileext_voxel ' & '  filename_image fileext_image '\n Voxel ' num2str(VoxelIndex)]);      
end

if ~MRSCont.flags.isGUI
        title(titleStr, 'Interpreter', 'none','FontSize', 16);
else
    title(titleStr, 'Interpreter', 'none','FontSize', 16,'Color', MRSCont.colormap.Foreground);
end

sz = get(gcf, 'Position');
set(gcf,'Position',[sz(1) sz(2) 850 400])

%%% 5. ADD OSPREY LOGO %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end

if exist(MRSCont.seg.overlapfile,'file')
    gzip(MRSCont.seg.overlapfile);
    delete(MRSCont.seg.overlapfile);
end

end

   