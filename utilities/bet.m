function [data_brain, brain_mask] = bet(varargin)
% This function does a  brain segmentation for diffusion data.
% It uses a median filter smoothing of the b0 volumes and an automatic histogram Otsu thresholding.

% Input
% ----------
% data: diffusion data, (X, Y, Z, T).
% bvals: b-values, (X, Y, Z, T).
% median_radius : radius (in voxels) of the applied median filter. Default: 4.
% num_pass: number of pass of the median filter. Default: 4.
% dilate: number of iterations for binary dilation. Default: 0.

% Output
% -------
% data_brain : (X, Y, Z, T) 
% brain_mask : (X, Y, Z) 

p = inputParser;
p.addParameter('data', []);
p.addParameter('bvals', []);
p.addParameter('median_radius', 4);
p.addParameter('num_pass', 4);
p.addParameter('dilate', 0);
p.parse(varargin{:});
diffusion_data = p.Results.data;
bvals = p.Results.bvals;
median_radius = p.Results.median_radius;
num_pass = p.Results.num_pass;
dilate = p.Results.dilate;

b0_idx = bvals < 50;
% Get the mean of b0 volumes.
b0_vol = mean(diffusion_data(:,:,:,b0_idx),4);

% Size of the median window in each dimension.
medarr = ones(1,3) * ((median_radius * 2) + 1);

% Apply median filter multiple times on data.
b0_vol_filtered = b0_vol;
for i = 1:num_pass
    b0_vol_filtered = medfilt3(b0_vol_filtered, medarr);
end


% Calculate threshold value based on Otsu's method.
nbins = 256;
[hist, bin_centers] = histcounts(b0_vol_filtered, nbins);

% class probabilities for all possible thresholds
weight1 = cumsum(hist);
weight2 = flip(cumsum(flip(hist)));

% class means for all possible thresholds
mean1 = cumsum(hist .* bin_centers(2:end)) ./ weight1;
mean2 = flip(cumsum(flip((hist .* bin_centers(2:end)))) ./ flip(weight2));

% Clip ends to align class 1 and class 2 variables:
% The last value of weight1/mean1 should pair with zero values in weight2/mean2, which do not exist.
variance12 = weight1(1:end-1) .* weight2(2:end) .* (mean1(1:end-1) - mean2(2:end)).^2;

[~, idx] = max(variance12);
thresh = bin_centers(idx);

% Get the brain mask by thresholding
brain_mask = b0_vol_filtered > thresh;

if dilate
    cross(:,:,1) = [0,0,0;0,1,0;0,0,0];
    cross(:,:,2) = [0,1,0;1,1,1;0,1,0];
    cross(:,:,3) = [0,0,0;0,1,0;0,0,0];
    for i = 1:dilate
        brain_mask = imdilate(brain_mask, cross);
    end
end

data_brain = diffusion_data .* brain_mask;


end

