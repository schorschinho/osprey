function [cpx_file, Ref, img_array, img_body ,noise_array, noise_body] = loadRefScan(cpx_file)
%% [cpx_file ,Ref, img_array, img_body, noise_array, noise_body] = loadRefScan(cpx_file)
%   Loads a Philips coil reference scan (*.cpx format).
%
%   Input:
%       cpx_file    Coil sensitivity scan (*.cpx).
%
%   Output:
%       ref_file    cpx file used
%       Ref         struct containing dimensional information
%       img_array   4D stack of receiver coil images
%                   [number of coils, x, y, z]
%       img_body    4D stack of body coil images
%       noise_array 2D stack of receiver coil noise signals
%                   [number of coils, number of pixels estimated for noise]
%       noise_body  2D stack of body coil noise signals
%
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-03-15)
%       goeltzs1@jhmi.edu
%   
%   Credits:    
%       This code is based on an initial PRIAM reconstruction routine.
%       Dr. Vincent O. Boer (vincentob@drcmr.dk)
%       Danish Research Centre for Magnetic Resonance (Hvidovre Hospital)
%
%   History:
%       2018-03-15: First version of the code.
%


% Check whether .cpx file exists
if ~isfile(cpx_file) 
    msg = 'Coil sensitivity file (*.cpx) not found.';
    error(msg);   
end

%% Read coil reference imaging scan (*.cpx)

% Load the selected cpx file
fprintf('Loading coil reference scan...\n%s', cpx_file);
img = read_cpx(cpx_file);

% Permute the coil reference image array to be in the following order:
% 1 - number of coils
% 2 - image dimension x
% 3 - image dimension y
% 4 - image dimension z
% 5 - number of stacks:
%   1st stack contains receiver coil images
%   2nd stack contains body coil images (and is empty if less body coil
%   elements than receiver coil images)
img = permute(img,[5 1 2 4 3]);
Ref.ncoils = size(img,1);
Ref.dimx = size(img,2);
Ref.dimy = size(img,3);
Ref.dimz = size(img,4);

% Separate into two separate arrays
img_array = double(img(:,:,:,:,1));
img_body = double(img(:,:,:,:,2));

% Normalize within each array
img_array = img_array/max(abs(img_array(:)));
img_body = img_body/max(abs(img_body(:)));

% If there is a *.raw file for the reference image, extract the
% noise from the header.
if exist([cpx_file(1:end-4) '.raw'],'file')
    disp('Found *.raw file for the coil reference scan. Extracting noise from the *.raw header.');
    [~, info] = read_noise([cpx_file(1:end-4) '.raw']);
    noise_array = info.FRC_NOISE_DATA(:,:,1);
    noise_body = info.FRC_NOISE_DATA(:,:,2);
else
    % If not, extract noise levels from the edge of the first slice.
    disp('No *.raw file for the coil reference scan found. Extracting noise from the edges of the image.');
    noise_array = squeeze(img_array(:,1:30,:,1));
    noise_array = reshape(noise_array,[Ref.ncoils size(noise_array,2)*size(noise_array,3)]);
    noise_body = squeeze(img_body(:,1:30,:,1));
    noise_body = reshape(noise_body,[Ref.ncoils size(noise_body,2)*size(noise_body,3)]);
end

disp('Loading coil reference scan finished!');