function [sin_info] = get_sin_info(sin_file)
%% [sin_info] = get_sin_info(sin_file)
%   Loads a Philips *.sin file and extracts the header information.
%
%   Input:
%       sin_file    *.sin file to be loaded
%
%   Output:
%       sin_info    struct containing the header info
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-03-15)
%       goeltzs1@jhmi.edu
%   
%   Credits:    
%       This code is inspired by code from:
%       Dr. Chenguang Zhao (chenguang.z.zhao@philips.com)
%       Philips Healthcare Suzhou
%
%   History:
%       2018-03-15: First version of the code.
%

% Open *.sin file
sin_info = [];
fid = fopen(sin_file, 'r');
data = fread(fid, '*char')';
fclose(fid);

% Search for regular expression coding for...

% channel numbers and channel names
expn = '\s(\d+)\s\d\d\s\d\d\:\schannel_names\s*:\s(\S+)[\r\n]';
tokens = regexp(data, expn, 'tokens');
for ii = 1:length(tokens)
    sin_info.acq_nos(ii) = str2double(tokens{1,ii}{1,1});
    sin_info.channel_names{ii} = tokens{1,ii}{1,2};
end

% voxel sizes
expn = '\s(\d\d)\s\d\d\s\d\d\:\svoxel_sizes\s*:(\s*[0-9]+.[0-9]+)+[\r\n]';
tokens = regexp(data, expn, 'tokens');
sin_info.voxel_size = str2num(tokens{1}{2}); % row - column - slice

% voxel sizes
expn = '\s\d\d\s\d\d\s(\d\d)\:\sloc_ap_rl_fh_offcentres\s*:(\s*-?\d+\.\d+)+[\r\n]';
tokens = regexp(data, expn, 'tokens');
sin_info.voxel_offsets = str2num(tokens{1}{2}); % AP - RL - FH

% slice orientations
expn = '\s\d\d\s\d\d\s(\d\d)\:\sslice_orientations\s*:(\s*[0-9]+)+[\r\n]';
tokens = regexp(data, expn, 'tokens');
sin_info.slice_orientation = str2num(tokens{1}{2}); % 0 = axial slices, 1 = sagittal slices, 2 = coronal slices

end
