function [off_spec_no_lb,on_spec_no_lb] = load_mrsi_data(FileName, PathName,lb, spec_zfill, k_zfill, seq_type, water_filter, k_fft2_wat_ref, k_fft2_wat_ref_no_k_zfill, k_sort_b4_2dfft, wat_file)
%% [off_spec_no_lb,on_spec_no_lb] = load_mrsi_data(MRSCont)
%   This functions load data/list MRSI data (k-space) and performs motion correction
%   Motion correciton is based on correlation coefficents between motion
%   and non-motion affected data
%
%   USAGE:
%       [off_spec_no_lb,on_spec_no_lb] = OspreyProcess(MRSCont);
%
%   INPUTS:
%       FileName     = name of the data file.
%       PathName     = name of the path to the data file.
%       lb           = exponential linebroadening value used during
%       correction for Osprey we are exporting spectra without LB
%       spec_zfill   = zero filling factor, for Osprey we are exporting
%       k_zfill      = zero filling facotor for k space
%       spectra without zero filling
%       seq_type      = sequence type used in the MRSI acqusiton 
%                       OPTIONS:    - PRESS
%                                   - Water Reference
%                                   - SE multislice
%                                   - water multislice
%                                   - MEGA-PRESS
%                                   - MEGA multislice 
%                                   - HERMES
%                                   - HERMES lip sup
%       water_filter  = flag for HSVD water removal
%       k_fft2_wat_ref = no idea what this is
%       k_fft2_wat_ref_no_k_zfill = no idea what this is
%       k_sort_b4_2dfft = no idea what this is
%       wat_file      = water reference scan file 
%
%   OUTPUTS:
%       [off_spec_no_lb,on_spec_no_lb]     = corrected on and off spectra
%
%   AUTHOR:
%       Dr. Kimberly Chan
%       
%       Dr.Helge Zoellner (Johns Hopkins University, 2021-01-01)
%       hzoelln2@jhmi.edu
%
%   CREDITS:
%       This code is based on numerous functions Dr. Kimberly Chan (UT
%       Southwestern)
%
%   HISTORY:
%       2021-02-01: Adaptions for Osprey.


tic
disp('Opening data.')
this_file = [PathName, FileName];



endian = 'l';
type = 'float';

disp('Reading scan parameters.')
fname_scan_params = [this_file(1:(end-4)),'list'];
scan_params = textread(fname_scan_params, '%s');
% Find the data lines and ignore the noise channels.
data_lines = find(strcmp(scan_params,'STD'));
tot_offset_idx = scan_params(data_lines(end) + 20); % last STD + 20 more offsets
data_lines = data_lines(3:end); % ignore first 4 STDs.
offset = str2num(scan_params{data_lines(1) + 20});

tot_offsets = (str2num(tot_offset_idx{1}) - offset)/8192 + 1; % All offsets divided by 8192 bytes.


disp('Reading data.')
fp=fopen(this_file, 'rb', endian);

fseek(fp, offset, -1); % Start reading from the offset

data_raw = fread(fp, 2048*tot_offsets, type); % 1024 real points, 2048 complex points
                                       % size of each floating point. (unsigned
                                       % integer?
data=data_raw(1:2:end,:)+1i*data_raw(2:2:end,:); % Points alternate between real and complex
data = reshape(data,[1024, tot_offsets]); % Reshape back to 1024 points by 8384 offsets.



coil = zeros(1,length(data_lines - 1));
kx = coil;
ky = coil;
avg = coil;
sign = coil;
loc = coil;
count = 0;
all_count = kx;
for dl = 1:(length(data_lines))
    if dl == length(data_lines)
        this_dl = cellfun(@str2num,{scan_params{(data_lines(dl)+1):(data_lines(dl)+20)}});
    else
        this_dl = cellfun(@str2num,{scan_params{(data_lines(dl)+1):(data_lines(dl + 1)-1)}});
    end
    
    count = count + 1;
    all_count(dl) = count;
    loc(dl) = this_dl(5) + 1;
    coil(dl) = this_dl(6);
    kx(dl) = this_dl(9);
    ky(dl) = this_dl(10);
    avg(dl) = this_dl(12);
    sign(dl) = this_dl(13);
end

kx_tot = max(kx) - min(kx) + 1;
ky_tot = max(ky) - min(ky) + 1;

x_tot = k_zfill*kx_tot;
y_tot = k_zfill*ky_tot;

if (strcmp(seq_type, 'MEGA multislice') || strcmp(seq_type, 'Water multislice') || strcmp(seq_type, 'SE multislice'))
    k_sort_on = zeros(max(loc), kx_tot, ky_tot, max(coil) + 1, 1024);
else
    k_sort_on = zeros(kx_tot, ky_tot, max(coil) + 1, 1024);
end

k_sort_off = k_sort_on;
k_sort_on2 = k_sort_on;
k_sort_off2 = k_sort_off;
count_sort = zeros(size(k_sort_off2,1), size(k_sort_off2,2));
on_k_count = zeros(kx_tot, ky_tot, 2);
off_k_count = zeros(kx_tot, ky_tot, 2);

if (strcmp(seq_type, 'HERMES') || strcmp(seq_type, 'HERMES lip sup'))
    k_sort_on1_on2 = zeros(kx_tot, ky_tot, max(coil) + 1, 1024);
    k_sort_on1_off2 = k_sort_on1_on2;
    k_sort_off1_on2 = k_sort_on1_on2;
    k_sort_off1_off2 = k_sort_on1_on2;
end


disp('Reorganizing k-space locations.')

k_count = 0;

last_kx = 500;
last_ky = 500;
avg_num_k = zeros(kx_tot, ky_tot);
k_tot = 0;
% Rearrange k-space values (so no negative indices)
for dl = 1:(length(data_lines))
    this_kx = kx(dl) + abs(min(kx)) + 1;
    this_ky = ky(dl) + abs(min(ky)) + 1; 
    
    this_avg = avg(dl);
    this_coil = coil(dl);
    this_loc = loc(dl);
%   
    if this_kx ~= last_kx || this_ky ~= last_ky
        k_count = k_count + 1;
        last_kx = this_kx;
        last_ky = this_ky;
        count_sort(this_kx, this_ky) = k_count;
    end

    if strcmp(seq_type, 'MEGA-PRESS')        % MEGA-PRESS
       
        if this_avg == 0 
            k_sort_on(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
            on_k_count(this_kx, this_ky, 1) = 1;
            if this_coil == 1
                k_tot = k_tot + 1;
                avg_num_k(this_kx, this_ky) = k_tot;
            end
        elseif this_avg == 1 
            k_sort_off(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
            off_k_count(this_kx, this_ky, 1) = 1;
        elseif this_avg == 2
            k_sort_on2(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
            on_k_count(this_kx, this_ky, 2) = 1;
        elseif this_avg == 3
            k_sort_off2(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
            off_k_count(this_kx, this_ky, 2) = 1;
        end
    elseif strcmp(seq_type, 'HERMES')        % HERMES
        switch this_avg
            case 0
               k_sort_off1_on2(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
            case 1
               k_sort_on1_off2(this_kx, this_ky, this_coil + 1, :) = data(:,dl); %x
            case 2
               k_sort_on1_on2(this_kx, this_ky, this_coil + 1, :) = data(:,dl); %x
            case 3
               k_sort_off1_off2(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        end
    elseif strcmp(seq_type, 'HERMES lip sup')        % HERMES lipid suppression
        switch this_avg
            case 0
               k_sort_on1_off2(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
            case 1
               k_sort_off1_off2(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
            case 2
               k_sort_on1_on2(this_kx, this_ky, this_coil + 1, :) = data(:,dl); 
            case 3
               k_sort_off1_on2(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        end
    elseif strcmp(seq_type, 'Water Reference')
        k_sort_off(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        k_sort_on(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        k_sort_off2(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        k_sort_on2(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
    elseif strcmp(seq_type, 'PRESS')
        switch this_avg
            case 0
                k_sort_off(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
                k_sort_on(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
            case 1
                k_sort_off2(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
                k_sort_on2(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        end
    elseif strcmp(seq_type, 'MEGA multislice')        % MEGA Multislice, default is 3 slices
        switch this_avg
            case 0
                k_sort_on(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
            case 1
                k_sort_off(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
            case 2
                k_sort_on2(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
            case 3
                k_sort_off2(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        end
    elseif strcmp(seq_type, 'SE multislice')        % SE Multislice, default is 3 slices   
        k_sort_off(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        k_sort_on(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        k_sort_off2(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        k_sort_on2(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);

    elseif strcmp(seq_type, 'Water multislice')        % MEGA Multislice, default is 3 slices
        k_sort_on(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        k_sort_off(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        k_sort_on2(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
        k_sort_off2(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
    end
end

% k_sort_on(7,8,:,:) = 0;
% k_sort_off(7,8,:,:) = 0;

% ------------------------------------------------------------------------------------------
% --------------- %
% Hanning Filter  %
% --------------- %
% Apply a Hanning filter on k-space data to improve the PSF

if (~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'Water multislice') && ~strcmp(seq_type, 'SE multislice'))
    hanning_x = repmat(hanning(kx_tot), [1 ky_tot max(coil) + 1 1024]);
    hanning_y = permute(repmat(hanning(ky_tot), [1 kx_tot max(coil) + 1 1024]), [2 1 3 4]);
else
    hanning_x = repmat(hanning(kx_tot), [1 max(loc) ky_tot max(coil) + 1 1024]);
    hanning_y = permute(repmat(hanning(ky_tot), [1 kx_tot max(loc) max(coil) + 1 1024]), [3 2 1 4 5]);
    hanning_x = permute(hanning_x, [2 1 3 4 5]);
end


k_sort_on = k_sort_on.*hanning_x;
k_sort_on = k_sort_on.*hanning_y;
k_sort_on2 = k_sort_on2.*hanning_x;
k_sort_on2 = k_sort_on2.*hanning_y;
k_sort_off = k_sort_off.*hanning_x;
k_sort_off = k_sort_off.*hanning_y;
k_sort_off2 = k_sort_off2.*hanning_x;
k_sort_off2 = k_sort_off2.*hanning_y;

if (strcmp(seq_type, 'HERMES') || strcmp(seq_type, 'HERMES lip sup'))
    k_sort_on1_on2 = k_sort_on1_on2.*hanning_x;
    k_sort_on1_off2 = k_sort_on1_off2.*hanning_x;
    k_sort_off1_on2 = k_sort_off1_on2.*hanning_x;
    k_sort_off1_off2 = k_sort_off1_off2.*hanning_x;
    
    k_sort_on1_on2 = k_sort_on1_on2.*hanning_y;
    k_sort_on1_off2 = k_sort_on1_off2.*hanning_y;
    k_sort_off1_on2 = k_sort_off1_on2.*hanning_y;
    k_sort_off1_off2 = k_sort_off1_off2.*hanning_y;
end

% ------------------------------------------------------------------------------------------
% ------------------%
% Analyzing K-space %
% ------------------%

if ~isempty(k_fft2_wat_ref_no_k_zfill)
    if (strcmp(seq_type, 'MEGA multislice') || strcmp(seq_type, 'Water multislice') || strcmp(seq_type, 'SE multislice'))

        wat_on_peak1 = squeeze(k_sort_b4_2dfft(:,:,:,:,1));
        wat_off_peak1 = squeeze(k_sort_b4_2dfft(:,:,:,:,1));
        wat_on_peak2 = squeeze(k_sort_b4_2dfft(:,:,:,:,1));
        wat_off_peak2 = squeeze(k_sort_b4_2dfft(:,:,:,:,1));
    else
        wat_on_peak1 = squeeze(k_fft2_wat_ref_no_k_zfill(:,:,:,1));
        wat_off_peak1 = squeeze(k_fft2_wat_ref_no_k_zfill(:,:,:,1));
        wat_on_peak2 = squeeze(k_fft2_wat_ref_no_k_zfill(:,:,:,1));
        wat_off_peak2 = squeeze(k_fft2_wat_ref_no_k_zfill(:,:,:,1));
    end
else
    wat_on_peak1 = squeeze(k_sort_on(:,:,:,1));
    wat_off_peak1 = squeeze(k_sort_off(:,:,:,1));
    wat_on_peak2 = squeeze(k_sort_on2(:,:,:,1));
    wat_off_peak2 = squeeze(k_sort_off2(:,:,:,1));
    if (strcmp(seq_type, 'MEGA multislice') || strcmp(seq_type, 'Water multislice') || strcmp(seq_type, 'SE multislice'))
        wat_on_peak1 = squeeze(k_sort_on(:,:,:,:,1));
        wat_off_peak1 = squeeze(k_sort_off(:,:,:,:,1));
        wat_on_peak2 = squeeze(k_sort_on2(:,:,:,:,1));
        wat_off_peak2 = squeeze(k_sort_off2(:,:,:,:,1));
    end   
end

if (~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'Water multislice') && ~strcmp(seq_type, 'SE multislice'))
    % Phase each coil before summing over the channels.
    wat_on_peak1 = repmat(wat_on_peak1, [1 1 1 1024]);
    wat_off_peak1 = repmat(wat_off_peak1, [1 1 1 1024]);
    wat_on_peak2 = repmat(wat_on_peak2, [1 1 1 1024]);
    wat_off_peak2 = repmat(wat_off_peak2, [1 1 1 1024]);
else
   % Phase each coil before summing over the channels.
    wat_on_peak1 = repmat(wat_on_peak1, [1 1 1 1 1024]);
    wat_off_peak1 = repmat(wat_off_peak1, [1 1 1 1 1024]);
    wat_on_peak2 = repmat(wat_on_peak2, [1 1 1 1 1024]);
    wat_off_peak2 = repmat(wat_off_peak2, [1 1 1 1 1024]);
end

k_sort_on_phased1 = k_sort_on.*conj(wat_off_peak1)./abs(wat_off_peak1); 
k_sort_off_phased1 = k_sort_off.*conj(wat_off_peak1)./abs(wat_off_peak1);
k_sort_on_phased2 = k_sort_on2.*conj(wat_off_peak2)./abs(wat_off_peak2); 
k_sort_off_phased2 = k_sort_off2.*conj(wat_off_peak2)./abs(wat_off_peak2);


if (~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'Water multislice') && ~strcmp(seq_type, 'SE multislice'))
    k_sort_on_phased1 = squeeze(sum(k_sort_on_phased1,3));
    k_sort_off_phased1 = squeeze(sum(k_sort_off_phased1,3));
    k_sort_on_phased2 = squeeze(sum(k_sort_on_phased2,3));
    k_sort_off_phased2 = squeeze(sum(k_sort_off_phased2,3));
else
    k_sort_on_phased1 = squeeze(sum(k_sort_on_phased1,4));
    k_sort_off_phased1 = squeeze(sum(k_sort_off_phased1,4));
    k_sort_on_phased2 = squeeze(sum(k_sort_on_phased2,4));
    k_sort_off_phased2 = squeeze(sum(k_sort_off_phased2,4));
end
% 

if (~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'Water multislice') && ~strcmp(seq_type, 'SE multislice'))
    k_sort_on_phased_k1 = k_sort_on_phased1;
    k_sort_off_phased_k1 = k_sort_off_phased1;
    k_sort_on_phased_k2 = k_sort_on_phased2;
    k_sort_off_phased_k2 = k_sort_off_phased2;
    
    motion_corr_lb = 5;
    exp_lb = permute(squeeze((repmat(exp(-(motion_corr_lb*pi*(1:1024))/1024), [1 1 size(k_sort_on_phased2,1)...
        size(k_sort_on_phased2,2)]))), [2 3 1]);
    
    % Form spectra in k-space.
    
    on_spec_k_1 = fftshift(fft(squeeze(k_sort_on_phased_k1).*exp_lb,1024*spec_zfill,3),3);
    off_spec_k_1 = fftshift(fft(squeeze(k_sort_off_phased_k1).*exp_lb,1024*spec_zfill,3),3);
    on_spec_k_2 = fftshift(fft(squeeze(k_sort_on_phased_k2).*exp_lb,1024*spec_zfill,3),3);
    off_spec_k_2 = fftshift(fft(squeeze(k_sort_off_phased_k2).*exp_lb,1024*spec_zfill,3),3);
    
    % -------------- Motion Correction Identify ------------%
    % Identify motion & phase corrections                   %
    % ------------------------------------------------------%
    
    if strcmp(seq_type, 'MEGA-PRESS')
        
%         % Remove only
%         [k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
%         on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2] = motion_correct_mrsi_full(k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
%                                                                                     on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,...
%                                                                                     spec_zfill, seq_type, kx_tot, ky_tot);
        
        % If I don't want to delete k-space data or correct something HZ
           k_ph_corr_on1 = zeros(size(k_sort_on2,1),size(k_sort_on2,2));
           k_ph_corr_on2 = zeros(size(k_sort_on2,1),size(k_sort_on2,2));
           k_ph_corr_off1 = zeros(size(k_sort_on2,1),size(k_sort_on2,2));
           k_ph_corr_off2 = zeros(size(k_sort_on2,1),size(k_sort_on2,2));   
        
        
        % Remove + phase - Used in paper
%         [k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
%         on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,...
%         k_ph_corr_on1, k_ph_corr_on2, k_ph_corr_off1, k_ph_corr_off2,...
%         on1_replace_track, off1_replace_track, on2_replace_track, ...
%         off2_replace_track, zero_replace_track, k_space_locs, corr_options]  = motion_correct_mrsi_full_ph(k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
%                                                                                                         on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,...
%                                                                                                         spec_zfill, seq_type, kx_tot, ky_tot);

%         % Remove + phase streamlined
%         [k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
%         on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,...
%         k_ph_corr_on1, k_ph_corr_on2, k_ph_corr_off1, k_ph_corr_off2,...
%         on1_replace_track, off1_replace_track, on2_replace_track, ...
%         off2_replace_track, zero_replace_track, k_space_locs]  = motion_correct_mrsi_ph_streamline2(k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
%                                                                                                         on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,...
%                                                                                                         spec_zfill, seq_type, kx_tot, ky_tot);
    end
else
    k_sort_on_phased_k1 = k_sort_on_phased1;
    k_sort_off_phased_k1 = k_sort_off_phased1;
    k_sort_on_phased_k2 = k_sort_on_phased2;
    k_sort_off_phased_k2 = k_sort_off_phased2;
    
    sz_k_sort = size(k_sort_off_phased_k2);
    exp_lb = permute(squeeze((repmat(exp(-(lb*pi*(1:1024))/1024), [1 1 1 sz_k_sort(1:3)]))), [2 3 4 1]);

    % fft the fIDs 
    on_spec_k_1 = (fftshift(fft(squeeze(k_sort_on_phased_k1).*exp_lb,1024*spec_zfill,4),4));
    off_spec_k_1 = (fftshift(fft(squeeze(k_sort_off_phased_k1).*exp_lb,1024*spec_zfill,4),4));
    on_spec_k_2 = (fftshift(fft(squeeze(k_sort_on_phased_k2).*exp_lb,1024*spec_zfill,4),4));
    off_spec_k_2 = (fftshift(fft(squeeze(k_sort_off_phased_k2).*exp_lb,1024*spec_zfill,4),4));
    if strcmp(seq_type, 'MEGA multislice')
        %k_ph_corr_on1, k_ph_corr_off1, k_ph_corr_on2, k_ph_corr_off2,...
        [k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
         on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2] = motion_correct_mrsi(k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
                                                                                    on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,...
                                                                                    spec_zfill, seq_type, kx_tot, ky_tot);
    end
end

ppm_axis = linspace(-1000,1000,1024*spec_zfill)/128 + 4.68;

if (strcmp(seq_type, 'HERMES') || strcmp(seq_type, 'HERMES lip sup'))
    
    % For each point in time and coil take the 2D fft
    sz_k_sort_on = size(k_sort_on1_on2);
    sz_k_sort_off = size(k_sort_off1_off2);

    k_fft2_on1_on2 = zeros([x_tot,y_tot,sz_k_sort_on(3:end)]);
    k_fft2_on1_off2 = zeros([x_tot,y_tot,sz_k_sort_off(3:end)]);
    k_fft2_off1_on2 = zeros([x_tot,y_tot,sz_k_sort_on(3:end)]);
    k_fft2_off1_off2 = zeros([x_tot,y_tot,sz_k_sort_off(3:end)]);

    sz_k_on = size(k_sort_on1_on2);
    k_fft2_on1_on2_wat = zeros([x_tot,y_tot, sz_k_on(3)]);
    k_fft2_on1_off2_wat = zeros([x_tot,y_tot, sz_k_on(3)]);
    k_fft2_off1_on2_wat = zeros([x_tot,y_tot, sz_k_on(3)]);
    k_fft2_off1_off2_wat = zeros([x_tot,y_tot, sz_k_on(3)]);
    
    wat_k_space_on1_on2 = zeros(sz_k_on(1:3));
    wat_k_space_on1_off2 = zeros(sz_k_on(1:3));
    wat_k_space_off1_on2 = zeros(sz_k_on(1:3));
    wat_k_space_off1_off2 = zeros(sz_k_on(1:3));

    disp('Taking Fourier transforms.')
    for c_idx = 1:size(k_sort_on1_on2,3) % Each coil
        wat_on1_on2_peak = squeeze(k_sort_on1_on2(:,:,c_idx,1));
        wat_on1_off2_peak = squeeze(k_sort_on1_off2(:,:,c_idx,1));
        wat_off1_on2_peak = squeeze(k_sort_off1_on2(:,:,c_idx,1));
        wat_off1_off2_peak = squeeze(k_sort_off1_off2(:,:,c_idx,1));
        
        k_fft2_on1_on2_wat(:,:,c_idx) = fft2(wat_on1_on2_peak, x_tot,y_tot);
        k_fft2_on1_off2_wat(:,:,c_idx) = fft2(wat_on1_off2_peak, x_tot,y_tot);
        k_fft2_off1_on2_wat(:,:,c_idx) = fft2(wat_off1_on2_peak, x_tot,y_tot);
        k_fft2_off1_off2_wat(:,:,c_idx) = fft2(wat_off1_off2_peak, x_tot,y_tot);

        wat_k_space_on1_on2(:,:,c_idx) = wat_on1_on2_peak.*conj(wat_on1_on2_peak)./abs(wat_on1_on2_peak);
        wat_k_space_on1_off2(:,:,c_idx) = wat_on1_off2_peak.*conj(wat_on1_off2_peak)./abs(wat_on1_off2_peak);
        wat_k_space_off1_on2(:,:,c_idx) = wat_off1_on2_peak.*conj(wat_off1_on2_peak)./abs(wat_off1_on2_peak);
        wat_k_space_off1_off2(:,:,c_idx) = wat_off1_off2_peak.*conj(wat_off1_off2_peak)./abs(wat_off1_off2_peak);
        
        for t_idx = 1:size(k_sort_on1_on2, 4) % Each point in time
            k_fft2_on1_on2(:,:,c_idx, t_idx) =  fft2(squeeze(k_sort_on1_on2(:,:,c_idx,t_idx)),x_tot,y_tot); % zerofill k-space
            k_fft2_on1_off2(:,:,c_idx, t_idx) = fft2(squeeze(k_sort_on1_off2(:,:,c_idx,t_idx)),x_tot,y_tot); 
            k_fft2_off1_on2(:,:,c_idx, t_idx) = fft2(squeeze(k_sort_off1_on2(:,:,c_idx,t_idx)),x_tot,y_tot);
            k_fft2_off1_off2(:,:,c_idx, t_idx) = fft2(squeeze(k_sort_off1_off2(:,:,c_idx,t_idx)),x_tot,y_tot); 
        end
    end


    wat_k_space_on1_on2 = sum(wat_k_space_on1_on2,3);
    wat_k_space_on1_off2 = sum(wat_k_space_on1_off2,3);
    wat_k_space_off1_on2 = sum(wat_k_space_off1_on2,3);
    wat_k_space_off1_off2 = sum(wat_k_space_off1_off2,3);


    %wat_map = squeeze(abs(fftshift(sum(k_fft2_on_wat,3))));

    % Phase each coil before summing over the channels.
    k_fft2_on1_on2_wat = repmat(k_fft2_on1_on2_wat, [1 1 1 1024]);
    k_fft2_on1_off2_wat = repmat(k_fft2_on1_off2_wat, [1 1 1 1024]);
    k_fft2_off1_on2_wat = repmat(k_fft2_off1_on2_wat, [1 1 1 1024]);
    k_fft2_off1_off2_wat = repmat(k_fft2_off1_off2_wat, [1 1 1 1024]);
    
    if isempty(k_fft2_wat_ref)
        k_fft2_on1_on2_phased = k_fft2_on1_on2.*conj(k_fft2_on1_on2_wat)./abs(k_fft2_on1_on2_wat);
        k_fft2_on1_off2_phased = k_fft2_on1_off2.*conj(k_fft2_on1_off2_wat)./abs(k_fft2_on1_off2_wat);
        k_fft2_off1_on2_phased = k_fft2_off1_on2.*conj(k_fft2_off1_on2_wat)./abs(k_fft2_off1_on2_wat);
        k_fft2_off1_off2_phased = k_fft2_off1_off2.*conj(k_fft2_off1_off2_wat)./abs(k_fft2_off1_off2_wat);
    else
        k_fft2_on1_on2_phased = k_fft2_on1_on2.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
        k_fft2_on1_off2_phased = k_fft2_on1_off2.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
        k_fft2_off1_on2_phased = k_fft2_off1_on2.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
        k_fft2_off1_off2_phased = k_fft2_off1_off2.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
    end

    k_fft2_on1_on2 = squeeze(sum(k_fft2_on1_on2_phased,3));
    k_fft2_on1_off2 = squeeze(sum(k_fft2_on1_off2_phased,3));
    k_fft2_off1_on2 = squeeze(sum(k_fft2_off1_on2_phased,3));
    k_fft2_off1_off2 = squeeze(sum(k_fft2_off1_off2_phased,3));

elseif (~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'Water multislice') && ~strcmp(seq_type, 'SE multislice'))
    
    % For each point in time and coil take the 2D fft
    sz_k_sort_on = size(k_sort_on);
    sz_k_sort_off = size(k_sort_off);

    k_fft2_on = zeros([x_tot,y_tot,sz_k_sort_on(3:end)]);
    k_fft2_on2 = zeros([x_tot,y_tot,sz_k_sort_on(3:end)]);
    k_fft2_off = zeros([x_tot,y_tot,sz_k_sort_off(3:end)]);
    k_fft2_off2 = zeros([x_tot,y_tot,sz_k_sort_off(3:end)]);

    sz_k_on = size(k_sort_on);
    k_fft2_on_wat = zeros([x_tot,y_tot, sz_k_on(3)]);
    k_fft2_on2_wat = zeros([x_tot,y_tot, sz_k_on(3)]);
    k_fft2_off_wat = zeros([x_tot,y_tot, sz_k_on(3)]);
    k_fft2_off2_wat = zeros([x_tot,y_tot, sz_k_on(3)]);
    
    wat_k_space_on = zeros(sz_k_on(1:3));
    wat_k_space_on2 = zeros(sz_k_on(1:3));
    wat_k_space_off = zeros(sz_k_on(1:3));
    wat_k_space_off2 = zeros(sz_k_on(1:3));
    
    % -------------- Motion Correction Phase ------------%
    % Phase correction (for motion) here                 %
    % ---------------------------------------------------%
% % 
    k_ph_corr_on1_rep = repmat(k_ph_corr_on1, [1 1 size(k_sort_on,3) size(k_sort_on, 4)]);
    k_ph_corr_on2_rep = repmat(k_ph_corr_on2, [1 1 size(k_sort_on,3) size(k_sort_on, 4)]);
    k_ph_corr_off1_rep = repmat(k_ph_corr_off1, [1 1 size(k_sort_on,3) size(k_sort_on, 4)]);
    k_ph_corr_off2_rep = repmat(k_ph_corr_off2, [1 1 size(k_sort_on,3) size(k_sort_on, 4)]);
%     
    k_sort_on = k_sort_on.*exp(1i*pi*k_ph_corr_on1_rep/180);
    k_sort_on2 = k_sort_on2.*exp(1i*pi*k_ph_corr_on2_rep/180);
    k_sort_off = k_sort_off.*exp(1i*pi*k_ph_corr_off1_rep/180);
    k_sort_off2 = k_sort_off2.*exp(1i*pi*k_ph_corr_off2_rep/180);

    disp('Taking Fourier transforms.')
    for c_idx = 1:size(k_sort_on,3) % Each coil
        wat_on_peak = squeeze(k_sort_on(:,:,c_idx,1));
        wat_on2_peak = squeeze(k_sort_on2(:,:,c_idx,1));
        wat_off_peak = squeeze(k_sort_off(:,:,c_idx,1));
        wat_off2_peak = squeeze(k_sort_off2(:,:,c_idx,1));
        
        k_fft2_on_wat(:,:,c_idx) = fft2(wat_on_peak, x_tot,y_tot);
        k_fft2_on2_wat(:,:,c_idx) = fft2(wat_on2_peak, x_tot,y_tot);
        k_fft2_off_wat(:,:,c_idx) = fft2(wat_off_peak, x_tot,y_tot);
        k_fft2_off2_wat(:,:,c_idx) = fft2(wat_off2_peak, x_tot,y_tot);

        wat_k_space_on(:,:,c_idx) = wat_on_peak.*conj(wat_on_peak)./abs(wat_on_peak);
        wat_k_space_on2(:,:,c_idx) = wat_on2_peak.*conj(wat_on2_peak)./abs(wat_on2_peak);
        wat_k_space_off(:,:,c_idx) = wat_off_peak.*conj(wat_off_peak)./abs(wat_off_peak);
        wat_k_space_off2(:,:,c_idx) = wat_off2_peak.*conj(wat_off2_peak)./abs(wat_off2_peak);
        
        for t_idx = 1:size(k_sort_on, 4) % Each point in time
            k_fft2_on(:,:,c_idx, t_idx) = fft2(squeeze(k_sort_on(:,:,c_idx,t_idx)),x_tot,y_tot); % zerofill k-space
            k_fft2_on2(:,:,c_idx, t_idx) = fft2(squeeze(k_sort_on2(:,:,c_idx,t_idx)),x_tot,y_tot); 
            k_fft2_off(:,:,c_idx, t_idx) = fft2(squeeze(k_sort_off(:,:,c_idx,t_idx)),x_tot,y_tot);
            k_fft2_off2(:,:,c_idx, t_idx) = fft2(squeeze(k_sort_off2(:,:,c_idx,t_idx)),x_tot,y_tot); 
        end
    end


    wat_k_space_on = sum(wat_k_space_on,3);
    wat_k_space_on2 = sum(wat_k_space_on2,3);
    wat_k_space_off = sum(wat_k_space_off,3);
    wat_k_space_off2 = sum(wat_k_space_off2,3);


    %wat_map = squeeze(abs(fftshift(sum(k_fft2_on_wat,3))));

    % Phase each coil before summing over the channels.
    k_fft2_on_wat = repmat(k_fft2_on_wat, [1 1 1 1024]);
    k_fft2_on2_wat = repmat(k_fft2_on2_wat, [1 1 1 1024]);
    k_fft2_off_wat = repmat(k_fft2_off_wat, [1 1 1 1024]);
    k_fft2_off2_wat = repmat(k_fft2_off2_wat, [1 1 1 1024]);
    
    if isempty(k_fft2_wat_ref)
        k_fft2_on_phased = k_fft2_on.*conj(k_fft2_on_wat)./abs(k_fft2_on_wat);
        k_fft2_on2_phased = k_fft2_on2.*conj(k_fft2_on2_wat)./abs(k_fft2_on2_wat);
        k_fft2_off_phased = k_fft2_off.*conj(k_fft2_off_wat)./abs(k_fft2_off_wat);
        k_fft2_off2_phased = k_fft2_off2.*conj(k_fft2_off2_wat)./abs(k_fft2_off2_wat);
        
        parameters.wat_file = 'none';
    else
        k_fft2_on_phased = k_fft2_on.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
        k_fft2_on2_phased = k_fft2_on2.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
        k_fft2_off_phased = k_fft2_off.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
        k_fft2_off2_phased = k_fft2_off2.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
        
        parameters.wat_file = wat_file;
    end

    k_fft2_on = squeeze(sum(k_fft2_on_phased,3));
    k_fft2_on2 = squeeze(sum(k_fft2_on2_phased,3));
    k_fft2_off = squeeze(sum(k_fft2_off_phased,3));
    k_fft2_off2 = squeeze(sum(k_fft2_off2_phased,3));
else % MEGA multislice, order: slice, kx, ky, coil, data
    
    sz_k_sort_on = size(k_sort_on);
    sz_k_sort_off = size(k_sort_off);

    k_fft2_on = zeros([3,x_tot,y_tot,sz_k_sort_on(4:end)]);
    k_fft2_on2 = zeros([3,x_tot,y_tot,sz_k_sort_on(4:end)]);
    k_fft2_off = zeros([3,x_tot,y_tot,sz_k_sort_off(4:end)]);
    k_fft2_off2 = zeros([3,x_tot,y_tot,sz_k_sort_off(4:end)]);

    sz_k_on = size(k_sort_on);
    k_fft2_on_wat = zeros([3,x_tot,y_tot, sz_k_on(4)]);
    k_fft2_on2_wat = zeros([3,x_tot,y_tot, sz_k_on(4)]);
    k_fft2_off_wat = zeros([3,x_tot,y_tot, sz_k_on(4)]);
    k_fft2_off2_wat = zeros([3,x_tot,y_tot, sz_k_on(4)]);
    
    wat_k_space_on = zeros(sz_k_on(1:4));
    wat_k_space_on2 = zeros(sz_k_on(1:4));
    wat_k_space_off = zeros(sz_k_on(1:4));
    wat_k_space_off2 = zeros(sz_k_on(1:4));

    disp('Taking Fourier transforms.')
    count_ft = 0;
    for c_idx = 1:size(k_sort_on,4) % Each coil
        wat_on_peak = squeeze(k_sort_on(:,:,:,c_idx,1));
        wat_on2_peak = squeeze(k_sort_on2(:,:,:,c_idx,1));
        wat_off_peak = squeeze(k_sort_off(:,:,:,c_idx,1));
        wat_off2_peak = squeeze(k_sort_off2(:,:,:,c_idx,1));
        
        for s_idx = 1:3
            k_fft2_on_wat(s_idx,:,:,c_idx) = fft2(squeeze(wat_on_peak(s_idx,:,:)), x_tot,y_tot);
            k_fft2_on2_wat(s_idx,:,:,c_idx) = fft2(squeeze(wat_on2_peak(s_idx,:,:)), x_tot,y_tot);
            k_fft2_off_wat(s_idx,:,:,c_idx) = fft2(squeeze(wat_off_peak(s_idx,:,:)), x_tot,y_tot);
            k_fft2_off2_wat(s_idx,:,:,c_idx) = fft2(squeeze(wat_off2_peak(s_idx,:,:)), x_tot,y_tot);
        end

        wat_k_space_on(:,:,:,c_idx) = wat_on_peak.*conj(wat_on_peak)./abs(wat_on_peak);
        wat_k_space_on2(:,:,:,c_idx) = wat_on2_peak.*conj(wat_on2_peak)./abs(wat_on2_peak);
        wat_k_space_off(:,:,:,c_idx) = wat_off_peak.*conj(wat_off_peak)./abs(wat_off_peak);
        wat_k_space_off2(:,:,:,c_idx) = wat_off2_peak.*conj(wat_off2_peak)./abs(wat_off2_peak);
        
        for t_idx = 1:size(k_sort_on, 5) % Each point in time
            for s_idx = 1:3
                count_ft = count_ft + 1;
                k_fft2_on(s_idx,:,:,c_idx, t_idx) = fft2(squeeze(k_sort_on(s_idx,:,:,c_idx,t_idx)),x_tot,y_tot); % zerofill k-space
                k_fft2_on2(s_idx,:,:,c_idx, t_idx) = fft2(squeeze(k_sort_on2(s_idx,:,:,c_idx,t_idx)),x_tot,y_tot); 
                k_fft2_off(s_idx,:,:,c_idx, t_idx) = fft2(squeeze(k_sort_off(s_idx,:,:,c_idx,t_idx)),x_tot,y_tot);
                k_fft2_off2(s_idx,:,:,c_idx, t_idx) = fft2(squeeze(k_sort_off2(s_idx,:,:,c_idx,t_idx)),x_tot,y_tot); 
            end
        end
        t2 = toc;
        disp(sprintf('Elapsed time is %0.1f minutes.', t2/60))
    end


    wat_k_space_on = sum(wat_k_space_on,4);
    wat_k_space_on2 = sum(wat_k_space_on2,4);
    wat_k_space_off = sum(wat_k_space_off,4);
    wat_k_space_off2 = sum(wat_k_space_off2,4);

    % Phase each coil before summing over the channels.
    k_fft2_on_wat = repmat(k_fft2_on_wat, [1 1 1 1 1024]);
    k_fft2_on2_wat = repmat(k_fft2_on2_wat, [1 1 1 1 1024]);
    k_fft2_off_wat = repmat(k_fft2_off_wat, [1 1 1 1 1024]);
    k_fft2_off2_wat = repmat(k_fft2_off2_wat, [1 1 1 1 1024]);
    
    if isempty(k_fft2_wat_ref)
        k_fft2_on_phased = k_fft2_on.*conj(k_fft2_on_wat)./abs(k_fft2_on_wat);
        k_fft2_on2_phased = k_fft2_on2.*conj(k_fft2_on2_wat)./abs(k_fft2_on2_wat);
        k_fft2_off_phased = k_fft2_off.*conj(k_fft2_off_wat)./abs(k_fft2_off_wat);
        k_fft2_off2_phased = k_fft2_off2.*conj(k_fft2_off2_wat)./abs(k_fft2_off2_wat);
        
        parameters.wat_file = [];
    else
        k_fft2_on_phased = k_fft2_on.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
        k_fft2_on2_phased = k_fft2_on2.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
        k_fft2_off_phased = k_fft2_off.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
        k_fft2_off2_phased = k_fft2_off2.*conj(k_fft2_wat_ref)./abs(k_fft2_wat_ref);
        
        parameters.wat_file = wat_file;
    end

    k_fft2_on = squeeze(sum(k_fft2_on_phased,4));
    k_fft2_on2 = squeeze(sum(k_fft2_on2_phased,4));
    k_fft2_off = squeeze(sum(k_fft2_off_phased,4));
    k_fft2_off2 = squeeze(sum(k_fft2_off2_phased,4));
end
t2 = toc;
disp(sprintf('Elapsed time is %0.1f minutes.', t2/60))
% % ------------------------------------------------------------------------------------------
% % ------------------- %
% % HLSVD water removal %
% % ------------------- %
% 
% Taken from Gannet 3.0, assuming 2000 Hz spectral width
if water_filter
    if ~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'SE multislice')
        if (strcmp(seq_type, 'HERMES') || strcmp(seq_type, 'HERMES lip sup'))
            k_fft2_on = k_fft2_on1_on2;
        end

        count = 0;
        for x = 1:size(k_fft2_on,1)
            for y = 1:size(k_fft2_on,2)

                count = count + 1;
                disp(sprintf('Removing water of %d of %d.', count, x_tot*y_tot))
                if (strcmp(seq_type, 'HERMES') || strcmp(seq_type, 'HERMES lip sup'))
                    k_fft2_on1_on2(x,y,:) = waterremovalSVD(squeeze(k_fft2_on1_on2(x,y,:)), 2, 8, -0.08, 0.08, 0, 1024);
                    k_fft2_on1_off2(x,y,:) = waterremovalSVD(squeeze(k_fft2_on1_off2(x,y,:)), 2, 8, -0.08, 0.08, 0, 1024);
                    k_fft2_off1_on2(x,y,:) = waterremovalSVD(squeeze(k_fft2_off1_on2(x,y,:)), 2, 8, -0.08, 0.08, 0, 1024);
                    k_fft2_off1_off2(x,y,:) = waterremovalSVD(squeeze(k_fft2_off1_off2(x,y,:)), 2, 8, -0.08, 0.08, 0, 1024);
                else
                    k_fft2_on(x,y,:) = waterremovalSVD(squeeze(k_fft2_on(x,y,:)), 2, 11, -0.11, 0.11, 0, 1024);
                    k_fft2_on2(x,y,:) = waterremovalSVD(squeeze(k_fft2_on2(x,y,:)), 2, 11, -0.11, 0.11, 0, 1024);
                    k_fft2_off(x,y,:) = waterremovalSVD(squeeze(k_fft2_off(x,y,:)), 2, 11, -0.11, 0.11, 0, 1024);
                    k_fft2_off2(x,y,:) = waterremovalSVD(squeeze(k_fft2_off2(x,y,:)), 2, 11, -0.11, 0.11, 0, 1024);  
                end
            end
        end
    else
        count = 0;
        for s_idx = 1:3
            for x = 1:size(k_fft2_on,2)
                for y = 1:size(k_fft2_on,3)

                    count = count + 1;
                    disp(sprintf('Removing water of %d of %d.', count, 3*x_tot*y_tot))

                    k_fft2_on(s_idx,x,y,:) = waterremovalSVD(squeeze(k_fft2_on(s_idx,x,y,:)), 2, 8, -0.08, 0.08, 0, 1024);
                    k_fft2_on2(s_idx,x,y,:) = waterremovalSVD(squeeze(k_fft2_on2(s_idx,x,y,:)), 2, 8, -0.08, 0.08, 0, 1024);
                    k_fft2_off(s_idx,x,y,:) = waterremovalSVD(squeeze(k_fft2_off(s_idx,x,y,:)), 2, 8, -0.08, 0.08, 0, 1024);
                    k_fft2_off2(s_idx,x,y,:) = waterremovalSVD(squeeze(k_fft2_off2(s_idx,x,y,:)), 2, 8, -0.08, 0.08, 0, 1024);
                end
                t2 = toc;
                disp(sprintf('Elapsed time is %0.1f minutes.', t2/60))
            end
        end
    end
        
end
% ------------------------------------------------------------------------------------------


exp_lb = permute(squeeze((repmat(exp(-(lb*pi*(1:1024))/1024), [1 1 x_tot y_tot]))), [2 3 1]);

if (strcmp(seq_type, 'HERMES') || strcmp(seq_type, 'HERMES lip sup'))
    on1_on2_spec = fftshift(fft(k_fft2_on1_on2.*exp_lb,1024*spec_zfill,3));
    on1_off2_spec = fftshift(fft(k_fft2_on1_off2.*exp_lb,1024*spec_zfill,3));
    off1_on2_spec = fftshift(fft(k_fft2_off1_on2.*exp_lb,1024*spec_zfill,3));
    off1_off2_spec = fftshift(fft(k_fft2_off1_off2.*exp_lb,1024*spec_zfill,3));
    cho_im_off = zeros(size(on1_on2_spec,1),size(on1_on2_spec,2));
    
else
    if (~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'Water multislice') && ~strcmp(seq_type, 'SE multislice'))
        on_spec1 = fftshift(fft(k_fft2_on.*exp_lb,1024*spec_zfill,3));
        off_spec1 = fftshift(fft(k_fft2_off.*exp_lb,1024*spec_zfill,3));
        on_spec2 = fftshift(fft(k_fft2_on2.*exp_lb,1024*spec_zfill,3));
        off_spec2 = fftshift(fft(k_fft2_off2.*exp_lb,1024*spec_zfill,3));
        
%         off_spec_nolb_1 = fftshift(fft(k_fft2_off,1024*spec_zfill,3));
%         off_spec_nolb_2 = fftshift(fft(k_fft2_off2,1024*spec_zfill,3));
%         on_spec_nolb_1 = fftshift(fft(k_fft2_off,1024*spec_zfill,3));
%         on_spec_nolb_2 = fftshift(fft(k_fft2_off2,1024*spec_zfill,3));

        % I don't want zerofiling to be happening
        % These shoudl be the corrected spectra, but this is not working
        % for whatever reason
        off_spec_nolb_1 = fftshift(fft(k_fft2_off,1024,3));
        off_spec_nolb_2 = fftshift(fft(k_fft2_off2,1024,3));
        on_spec_nolb_1 = fftshift(fft(k_fft2_off,1024,3));
        on_spec_nolb_2 = fftshift(fft(k_fft2_off2,1024,3));
        
        %Therefore I'm overwriting the recently 'corrected' specctra with
        %the uncorrected ones
        on_spec_nolb_1 = fftshift(fft(squeeze(k_sort_on_phased_k1),1024,3),3);
        off_spec_nolb_1 = fftshift(fft(squeeze(k_sort_off_phased_k1),1024,3),3);
        on_spec_nolb_2 = fftshift(fft(squeeze(k_sort_on_phased_k2),1024,3),3);
        off_spec_nolb_2 = fftshift(fft(squeeze(k_sort_off_phased_k2),1024,3),3);
            
        on_spec_nolb_1(isnan(on_spec_nolb_1)) = 0 + 1i*0;
        off_spec_nolb_1(isnan(off_spec_nolb_1)) = 0 + 1i*0;
        on_spec_nolb_2(isnan(on_spec_nolb_2)) = 0 + 1i*0;
        off_spec_nolb_2(isnan(off_spec_nolb_2)) = 0 + 1i*0;
        
        
        on_spec = (on_spec1 + on_spec2)/2;
        off_spec = (off_spec1 + off_spec2)/2;
        off_spec_no_lb = (off_spec_nolb_1 + off_spec_nolb_2)/2;
        on_spec_no_lb = (on_spec_nolb_1 + on_spec_nolb_2)/2;
    else
        exp_lb = permute(squeeze((repmat(exp(-(lb*pi*(1:1024))/1024), [1 1 3 x_tot y_tot]))), [2 3 4 1]);
        on_spec1 = fftshift(fft(k_fft2_on.*exp_lb,1024*spec_zfill,4));
        off_spec1 = fftshift(fft(k_fft2_off.*exp_lb,1024*spec_zfill,4));
        on_spec2 = fftshift(fft(k_fft2_on2.*exp_lb,1024*spec_zfill,4));
        off_spec2 = fftshift(fft(k_fft2_off2.*exp_lb,1024*spec_zfill,4));
        
        off_spec_nolb_1 = fftshift(fft(k_fft2_off,1024*spec_zfill,4));
        off_spec_nolb_2 = fftshift(fft(k_fft2_off2,1024*spec_zfill,4));
        on_spec_nolb_1 = fftshift(fft(k_fft2_off,1024*spec_zfill,4));
        on_spec_nolb_2 = fftshift(fft(k_fft2_off2,1024*spec_zfill,4));
        
        on_spec_nolb_1(isnan(on_spec_nolb_1)) = 0 + 1i*0;
        off_spec_nolb_1(isnan(off_spec_nolb_1)) = 0 + 1i*0;
        on_spec_nolb_2(isnan(on_spec_nolb_2)) = 0 + 1i*0;
        off_spec_nolb_2(isnan(off_spec_nolb_2)) = 0 + 1i*0;

        on_spec = (on_spec1 + on_spec2)/2;
        off_spec = (off_spec1 + off_spec2)/2;
        off_spec_no_lb = (off_spec_nolb_1 + off_spec_nolb_2)/2;
        on_spec_no_lb = (on_spec_nolb_1 + on_spec_nolb_2)/2;
        
    end    
end

end