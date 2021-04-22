function [k_fft2_wat_ref, k_fft2_wat_ref_no_k_zfill, k_sort_b4_2dfft] = process_wat_ref(this_file, k_zfill, seq_type)

disp('Opening water reference data.')

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


coil = zeros(1,length(data_lines));
kx = coil;
ky = coil;
avg = coil;
sign = coil;
loc = coil;

for dl = 1:(length(data_lines))
    if dl == length(data_lines)
        this_dl = cellfun(@str2num,{scan_params{(data_lines(dl)+1):(data_lines(dl)+20)}});
    else
        this_dl = cellfun(@str2num,{scan_params{(data_lines(dl)+1):(data_lines(dl + 1)-1)}});
    end
    loc(dl) = this_dl(5) + 1;
    coil(dl) = this_dl(6);
%     kx(dl) = this_dl(9);
%     ky(dl) = this_dl(10);
    kx(dl) = this_dl(10);
    ky(dl) = this_dl(9);
    avg(dl) = this_dl(12);
    sign(dl) = this_dl(13);
end

kx_tot = max(kx) - min(kx) + 1;
ky_tot = max(ky) - min(ky) + 1;

x_tot = k_zfill*kx_tot;
y_tot = k_zfill*ky_tot;

if strcmp(seq_type, 'MEGA multislice') || strcmp(seq_type, 'SE multislice')
    k_sort = zeros(max(loc), kx_tot, ky_tot, max(coil) + 1, 1024);
else
    k_sort = zeros(kx_tot, ky_tot, max(coil) + 1, 1024);
end

disp('Reorganizing water reference k-space locations.')
% Rearrange k-space values (so no negative indices)
for dl = 1:(length(data_lines) - 1)
    this_kx = kx(dl) + abs(min(kx)) + 1;
    this_ky = ky(dl) + abs(min(ky)) + 1;
    this_avg = avg(dl);
    this_coil = coil(dl);
    this_loc = loc(dl);
    
    if ~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'SE multislice')
        k_sort(this_kx, this_ky, this_coil + 1, :) = data(:,dl);
    else
        k_sort(this_loc,this_kx, this_ky, this_coil + 1, :) = data(:,dl);
    end
    
end

% ------------------------------------------------------------------------------------------
% --------------- %
% Hanning Filter  %
% --------------- %
% Apply a Hanning filter on k-space data to filter out noise that are
% at higher k-space values

if ~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'SE multislice')
    hanning_x = repmat(hanning(kx_tot), [1 ky_tot max(coil) + 1 1024]);
    hanning_y = permute(repmat(hanning(ky_tot), [1 kx_tot max(coil) + 1 1024]), [2 1 3 4]);
else
    hanning_x = repmat(hanning(kx_tot), [1 max(loc) ky_tot max(coil) + 1 1024]);
    hanning_y = permute(repmat(hanning(ky_tot), [1 kx_tot max(loc) max(coil) + 1 1024]), [3 2 1 4 5]);
    hanning_x = permute(hanning_x, [2 1 3 4 5]);
end

k_sort = k_sort.*hanning_x;
k_sort = k_sort.*hanning_y;
k_sort_b4_2dfft = k_sort;

% ------------------------------------------------------------------------------------------
% For each point in time and coil take the 2D fft
if ~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'SE multislice')
    sz_k_sort = size(k_sort);
    k_fft2_wat_ref = zeros([x_tot,y_tot,sz_k_sort(3:end)]);
    k_fft2_wat = zeros([x_tot,y_tot, sz_k_sort(3)]);
    k_fft2_wat_no_k_zfill = zeros([kx_tot,ky_tot, sz_k_sort(3)]);
    k_fft2_wat_ref_no_k_zfill = zeros([kx_tot,ky_tot, sz_k_sort(3)]);

    disp('Taking Fourier transforms of water reference data.')
    for c_idx = 1:size(k_sort,3) % Each coil
        wat_on_peak = squeeze(k_sort(:,:,c_idx,1));
        k_fft2_wat(:,:,c_idx) = fft2(wat_on_peak, x_tot,y_tot);
        k_fft2_wat_no_k_zfill(:,:,c_idx) = fft2(wat_on_peak);

        for t_idx = 1:size(k_sort, 4) % Each point in time
            k_fft2_wat_ref(:,:,c_idx, t_idx) = fft2(squeeze(k_sort(:,:,c_idx,t_idx)),x_tot,y_tot);
            k_fft2_wat_ref_no_k_zfill(:,:,c_idx, t_idx) = fft2(squeeze(k_sort(:,:,c_idx,t_idx)));
        end
    end
else
    sz_k_sort = size(k_sort);
    k_fft2_wat_ref = zeros([3,x_tot,y_tot,sz_k_sort(4:end)]);
    k_fft2_wat = zeros([3,x_tot,y_tot, sz_k_sort(4)]);
    k_fft2_wat_no_k_zfill = zeros([3,kx_tot,ky_tot, sz_k_sort(4)]);
    k_fft2_wat_ref_no_k_zfill = zeros([3,kx_tot,ky_tot, sz_k_sort(4)]);

    disp('Taking Fourier transforms of water reference data.')
    for c_idx = 1:size(k_sort,4) % Each coil
        wat_on_peak = squeeze(k_sort(:,:,:,c_idx,1));

        for s_idx = 1:3
            k_fft2_wat(s_idx,:,:,c_idx) = fft2(squeeze(wat_on_peak(s_idx,:,:)), x_tot,y_tot);
            k_fft2_wat_no_k_zfill(s_idx,:,:,c_idx) = fft2(squeeze(wat_on_peak(s_idx,:,:)));
        end

        for t_idx = 1:size(k_sort, 5) % Each point in time
            for s_idx = 1:3
                k_fft2_wat_ref(s_idx,:,:,c_idx, t_idx) = fft2(squeeze(k_sort(s_idx,:,:,c_idx,t_idx)),x_tot,y_tot);
                k_fft2_wat_ref_no_k_zfill(s_idx,:,:,c_idx, t_idx) = fft2(squeeze(k_sort(s_idx,:,:,c_idx,t_idx)));
            end
        end
    end
end
