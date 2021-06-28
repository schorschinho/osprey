function [k_fft2_wat_ref, k_fft2_wat_ref_no_k_zfill, k_sort_b4_2dfft] = process_wat_ref(this_file, k_zfill, seq_type)

    fprintf('\nReading water scan parameters.')
    [data] = loadRawKspace(this_file);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine dimensions of the acquisition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Determine number of channels
    n_coils = length(unique(data.chan));
    
    % Determine number of mixes
    n_mixes = length(unique(data.mix));
    
    % Determine number of averages per mix
    n_averages = data.kspace_properties.number_of_signal_averages;

    % Determine number of data points per scan
    n_points = data.kspace_properties.F_resolution(1);
     % Determine number of data points per scan
    kz_tot = data.kspace_properties.number_of_locations(1);
    % Determine number of data points per scan
    kx_tot = abs(data.kspace_properties.kx_range(1)) + abs(data.kspace_properties.kx_range(2)) +1 ;
     % Determine number of data points per scan
    ky_tot = abs(data.kspace_properties.ky_range(1)) + abs(data.kspace_properties.ky_range(2)) +1 ;
    
    if kz_tot > 1
        seq_type = 'MEGA multislice';
    end
    
    data.kx = data.kx + abs(min(data.kx)) + 1;
    data.ky = data.ky + abs(min(data.ky)) + 1;
    data.loca = data.loca + abs(min(data.kz)) + 1;
    data.aver = data.aver + 1;
    data.chan = data.chan + 1;
    
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start splitting the list of total scans into its parts:
    % Noise scans, water-suppressed scans, and water-unsuppressed scans.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    noise_line = sum(strcmp(data.typ,'NOI'));
    
    % Separate the water-suppressed from the water-unsuppressed scans.
    % Water-suppressed scans have the data type 'STD' and mix index 0:
    isdata = strcmp(data.typ,'STD') & (data.mix == 0);
    data_matrix = cell2mat(data.complexdata(isdata));
    
if kz_tot > 1
    seq_type = 'MEGA multislice';
end
  

if strcmp(seq_type, 'MEGA multislice') || strcmp(seq_type, 'SE multislice')
    k_sort = zeros(kz_tot, kx_tot, ky_tot, n_coils, n_points);
else
    k_sort = zeros(kx_tot, ky_tot, n_coils, n_points);
end

disp('Reorganizing water reference k-space locations.')
% Rearrange k-space values (so no negative indices)

        
for dl = 1:size(data_matrix,1)    
    if ~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'SE multislice')
        k_sort(data.kx(dl+noise_line), data.ky(dl+noise_line), data.chan(dl+noise_line), :) = data_matrix(dl,:);
    else
        k_sort(data.loca(dl+noise_line),data.kx(dl+noise_line), data.ky(dl+noise_line), data.chan(dl+noise_line), :) = data_matrix(dl,:);
    end
    
end

% ------------------------------------------------------------------------------------------
% --------------- %
% Hanning Filter  %
% --------------- %
% Apply a Hanning filter on k-space data to filter out noise that are
% at higher k-space values

if ~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'SE multislice')
    hanning_x = repmat(hanning(kx_tot), [1 ky_tot n_coils n_points]);
    hanning_y = permute(repmat(hanning(ky_tot), [1 kx_tot n_coils n_points]), [2 1 3 4]);
else
    hanning_x = repmat(hanning(kx_tot), [1 kz_tot ky_tot n_coils n_points]);
    hanning_y = permute(repmat(hanning(ky_tot), [1 kx_tot kz_tot n_coils n_points]), [3 2 1 4 5]);
    hanning_x = permute(hanning_x, [2 1 3 4 5]);
end

k_sort = k_sort.*hanning_x;
k_sort = k_sort.*hanning_y;
k_sort_b4_2dfft = k_sort;

% ------------------------------------------------------------------------------------------
% For each point in time and coil take the 2D fft
if ~strcmp(seq_type, 'MEGA multislice') && ~strcmp(seq_type, 'SE multislice')
    sz_k_sort = size(k_sort);
    k_fft2_wat_ref = zeros([kx_tot,ky_tot,sz_k_sort(3:end)]);
    k_fft2_wat = zeros([kx_tot,ky_tot, sz_k_sort(3)]);
    k_fft2_wat_no_k_zfill = zeros([kx_tot,ky_tot, sz_k_sort(3)]);
    k_fft2_wat_ref_no_k_zfill = zeros([kx_tot,ky_tot, sz_k_sort(3)]);

    disp('Taking Fourier transforms of water reference data.')
    for c_idx = 1:size(k_sort,3) % Each coil
        wat_on_peak = squeeze(k_sort(:,:,c_idx,1));
        k_fft2_wat(:,:,c_idx) = fft2(wat_on_peak, kx_tot,ky_tot);
        k_fft2_wat_no_k_zfill(:,:,c_idx) = fft2(wat_on_peak);

        for t_idx = 1:size(k_sort, 4) % Each point in time
            k_fft2_wat_ref(:,:,c_idx, t_idx) = fft2(squeeze(k_sort(:,:,c_idx,t_idx)),kx_tot,ky_tot);
            k_fft2_wat_ref_no_k_zfill(:,:,c_idx, t_idx) = fft2(squeeze(k_sort(:,:,c_idx,t_idx)));
        end
    end
else
    sz_k_sort = size(k_sort);
    k_fft2_wat_ref = zeros([3,kx_tot,ky_tot,sz_k_sort(4:end)]);
    k_fft2_wat = zeros([3,kx_tot,ky_tot, sz_k_sort(4)]);
    k_fft2_wat_no_k_zfill = zeros([3,kx_tot,ky_tot, sz_k_sort(4)]);
    k_fft2_wat_ref_no_k_zfill = zeros([3,kx_tot,ky_tot, sz_k_sort(4)]);

    disp('Taking Fourier transforms of water reference data.')
    for c_idx = 1:size(k_sort,4) % Each coil
        wat_on_peak = squeeze(k_sort(:,:,:,c_idx,1));

        for s_idx = 1:3
            k_fft2_wat(s_idx,:,:,c_idx) = fft2(squeeze(wat_on_peak(s_idx,:,:)), kx_tot,ky_tot);
            k_fft2_wat_no_k_zfill(s_idx,:,:,c_idx) = fft2(squeeze(wat_on_peak(s_idx,:,:)));
        end

        for t_idx = 1:size(k_sort, 5) % Each point in time
            for s_idx = 1:3
                k_fft2_wat_ref(s_idx,:,:,c_idx, t_idx) = fft2(squeeze(k_sort(s_idx,:,:,c_idx,t_idx)),kx_tot,ky_tot);
                k_fft2_wat_ref_no_k_zfill(s_idx,:,:,c_idx, t_idx) = fft2(squeeze(k_sort(s_idx,:,:,c_idx,t_idx)));
            end
        end
    end
end
