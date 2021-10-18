function [k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
    on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,...
    k_ph_corr_on1, k_ph_corr_on2, k_ph_corr_off1, k_ph_corr_off2,...
    on1_replace_track, off1_replace_track, on2_replace_track, ...
    off2_replace_track, zero_replace_track, k_space_locs, corr_options] = motion_correct_mrsi(k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
                                                                                                on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,...
                                                                                                spec_zfill, seq_type, kx_tot, ky_tot,thresholds)
% Written by Kimberly Chan
% Johns Hopkins School of Medicine
% Last updated 03/28/19


% -------------- %
% --Thresholds-- %
% -------------- %
thresh = thresholds.thresh;
ph_thresh = thresholds.ph_thresh;
last_resort_thresh = thresholds.last_resort_thresh;

k_ph_corr_on1 = zeros(size(k_sort_on2,1),size(k_sort_on2,2),size(k_sort_on2,3));
k_ph_corr_on2 = zeros(size(k_sort_on2,1),size(k_sort_on2,2),size(k_sort_on2,3));
k_ph_corr_off1 = zeros(size(k_sort_on2,1),size(k_sort_on2,2),size(k_sort_on2,3));
k_ph_corr_off2 = zeros(size(k_sort_on2,1),size(k_sort_on2,2),size(k_sort_on2,3));

% Keep track of which averages are replaced or zeroed
on1_replace_track = zeros(size(k_sort_on2,1),size(k_sort_on2,2),size(k_sort_on2,3));
off1_replace_track = zeros(size(k_sort_on2,1),size(k_sort_on2,2),size(k_sort_on2,3));
on2_replace_track = zeros(size(k_sort_on2,1),size(k_sort_on2,2),size(k_sort_on2,3));
off2_replace_track = zeros(size(k_sort_on2,1),size(k_sort_on2,2),size(k_sort_on2,3));

zero_replace_track = zeros(size(k_sort_on2,1),size(k_sort_on2,2),size(k_sort_on2,3));

k_space_locs = zeros(size(k_sort_on2,1),size(k_sort_on2,2),size(k_sort_on2,3)); 

% Keep track of which options get triggered:
corr_options = [];

dims = size(k_sort_on);
n_points = dims(end);
    
% Eliminate bad averages.
count_k = 0;
track_all_p_thresh = [];
disp('Eliminating bad k-space points')
for this_kz = 1:size(on_spec_k_1,1)
    for this_kx = 1:size(on_spec_k_1,2)
        for this_ky = 1:size(on_spec_k_1,3)
            count_k = count_k + 1;
            %disp(sprintf('Eliminating bad k-space points of %d of %d.', count_k, 3*size(on_spec_k_1,2)*size(on_spec_k_1,3)))
            ppm_axis = linspace(-1000,1000,n_points*spec_zfill)/128 + 4.68;
            if this_kx == 7 && this_ky == 6
                this_kx;
                this_ky;
            end
            % Water approximate range
            met_start_wat = find(round(ppm_axis*100)/100 >= 4.4);
            met_end_wat = find(round(ppm_axis*100)/100 <= 5.2);
            if length(met_start_wat) > 1
                met_start_wat = met_start_wat(1);
            end
            if length(met_end_wat) > 1
                met_end_wat = met_end_wat(end);
            end

            % Lipid approximate range
            met_start_lip = find(round(ppm_axis*100)/100 >= 0);
            met_end_lip = find(round(ppm_axis*100)/100 <= 2.5);
            if length(met_start_lip) > 1
                met_start_lip = met_start_lip(1);
            end
            if length(met_end_lip) > 1
                met_end_lip = met_end_lip(end);
            end


            % Cr Cho approximate range
            met_start_crcho = find(round(ppm_axis*100)/100 >= 2.8);
            met_end_crcho = find(round(ppm_axis*100)/100 <= 3.5);
            if length(met_start_crcho) > 1
                met_start_crcho = met_start_crcho(1);
            end
            if length(met_end_crcho) > 1
                met_end_crcho = met_end_crcho(end);
            end


            % Reduce ON/OFF spectra into water, lipid, and Cr/Cho
            % ranges.
            this_on1_wat = real(squeeze(on_spec_k_1( this_kz, this_kx, this_ky, met_start_wat:met_end_wat)));
            this_off1_wat = real(squeeze(off_spec_k_1( this_kz, this_kx, this_ky, met_start_wat:met_end_wat)));
            this_on2_wat = real(squeeze(on_spec_k_2( this_kz, this_kx, this_ky, met_start_wat:met_end_wat)));
            this_off2_wat = real(squeeze(off_spec_k_2( this_kz, this_kx, this_ky, met_start_wat:met_end_wat)));

            this_on1_wat_im = squeeze(on_spec_k_1( this_kz, this_kx, this_ky, met_start_wat:met_end_wat));
            this_off1_wat_im = squeeze(off_spec_k_1( this_kz, this_kx, this_ky, met_start_wat:met_end_wat));
            this_on2_wat_im = squeeze(on_spec_k_2( this_kz, this_kx, this_ky, met_start_wat:met_end_wat));
            this_off2_wat_im = squeeze(off_spec_k_2( this_kz, this_kx, this_ky, met_start_wat:met_end_wat));

            this_on1_lip = real(squeeze(on_spec_k_1( this_kz, this_kx, this_ky, met_start_lip:met_end_lip)));
            this_off1_lip = real(squeeze(off_spec_k_1( this_kz, this_kx, this_ky, met_start_lip:met_end_lip)));
            this_on2_lip = real(squeeze(on_spec_k_2( this_kz, this_kx, this_ky, met_start_lip:met_end_lip)));
            this_off2_lip = real(squeeze(off_spec_k_2( this_kz, this_kx, this_ky, met_start_lip:met_end_lip)));

            this_on1_lip_im = squeeze(on_spec_k_1( this_kz, this_kx, this_ky, met_start_lip:met_end_lip));
            this_off1_lip_im = squeeze(off_spec_k_1( this_kz, this_kx, this_ky, met_start_lip:met_end_lip));
            this_on2_lip_im = squeeze(on_spec_k_2( this_kz, this_kx, this_ky, met_start_lip:met_end_lip));
            this_off2_lip_im = squeeze(off_spec_k_2( this_kz, this_kx, this_ky, met_start_lip:met_end_lip));

            this_on1_crcho = real(squeeze(on_spec_k_1( this_kz, this_kx, this_ky, met_start_crcho:met_end_crcho)));
            this_off1_crcho = real(squeeze(off_spec_k_1( this_kz, this_kx, this_ky, met_start_crcho:met_end_crcho)));
            this_on2_crcho = real(squeeze(on_spec_k_2( this_kz, this_kx, this_ky, met_start_crcho:met_end_crcho)));
            this_off2_crcho = real(squeeze(off_spec_k_2( this_kz, this_kx, this_ky, met_start_crcho:met_end_crcho)));

            if all(this_on1_wat)% && all(sum([this_kx this_ky] ~= [round(kx_tot/2) round(ky_tot/2)]))
                % Find all the pairs of correlation coefficients
                % (pearson's)
                p_on1_on2_wat = corrcoef(this_on1_wat,this_on2_wat);
                p_on1_on2_wat = p_on1_on2_wat(2);
                p_off1_off2_wat = corrcoef(this_off1_wat,this_off2_wat);
                p_off1_off2_wat = p_off1_off2_wat(2);
                p_on1_off1_wat = corrcoef(this_on1_wat,this_off1_wat);
                p_on1_off1_wat = p_on1_off1_wat(2);
                p_on2_off2_wat = corrcoef(this_on2_wat,this_off2_wat);
                p_on2_off2_wat = p_on2_off2_wat(2);
                p_on1_off2_wat = corrcoef(this_on1_wat,this_off2_wat);
                p_on1_off2_wat = p_on1_off2_wat(2);
                p_on2_off1_wat = corrcoef(this_on2_wat,this_off1_wat);
                p_on2_off1_wat = p_on2_off1_wat(2);


                p_on1_on2_lip = corrcoef(this_on1_lip,this_on2_lip);
                p_on1_on2_lip = p_on1_on2_lip(2);
                p_off1_off2_lip = corrcoef(this_off1_lip,this_off2_lip);
                p_off1_off2_lip = p_off1_off2_lip(2);


                p_on1_on2_crcho = corrcoef(this_on1_crcho,this_on2_crcho);
                p_on1_on2_crcho = p_on1_on2_crcho(2);
                p_off1_off2_crcho = corrcoef(this_off1_crcho,this_off2_crcho);
                p_off1_off2_crcho = p_off1_off2_crcho(2);
                p_on1_off1_crcho = corrcoef(this_on1_crcho,this_off1_crcho);
                p_on1_off1_crcho = p_on1_off1_crcho(2);
                p_on2_off2_crcho = corrcoef(this_on2_crcho,this_off2_crcho);
                p_on2_off2_crcho = p_on2_off2_crcho(2);
                p_on1_off2_crcho = corrcoef(this_on1_crcho,this_off2_crcho);
                p_on1_off2_crcho = p_on1_off2_crcho(2);
                p_on2_off1_crcho = corrcoef(this_on2_crcho,this_off1_crcho);
                p_on2_off1_crcho = p_on2_off1_crcho(2);

                % For water, only compare between ONs or just between
                % OFFs because water is saturated in GSH.
                p_wat = [p_on1_on2_wat p_on1_off1_wat p_on1_off2_wat p_on2_off1_wat p_on2_off2_wat p_off1_off2_wat];
                p_wat(isnan(p_wat)) = 0;
                all_p_thresh_wat = p_wat < thresh;
                p_lip = [p_on1_on2_lip p_off1_off2_lip]; %p_on1_off1_lip p_on2_off2_lip
                p_lip(isnan(p_lip)) = 0;
                all_p_thresh_lip = p_lip < thresh;
                p_crcho = [p_on1_on2_crcho p_on1_off1_crcho p_on1_off2_crcho p_on2_off1_crcho p_on2_off2_crcho p_off1_off2_crcho];
                p_crcho(isnan(p_crcho)) = 0;

                all_p = [mean([p_lip(1) p_wat(1)]), mean([p_wat(2)]),...
                            mean([p_wat(3)]), mean([p_wat(4)]),... % p_crcho(4)
                            mean([p_wat(5)]), mean([p_lip(2) p_wat(6)])];

                all_p_thresh = [all_p_thresh_lip(1) | all_p_thresh_wat(1),...% | all_p_thresh_crcho(1),...
                                all_p_thresh_wat(2:5),...% | all_p_thresh_crcho(2:5),...
                                all_p_thresh_lip(2) | all_p_thresh_wat(6)];% | all_p_thresh_crcho(6)];% | all_p_thresh_lip;

                % Last resort, all of them don't match, unless one pair is > 0.6          
                all_p_thresh_wat_zero = p_wat < last_resort_thresh;
                all_p_thresh_lip_zero = p_lip < last_resort_thresh;
                all_p_thresh_crcho_zero = p_crcho < last_resort_thresh;

                all_p_thresh_zero = [all_p_thresh_lip_zero(1) | all_p_thresh_wat_zero(1),...% | all_p_thresh_crcho_zero(1),...
                                all_p_thresh_wat_zero(2:5),...% | all_p_thresh_crcho_zero(2:5),...
                                all_p_thresh_lip_zero(2) | all_p_thresh_wat_zero(6)];% | all_p_thresh_crcho_zero(6)];% | all_p_thresh_lip;

                track_all_p_thresh = [track_all_p_thresh; all_p_thresh];

                if (this_kx == 6) && (this_ky == 6)
                    xes = 3;
                end

                all_spec = [squeeze(on_spec_k_1( this_kz, this_kx, this_ky, :)), squeeze(on_spec_k_2( this_kz, this_kx, this_ky, :)),...
                            squeeze(off_spec_k_1( this_kz, this_kx, this_ky, :)), squeeze(off_spec_k_2( this_kz, this_kx, this_ky, :))];

                all_spec_wat = permute([this_on1_wat_im, this_on2_wat_im, this_off1_wat_im, this_off2_wat_im], [2 1]);
                all_spec_lip = permute([this_on1_lip_im, this_on2_lip_im, this_off1_lip_im, this_off2_lip_im], [2 1]);
                % Set to zero if all are different, ON1 and ON2 are
                % different, OFF1 and OFF2 are different.
                if all(all_p_thresh_zero) % All sub-acquisitions are different from each other
                    % If all the sub-acquisitions are different from
                    % one another, just replace with zeros (for now,
                    % interpolation in the future?)
                    k_sort_on(this_kz, this_kx,this_ky,:,:) = zeros([size(k_sort_on,4), size(k_sort_on,5)]);
                    k_sort_on2(this_kz, this_kx,this_ky,:,:) = zeros([size(k_sort_on,4), size(k_sort_on,5)]);
                    k_sort_off(this_kz, this_kx,this_ky,:,:) = zeros([size(k_sort_on,4), size(k_sort_on,5)]);
                    k_sort_off2(this_kz, this_kx,this_ky,:,:) = zeros([size(k_sort_on,4), size(k_sort_on,5)]);

                    on_spec_k_1(this_kz, this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,5));
                    on_spec_k_2(this_kz, this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,5));
                    off_spec_k_1(this_kz, this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,5));
                    off_spec_k_2(this_kz, this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,5));

                    zero_replace_track(this_kz, this_kx,this_ky) = 1;
                else
                    if all(all_p_thresh == [1 0 0 0 0 0]) % [1] Phase or replace ON1 OR ON2
                        this_comp = [sum(all_p(2:3)), sum(all_p(4:5))];
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 1];

                        if keep == 1 %Phase ON2 or replace ON2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 0], all_spec_wat, all_spec_lip, all_spec);

                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Phase ON1 or replace ON1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 0], all_spec_wat, all_spec_lip, all_spec);

                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 1 0 0 0 0]) % [2] Phase ON1 or replace ON1 w/ON2 %K L C
                        [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 2];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                            on1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [1 1 1 0 0 0]) % [3] Phase ON1 or replace ON1 w/ON2
                        [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 3];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                            on1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [1 1 1 1 0 0]) % [4] Replace ON1 w/ON2, and replace OFF1 with OFF2
                        [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 4];

                        if ph_corr_coeff(1) > ph_thresh
                            k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                            on1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end

                        if ph_corr_coeff(2) > ph_thresh
                            k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                        else
                            k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                            off1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end

                    elseif all(all_p_thresh == [1 1 1 1 1 0]) % [5] Phase or replace ON1 or ON2
                        this_comp = [sum(all_p(2:3)), sum(all_p(4:5))];
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 5];

                        if keep == 1 %Phase or replace ON2 with ON1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Phase or replace ON1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 1 1 1 1 1]) % [6] Pick the best ON/OFF pair
                        this_comp = all_p(2:5);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 6];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                        elseif keep == 2 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                        elseif keep == 3 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 4 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 1 0 0 0 0]) % [7]. Replace ON1 or OFF1
                        this_comp = [sum(all_p([1 3])), sum(all_p([4 6]))];
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 7];

                        if keep == 1 %Replace OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 1 1 0 0 0]) % [8] Replace ON1 w/ON2
                        [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 8];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                            on1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end

                    elseif all(all_p_thresh == [0 1 1 1 0 0]) % [9] Replace ON1 w/ON2
                        [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 9];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                            on1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [0 1 1 1 1 0]) % [10] Pick the best ON/OFF pair
                        this_comp = all_p(2:5);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 10];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                        elseif keep == 2 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                        elseif keep == 3 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 4 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 1 1 1 1 1]) % [11] Pick the best ON/OFF pair
                        this_comp = all_p(2:5);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 11];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                        elseif keep == 2 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                        elseif keep == 3 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 4 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 0 1 0 0 0]) % [12] Replace ON1 w/ON2
                        [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 12];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                            on1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [1 0 0 1 0 0]) % [13] Replace ON2 w/ON1
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 13];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                            on2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [1 0 0 0 1 0]) % [14] Replace ON2 w/ON1
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 14];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                            on2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [1 0 0 0 0 1]) % [15] Pick the best ON/OFF pair
                        this_comp = all_p(2:5);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 15];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                        elseif keep == 2 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                        elseif keep == 3 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 4 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 1 0 1 0 0]) % [16] Find max of 3, 5
                        this_comp = all_p([3 5]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 16];

                        if keep == 1 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 1 0 0 1 0]) % [17] Choose best ON/OFF pair other than #2,5
                        this_comp = all_p(3:4);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 17];

                        if keep == 1 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 1 0 0 0 1]) % [18] Replace OFF1
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 18];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                            off1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end

                    elseif all(all_p_thresh == [0 0 1 1 0 0]) % [19] Pick the best ON/OFF pair except #3,4
                        this_comp = all_p([2 5]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 19];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 0 1 0 1 0]) % [20] Replace OFF2
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 20];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                            off2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [0 0 1 0 0 1]) % [21] Replace OFF2
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 21];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                            off2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [0 0 0 1 1 0]) % [22] Replace ON2 w/ON1
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 22];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                            on1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [0 0 0 1 0 1]) % [23] Replace OFF1
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 23];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                            off1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [0 0 0 0 1 1]) % [24] Replace OFF2
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 24];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                            off2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [0 0 1 0 0 0]) % [25] Replace ON1 or OFF2
                        this_comp = [sum(all_p([1 2])), sum(all_p([5 6]))];
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 25];

                        if keep == 1 %Replace OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 0 0 1 0 0]) % [26] Replace ON2 or OFF1
                        this_comp = [sum(all_p([1 5])), sum(all_p([2 6]))];
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 26];

                        if keep == 1 %Replace OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 0 0 0 1 0]) % [27] Replace ON1 or OFF2
                        this_comp = [sum(all_p([1 4])), sum(all_p([3 6]))];
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 27];

                        if keep == 1 %Replace OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 0 0 0 0 1]) % [28] Replace OFF1 or OFF2
                        this_comp = [sum(all_p([2 4])), sum(all_p([3 5]))];
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 28];

                        if keep == 1 %Replace OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 1 0 1 0 0]) % [29] Replace ON1 or ON2, OFF1
                        this_comp = all_p([3 5]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 29];

                        if keep == 1 %Replace ON2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                        % Replace OFF1
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                            off1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [1 1 0 0 1 0]) % [30] Replace ON1 or ON2, replace OFF1 or OFF2
                        this_comp = all_p([3 4]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 30];

                        if keep == 1 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 1 0 0 0 1]) % [31] Replace ON2 and OFF1, etc. keep good pair
                        this_comp = all_p(3:5);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 31];

                        if keep == 1 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 3 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 0 1 1 0 0]) % [32] Replace ON2 and OFF2, replace ON1 and OFF1
                        this_comp = all_p([2 5]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 32];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end 
                        elseif keep == 2 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end

                    elseif all(all_p_thresh == [1 0 1 0 1 0]) % [33] Replace ON1 or ON2, replace OFF1 or OFF2
                        this_comp = all_p([2 4]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 33];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 0 1 0 0 1]) % [34] Keep good pair
                        this_comp = all_p([2 4 5]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 34];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 3 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 0 0 1 1 0]) % [35] Replace ON2s
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 35];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                            on2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [1 0 0 1 0 1]) % [36] Pick best pair
                        this_comp = all_p([2 3 5]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 36];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 3 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 0 0 0 1 1]) % [37] Pick the best ON/OFF pair
                        this_comp = all_p(2:4);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 37];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 3 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 1 1 0 1 0]) % [38] Replace ON1, OFF2
                        [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 38];

                        if ph_corr_coeff(1) > ph_thresh
                            k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                            on1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end

                        if ph_corr_coeff(2) > ph_thresh
                            k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                        else
                            k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                            off2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [0 1 1 0 0 1]) % [39] Pick best ON, OFF pair
                        this_comp = all_p([4 5]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 39];

                        if keep == 1 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 1 0 1 1 0]) % [40] Replace ON2, OFF1
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 40];

                        if ph_corr_coeff(1) > ph_thresh
                            k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                            on2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end

                        if ph_corr_coeff(2) > ph_thresh
                            k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                        else
                            k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                            off1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [0 1 0 1 0 1]) % [41] Pick the best pair
                        this_comp = all_p([3 5]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 41];

                        if keep == 1 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 1 0 0 1 1]) % [42] Pick the best ON/OFF pair
                        this_comp = all_p(3:4);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 42];

                        if keep == 1 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 0 1 1 1 0]) % [43] Replace ON2, OFF2
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                        if ph_corr_coeff(1) > ph_thresh
                            k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                            on2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end

                        if ph_corr_coeff(2) > ph_thresh
                            k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                        else
                            k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                            off2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [0 0 1 1 0 1]) % [44] Pick the best ON/OFF pair
                        this_comp = all_p([2 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [0 0 1 0 1 1]) % [45] Phase or replace OFF2
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 45];

                        if ph_corr_coeff > ph_thresh
                            k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                            off2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [0 0 0 1 1 1]) % [46] Pick the best ON/OFF pair
                        this_comp = all_p([2 3]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 46];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 1 1 0 1 0]) % [47] Replace ON1, OFF2
                        [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 47];

                        if ph_corr_coeff(1) > ph_thresh
                            k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                            on1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end

                        if ph_corr_coeff(2) > ph_thresh
                            k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                        else
                            k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                            off2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [1 1 1 0 0 1]) % [48] Pick best ON, OFF pair
                        this_comp = all_p([4 5]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 48];

                        if keep == 1 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 1 0 1 1 0]) % [49] Replace ON2, OFF1
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 49];

                        if ph_corr_coeff(1) > ph_thresh
                            k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                            on2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end

                        if ph_corr_coeff(2) > ph_thresh
                            k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                        else
                            k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                            off1_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [1 1 0 1 0 1]) % [50] Pick the best pair
                        this_comp = all_p([3 5]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 50];

                        if keep == 1 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 1 0 0 1 1]) % [51] Replace ON1 or ON2, replace OFF1 or OFF2
                        this_comp = all_p([3 4]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 51];

                        if keep == 1 %Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end 
                        end
                    elseif all(all_p_thresh == [1 0 1 1 1 0]) % [52] Replace ON2, OFF2
                        [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                        corr_options = [corr_options, 52];

                        if ph_corr_coeff(1) > ph_thresh
                            k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                        else
                            k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                            on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                            on2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end

                        if ph_corr_coeff(2) > ph_thresh
                            k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                        else
                            k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                            off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                            off2_replace_track(this_kz, this_kx,this_ky) = 1;
                        end
                    elseif all(all_p_thresh == [1 0 1 1 0 1]) % [53] Pick best pair
                        this_comp = all_p([2 5]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 53];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 0 1 0 1 1]) % [54] Pick the best ON/OFF pair
                        this_comp = all_p([2 4]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 54];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end
                    elseif all(all_p_thresh == [1 0 0 1 1 1]) % [55] Pick the best ON/OFF pair
                        this_comp = all_p([2 3]);
                        keep = find(this_comp == max(this_comp));
                        corr_options = [corr_options, 55];

                        if keep == 1 %Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        elseif keep == 2 %Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                        end

                    elseif all(all_p_thresh == [0 1 1 1 0 1]) % [56] Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            corr_options = [corr_options, 56];

                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                    elseif all(all_p_thresh == [0 1 1 0 1 1]) % [57] Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            corr_options = [corr_options, 57];

                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                    elseif all(all_p_thresh == [0 1 0 1 1 1]) % [58] Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            corr_options = [corr_options, 58];

                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                    elseif all(all_p_thresh == [0 0 1 1 1 1]) % [59] Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            corr_options = [corr_options, 59];

                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                    elseif all(all_p_thresh == [1 0 1 1 1 1]) % [60] Replace ON2, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 0 1], all_spec_wat, all_spec_lip, all_spec);
                            corr_options = [corr_options, 60];

                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off2(this_kz, this_kx,this_ky,:,:) = k_sort_off(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_1(this_kz, this_kx,this_ky,:,:); 
                                off2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                    elseif all(all_p_thresh == [1 1 0 1 1 1]) % [61] Replace ON2, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([0 1 1 0], all_spec_wat, all_spec_lip, all_spec);
                            corr_options = [corr_options, 61];

                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on2(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on2(this_kz, this_kx,this_ky,:,:) = k_sort_on(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_2(this_kz, this_kx,this_ky,:) = on_spec_k_1(this_kz, this_kx,this_ky,:,:);
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                    elseif all(all_p_thresh == [1 1 1 0 1 1]) % [62] Replace ON1, OFF2
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 0 1], all_spec_wat, all_spec_lip, all_spec);
                            corr_options = [corr_options, 62];

                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off2(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_2(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                                on2_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                    elseif all(all_p_thresh == [1 1 1 1 0 1]) % [63] Replace ON1, OFF1
                            [corr_ph, ph_corr_coeff] = ph_k_space([1 0 1 0], all_spec_wat, all_spec_lip, all_spec);
                            corr_options = [corr_options, 63];

                            if ph_corr_coeff(1) > ph_thresh
                                k_ph_corr_on1(this_kz, this_kx, this_ky) = corr_ph(1,1);
                            else
                                k_sort_on(this_kz, this_kx,this_ky,:,:) = k_sort_on2(this_kz, this_kx,this_ky,:,:);
                                on_spec_k_1(this_kz, this_kx,this_ky,:) = on_spec_k_2(this_kz, this_kx,this_ky,:,:);
                                on1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end

                            if ph_corr_coeff(2) > ph_thresh
                                k_ph_corr_off1(this_kz, this_kx, this_ky) = corr_ph(1,2);
                            else
                                k_sort_off(this_kz, this_kx,this_ky,:,:) = k_sort_off2(this_kz, this_kx,this_ky,:,:);
                                off_spec_k_1(this_kz, this_kx,this_ky,:) = off_spec_k_2(this_kz, this_kx,this_ky,:,:); 
                                off1_replace_track(this_kz, this_kx,this_ky) = 1;
                            end
                    end

                end
            end
        end
    end
end