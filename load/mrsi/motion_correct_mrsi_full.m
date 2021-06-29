function [k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
    on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,track_all_p_thresh] = motion_correct_mrsi_full(k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
                                                                                on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,...
                                                                                spec_zfill, seq_type, kx_tot, ky_tot,thresholds)
% -------------- %
% --Thresholds-- %
% -------------- %
thresh = thresholds.thresh; %0.8
alt_thresh = thresholds.ph_thresh; %0.7
last_resort_thresh = thresholds.last_resort_thresh; %0.6

track_all_p_thresh = [];
dims = size(k_sort_on);
n_points = dims(end);

if strcmp(seq_type, 'MEGA multislice')

    % Eliminate bad averages.
    count_k = 0;
    disp('Eliminating bad k-space points')
    for s_idx = 1:1:size(on_spec_k_1,1)
        for this_kx = 1:size(on_spec_k_1,2)
            for this_ky = 1:size(on_spec_k_1,3)
                count_k = count_k + 1;
                %disp(sprintf('Eliminating bad k-space points of %d of %d.', count_k, 3*size(on_spec_k_1,2)*size(on_spec_k_1,3)))
                ppm_axis = linspace(-1000,1000,n_points*spec_zfill)/128 + 4.68;

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
                this_on1_wat = squeeze(on_spec_k_1(s_idx, this_kx, this_ky, met_start_wat:met_end_wat));
                this_off1_wat = squeeze(off_spec_k_1(s_idx, this_kx, this_ky, met_start_wat:met_end_wat));
                this_on2_wat = squeeze(on_spec_k_2(s_idx, this_kx, this_ky, met_start_wat:met_end_wat));
                this_off2_wat = squeeze(off_spec_k_2(s_idx, this_kx, this_ky, met_start_wat:met_end_wat));

                this_on1_lip = squeeze(on_spec_k_1(s_idx, this_kx, this_ky, met_start_lip:met_end_lip));
                this_off1_lip = squeeze(off_spec_k_1(s_idx, this_kx, this_ky, met_start_lip:met_end_lip));
                this_on2_lip = squeeze(on_spec_k_2(s_idx, this_kx, this_ky, met_start_lip:met_end_lip));
                this_off2_lip = squeeze(off_spec_k_2(s_idx, this_kx, this_ky, met_start_lip:met_end_lip));
                
                this_on1_crcho = squeeze(on_spec_k_1(s_idx, this_kx, this_ky, met_start_crcho:met_end_crcho));
                this_off1_crcho = squeeze(off_spec_k_1(s_idx, this_kx, this_ky, met_start_crcho:met_end_crcho));
                this_on2_crcho = squeeze(on_spec_k_2(s_idx, this_kx, this_ky, met_start_crcho:met_end_crcho));
                this_off2_crcho = squeeze(off_spec_k_2(s_idx, this_kx, this_ky, met_start_crcho:met_end_crcho));

                if all(this_on1_wat) && all(sum([this_kx this_ky] ~= [round(kx_tot/2) round(ky_tot/2)]))

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

                    p_on1_on2_lip = corrcoef(this_on1_lip,this_on2_lip);
                    p_on1_on2_lip = p_on1_on2_lip(2);
                    p_off1_off2_lip = corrcoef(this_off1_lip,this_off2_lip);
                    p_off1_off2_lip = p_off1_off2_lip(2);
                    p_on1_off1_lip = corrcoef(this_on1_lip,this_off1_lip);
                    p_on1_off1_lip = p_on1_off1_lip(2);
                    p_on2_off2_lip = corrcoef(this_on2_lip,this_off2_lip);
                    p_on2_off2_lip = p_on2_off2_lip(2);
                    
                    p_on1_on2_crcho = corrcoef(this_on1_crcho,this_on2_crcho);
                    p_on1_on2_crcho = p_on1_on2_crcho(2);
                    p_off1_off2_crcho = corrcoef(this_off1_crcho,this_off2_crcho);
                    p_off1_off2_crcho = p_off1_off2_crcho(2);
                    p_on1_off1_crcho = corrcoef(this_on1_crcho,this_off1_crcho);
                    p_on1_off1_crcho = p_on1_off1_crcho(2);
                    p_on2_off2_crcho = corrcoef(this_on2_crcho,this_off2_crcho);
                    p_on2_off2_crcho = p_on2_off2_crcho(2);

                    % alternatives, Lip, ON1 & OFF2, ON2 & OFF1
                    p_on1_off2_lip = corrcoef(this_on1_lip,this_off2_lip);
                    p_on1_off2_lip = p_on1_off2_lip(2);
                    p_on2_off1_lip = corrcoef(this_on2_lip,this_off1_lip);
                    p_on2_off1_lip = p_on2_off1_lip(2);
                    
                    % alternatives, wat, ON1 & OFF2, ON2 & OFF1
                    p_on1_off2_wat = corrcoef(this_on1_wat,this_off2_wat);
                    p_on1_off2_wat = p_on1_off2_wat(2);
                    p_on2_off1_wat = corrcoef(this_on2_wat,this_off1_wat);
                    p_on2_off1_wat = p_on2_off1_wat(2);
                    
                    % alternatives, Cr Cho, ON1 & OFF2, ON2 & OFF1
                    p_on1_off2_crcho = corrcoef(this_on1_crcho,this_off2_crcho);
                    p_on1_off2_crcho = p_on1_off2_crcho(2);
                    p_on2_off1_crcho = corrcoef(this_on2_crcho,this_off1_crcho);
                    p_on2_off1_crcho = p_on2_off1_crcho(2);
                    
                    % For water, only compare between ONs or just between
                    % OFFs because water is saturated in GSH.
                    p_wat = [p_on1_on2_wat p_off1_off2_wat p_on1_off1_wat p_on2_off2_wat];
                    all_p_thresh_wat = p_wat < thresh;
                    p_lip = [p_on1_on2_lip p_off1_off2_lip p_on1_off1_lip p_on2_off2_lip];
                    all_p_thresh_lip = p_lip < thresh;
                    p_crcho = [p_on1_on2_crcho p_off1_off2_crcho p_on1_off1_crcho p_on2_off2_crcho];
                    all_p_thresh_crcho = p_crcho < thresh;

                    all_p_thresh = all_p_thresh_wat | all_p_thresh_lip | all_p_thresh_crcho;
                    
                    all_p_thresh_alt_lip = [p_on1_off2_lip p_on2_off1_lip] > thresh;
                    all_p_thresh_alt_wat = [p_on1_off2_wat p_on2_off1_wat] > thresh;
                    all_p_thresh_alt_crcho = [p_on1_off2_crcho p_on2_off1_crcho] > thresh;

                    all_p_thresh_alt_lip2 = [p_on1_off2_lip p_on2_off1_lip] > alt_thresh;
                    all_p_thresh_alt_wat2 = [p_on1_off2_wat p_on2_off1_wat] > alt_thresh;
                    all_p_thresh_alt_crcho2 = [p_on1_off2_crcho p_on2_off1_crcho] > alt_thresh;

                    all_p_thresh_lip_zero = [p_on1_on2_lip p_off1_off2_lip p_on1_off1_lip p_on2_off2_lip] < last_resort_thresh;
                    all_p_thresh_wat_zero = [p_on1_on2_wat p_off1_off2_wat p_on1_off1_wat p_on2_off2_wat] < last_resort_thresh;
                    all_p_thresh_crcho_zero = [p_on1_on2_crcho p_off1_off2_crcho p_on1_off1_crcho p_on2_off2_crcho] < last_resort_thresh;

                    all_p_thresh_zero = all_p_thresh_wat_zero | all_p_thresh_lip_zero | all_p_thresh_crcho_zero;

                    all_p_thresh_alt = all_p_thresh_alt_lip & all_p_thresh_alt_wat & all_p_thresh_alt_crcho;
                    all_p_thresh_alt2 = all_p_thresh_alt_lip2 & all_p_thresh_alt_wat2 & all_p_thresh_alt_crcho2;
                    
                    if (s_idx == 3)
                        dkdfd = 3;
                    end


                    % Set to zero if all are different, ON1 and ON2 are
                    % different, OFF1 and OFF2 are different.
                    if all(all_p_thresh_zero) % All sub-acquisitions are different from each other
                        % If all the sub-acquisitions are different from
                        % one another, just replace with zeros (for now,
                        % interpolation in the future?)
                        k_sort_on(s_idx,this_kx,this_ky,:,:) = zeros([size(k_sort_on,4), size(k_sort_on,5)]);
                        k_sort_on2(s_idx,this_kx,this_ky,:,:) = zeros([size(k_sort_on,4), size(k_sort_on,5)]);
                        k_sort_off(s_idx,this_kx,this_ky,:,:) = zeros([size(k_sort_on,4), size(k_sort_on,5)]);
                        k_sort_off2(s_idx,this_kx,this_ky,:,:) = zeros([size(k_sort_on,4), size(k_sort_on,5)]);

                        on_spec_k_1(s_idx,this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,5));
                        on_spec_k_2(s_idx,this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,5));
                        off_spec_k_1(s_idx,this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,5));
                        off_spec_k_2(s_idx,this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,5));
                    else
                        if all(all_p_thresh == [0 1 0 1]) % OFF2 is bad
                            k_sort_off2(s_idx,this_kx,this_ky,:,:) = k_sort_off(s_idx,this_kx,this_ky,:,:);
                            off_spec_k_2(s_idx,this_kx,this_ky,:) = off_spec_k_1(s_idx,this_kx,this_ky,:,:);
                            % KLC 04/12/18
                         % ---------------------------------------------- %     
                         % ----- Difficult to determine thresholds ------ %
                         % ---------------------------------------------- %
                         elseif all(all_p_thresh == [0 1 1 1])   % OFF1 and OFF2 are bad, keep the least bad one, still CC > 0.8    %KLC 04/12/18
                            this_p_thresh = [p_lip(3) < alt_thresh | p_wat(3) < alt_thresh | p_crcho(3) < alt_thresh,...
                                             p_lip(4) < alt_thresh | p_wat(4) < alt_thresh | p_crcho(4) < alt_thresh];

                            if all(this_p_thresh == [1 0])  % Replace OFF1
                                k_sort_off(s_idx,this_kx,this_ky,:,:) = k_sort_off2(s_idx,this_kx,this_ky,:,:);
                                off_spec_k_1(s_idx,this_kx,this_ky,:) = off_spec_k_2(s_idx,this_kx,this_ky,:);
                            elseif all(this_p_thresh == [0 1]) % Replace OFF2
                                k_sort_off2(s_idx,this_kx,this_ky,:,:) = k_sort_off(s_idx,this_kx,this_ky,:,:);
                                off_spec_k_2(s_idx,this_kx,this_ky,:) = off_spec_k_1(s_idx,this_kx,this_ky,:,:);
                                
                            end
                         elseif all(all_p_thresh == [1 0 1 1])   % ON2 and ON1 are bad(ish?), keep the least bad one, still CC > 0.8    %KLC 04/12/18
                            this_p_thresh = [p_lip(3) < alt_thresh | p_wat(3) < alt_thresh | p_crcho(3) < alt_thresh,...
                                             p_lip(4) < alt_thresh | p_wat(4) < alt_thresh | p_crcho(4) < alt_thresh];

                            if all(this_p_thresh == [1 0])  % Replace ON1
                                k_sort_on(s_idx,this_kx,this_ky,:,:) = k_sort_on2(s_idx,this_kx,this_ky,:,:);
                                on_spec_k_1(s_idx,this_kx,this_ky,:) = on_spec_k_2(s_idx,this_kx,this_ky,:,:);
                            elseif all(this_p_thresh == [0 1]) % Replace ON2
                                k_sort_on2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                                on_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                            end

                         elseif all(all_p_thresh == [1 1 0 0]) % One of the two ON/OFF pairs is bad.                %KLC 04/12/18
                            this_p_thresh = [p_wat(1) < alt_thresh | p_lip(1) < alt_thresh | p_crcho(1) < alt_thresh,...
                                             p_wat(2) < alt_thresh | p_lip(2) < alt_thresh | p_crcho(2) < alt_thresh];
                            if all(this_p_thresh == [1 0])
                                if all(all_p_thresh_alt2 == [1 0])   % ON2 is the bad one
                                    k_sort_on2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                                    on_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                                elseif all(all_p_thresh_alt2 == [0 1])   % ON1 is the bad one
                                    k_sort_on(s_idx,this_kx,this_ky,:,:) = k_sort_on2(s_idx,this_kx,this_ky,:,:);
                                    on_spec_k_1(s_idx,this_kx,this_ky,:) = on_spec_k_2(s_idx,this_kx,this_ky,:);
                                end
                            elseif all(this_p_thresh == [0 1])
                                if all(all_p_thresh_alt2 == [1 0])   % Replace OFF1
                                    k_sort_on2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                                    on_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                                elseif all(all_p_thresh_alt2 == [0 1])   % OFF2 is the bad one
                                    k_sort_off2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                                    off_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                                end
                            end
                        elseif all(all_p_thresh == [0 0 1 1])                                                   % KLC 04/12/18
                            this_p_thresh = [p_lip(3) < alt_thresh | p_wat(3) < alt_thresh | p_crcho(3) < alt_thresh,...
                                             p_lip(4) < alt_thresh | p_wat(4) < alt_thresh | p_crcho(4) < alt_thresh];
                            if all(this_p_thresh == [0 1])
                                if all(all_p_thresh_alt2 == [1 0])   % ON2 is the bad one
                                    k_sort_on2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                                    on_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                                elseif all(all_p_thresh_alt2 == [0 1])   % OFF2 is the bad one
                                    k_sort_off2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                                    off_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                                end
                            elseif all(this_p_thresh == [1 0])
                                if all(all_p_thresh_alt2 == [1 0])   % ON1 is the bad one
                                    k_sort_on2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                                    on_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                                elseif all(all_p_thresh_alt2 == [0 1])   % OFF1 is the bad one
                                    k_sort_off2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                                    off_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                                end
                            end
                                                    % KLC 04/12/18
                         % ------------------------------------------------------ %     
                         % ----- END difficult to determine thresholds END ------ %
                         % ------------------------------------------------------ %
                        elseif all(all_p_thresh == [1 1 0 1]) % ON2 and OFF2 are bad
                            k_sort_on2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                            k_sort_off2(s_idx,this_kx,this_ky,:,:) = k_sort_off(s_idx,this_kx,this_ky,:,:);
                            on_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:,:);
                            off_spec_k_2(s_idx,this_kx,this_ky,:) = off_spec_k_1(s_idx,this_kx,this_ky,:,:);

                        elseif all(all_p_thresh == [1 1 1 0]) % ON1 and OFF1 are bad
                            k_sort_on(s_idx,this_kx,this_ky,:,:) = k_sort_on2(s_idx,this_kx,this_ky,:,:);
                            k_sort_off(s_idx,this_kx,this_ky,:,:) = k_sort_off2(s_idx,this_kx,this_ky,:,:);
                            on_spec_k_1(s_idx,this_kx,this_ky,:) = on_spec_k_2(s_idx,this_kx,this_ky,:,:);
                            off_spec_k_1(s_idx,this_kx,this_ky,:) = off_spec_k_2(s_idx,this_kx,this_ky,:,:);

                        elseif all(all_p_thresh == [1 0 1 0])   % ON1 is bad
    %                           
                            k_sort_on(s_idx,this_kx,this_ky,:,:) = k_sort_on2(s_idx,this_kx,this_ky,:,:);
                            on_spec_k_1(s_idx,this_kx,this_ky,:) = on_spec_k_2(s_idx,this_kx,this_ky,:,:);
                        elseif all(all_p_thresh == [0 0 1 0]) % Either ON1 or OFF1 is the bad one, query which one it is.

                            if all(all_p_thresh_alt == [1 0])   % OFF1 is the bad one
                                k_sort_off(s_idx,this_kx,this_ky,:,:) = k_sort_off2(s_idx,this_kx,this_ky,:,:);
                                off_spec_k_1(s_idx,this_kx,this_ky,:) = off_spec_k_2(s_idx,this_kx,this_ky,:);
                            elseif all(all_p_thresh_alt == [0 1])   % ON1 is the bad one
                                k_sort_on(s_idx,this_kx,this_ky,:,:) = k_sort_on2(s_idx,this_kx,this_ky,:,:);
                                on_spec_k_1(s_idx,this_kx,this_ky,:) = on_spec_k_2(s_idx,this_kx,this_ky,:);
                            elseif all(all_p_thresh_alt == [0 0])   % Just replace with ON2 and OFF2
                                k_sort_off(s_idx,this_kx,this_ky,:,:) = k_sort_off2(s_idx,this_kx,this_ky,:,:);
                                k_sort_on(s_idx,this_kx,this_ky,:,:) = k_sort_on2(s_idx,this_kx,this_ky,:,:);

                                off_spec_k_1(s_idx,this_kx,this_ky,:) = off_spec_k_2(s_idx,this_kx,this_ky,:);
                                on_spec_k_1(s_idx,this_kx,this_ky,:) = on_spec_k_2(s_idx,this_kx,this_ky,:);
                            end
                        elseif all(all_p_thresh == [0 0 0 1]) % ON2 or OFF2 is bad.
                            if all(all_p_thresh_alt == [1 0])   % ON2 is the bad one
                                k_sort_on2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                                on_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                            elseif all(all_p_thresh_alt == [0 1])   % OFF2 is the bad one
                                k_sort_off2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                                off_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                            elseif all(all_p_thresh_alt == [0 0])  % Just replace with ON1 and OFF1
                                k_sort_off2(s_idx,this_kx,this_ky,:,:) = k_sort_off(s_idx,this_kx,this_ky,:,:);
                                k_sort_on2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);

                                off_spec_k_2(s_idx,this_kx,this_ky,:) = off_spec_k_1(s_idx,this_kx,this_ky,:);
                                on_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                            end
                        elseif all(all_p_thresh == [1 0 0 1]) % ON2 is bad

                            k_sort_on2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                            on_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                        elseif all(all_p_thresh == [0 1 1 0])   % OFF1 is bad
                            k_sort_off(s_idx,this_kx,this_ky,:,:) = k_sort_off2(s_idx,this_kx,this_ky,:,:);
                            off_spec_k_1(s_idx,this_kx,this_ky,:) = off_spec_k_2(s_idx,this_kx,this_ky,:);
                        elseif all(all_p_thresh == [1 0 0 0])   % Either ON1 or ON2 is bad
                            if all(all_p_thresh_alt == [1 0])   % ON2 is the bad one
                                k_sort_on2(s_idx,this_kx,this_ky,:,:) = k_sort_on(s_idx,this_kx,this_ky,:,:);
                                on_spec_k_2(s_idx,this_kx,this_ky,:) = on_spec_k_1(s_idx,this_kx,this_ky,:);
                            elseif all(all_p_thresh_alt == [0 1])   % ON1 is the bad one
                                k_sort_on(s_idx,this_kx,this_ky,:,:) = k_sort_on2(s_idx,this_kx,this_ky,:,:);
                                on_spec_k_1(s_idx,this_kx,this_ky,:) = on_spec_k_2(s_idx,this_kx,this_ky,:);
                            end
                        elseif all(all_p_thresh == [0 1 0 0])   % Either OFF1 or OFF2 is bad
                            if all(all_p_thresh_alt == [1 0])   % OFF1 is bad.
                                k_sort_off(s_idx,this_kx,this_ky,:,:) = k_sort_off2(s_idx,this_kx,this_ky,:,:);
                                off_spec_k_1(s_idx,this_kx,this_ky,:) = off_spec_k_2(s_idx,this_kx,this_ky,:);
                            elseif all(all_p_thresh_alt == [0 1])   % OFF2 is bad
                                k_sort_off2(s_idx,this_kx,this_ky,:,:) = k_sort_off(s_idx,this_kx,this_ky,:,:);
                                off_spec_k_2(s_idx,this_kx,this_ky,:) = off_spec_k_1(s_idx,this_kx,this_ky,:);
                            end
                        end
                    end

                end
            end
        end
    end
else
    
    % Eliminate bad averages.
    count_k = 0;    
    disp('Eliminating bad k-space points')
    for this_kx = 1:size(on_spec_k_1,1)
        for this_ky = 1:size(on_spec_k_1,2)
            count_k = count_k + 1;
            %disp(sprintf('Eliminating bad k-space points of %d of %d.', count_k, 3*size(on_spec_k_1,2)*size(on_spec_k_1,3)))
            ppm_axis = linspace(-1000,1000,n_points*spec_zfill)/128 + 4.68;

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
            this_on1_wat = squeeze(on_spec_k_1( this_kx, this_ky, met_start_wat:met_end_wat));
            this_off1_wat = squeeze(off_spec_k_1( this_kx, this_ky, met_start_wat:met_end_wat));
            this_on2_wat = squeeze(on_spec_k_2( this_kx, this_ky, met_start_wat:met_end_wat));
            this_off2_wat = squeeze(off_spec_k_2( this_kx, this_ky, met_start_wat:met_end_wat));

            this_on1_lip = squeeze(on_spec_k_1( this_kx, this_ky, met_start_lip:met_end_lip));
            this_off1_lip = squeeze(off_spec_k_1( this_kx, this_ky, met_start_lip:met_end_lip));
            this_on2_lip = squeeze(on_spec_k_2( this_kx, this_ky, met_start_lip:met_end_lip));
            this_off2_lip = squeeze(off_spec_k_2( this_kx, this_ky, met_start_lip:met_end_lip));

            this_on1_crcho = squeeze(on_spec_k_1( this_kx, this_ky, met_start_crcho:met_end_crcho));
            this_off1_crcho = squeeze(off_spec_k_1( this_kx, this_ky, met_start_crcho:met_end_crcho));
            this_on2_crcho = squeeze(on_spec_k_2( this_kx, this_ky, met_start_crcho:met_end_crcho));
            this_off2_crcho = squeeze(off_spec_k_2( this_kx, this_ky, met_start_crcho:met_end_crcho));

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
%                 p_on1_off1_lip = corrcoef(this_on1_lip,this_off1_lip);
%                 p_on1_off1_lip = p_on1_off1_lip(2);
%                 p_on2_off2_lip = corrcoef(this_on2_lip,this_off2_lip);
%                 p_on2_off2_lip = p_on2_off2_lip(2);

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
                all_p_thresh_wat = p_wat < thresh;
                p_lip = [p_on1_on2_lip p_off1_off2_lip]; %p_on1_off1_lip p_on2_off2_lip
                all_p_thresh_lip = p_lip < thresh;
                p_crcho = [p_on1_on2_crcho p_on1_off1_crcho p_on1_off2_crcho p_on2_off1_crcho p_on2_off2_crcho p_off1_off2_crcho];
                all_p_thresh_crcho = p_crcho < thresh;
                
                all_p = [mean([p_lip(1) p_wat(1)]), mean([p_wat(2)]),...
                        mean([p_wat(3)]), mean([p_wat(4)]),... % p_crcho(4)
                        mean([p_wat(5)]), mean([p_lip(2) p_wat(6)])];

                all_p_thresh = [all_p_thresh_lip(1) | all_p_thresh_wat(1),...% | all_p_thresh_crcho(1),...
                                all_p_thresh_wat(2:5),...% | all_p_thresh_crcho(2:5),...
                                all_p_thresh_lip(2) | all_p_thresh_wat(6)];% | all_p_thresh_crcho(6)];% | all_p_thresh_lip;

                % Last resort, all of them don't match, unless one pair is > 0.7          
                all_p_thresh_wat_zero = p_wat < last_resort_thresh;
                all_p_thresh_lip_zero = p_lip < last_resort_thresh;
                all_p_thresh_crcho_zero = p_crcho < last_resort_thresh;

                all_p_thresh_zero = [all_p_thresh_lip_zero(1) | all_p_thresh_wat_zero(1) | all_p_thresh_crcho_zero(1),...
                                all_p_thresh_wat_zero(2:5) | all_p_thresh_crcho_zero(2:5),...
                                all_p_thresh_lip_zero(2) | all_p_thresh_wat_zero(6) | all_p_thresh_crcho_zero(6)];% | all_p_thresh_lip;
                
                track_all_p_thresh = [track_all_p_thresh; all_p_thresh];
                
                if (this_kx == 6) && (this_ky == 6)
                    xes = 3;
                end
                            
                % Set to zero if all are different, ON1 and ON2 are
                % different, OFF1 and OFF2 are different.
                if all(all_p_thresh_zero) % All sub-acquisitions are different from each other
                    % If all the sub-acquisitions are different from
                    % one another, just replace with zeros (for now,
                    % interpolation in the future?)
                    k_sort_on(this_kx,this_ky,:,:) = zeros([size(k_sort_on,3), size(k_sort_on,4)]);
                    k_sort_on2(this_kx,this_ky,:,:) = zeros([size(k_sort_on,3), size(k_sort_on,4)]);
                    k_sort_off(this_kx,this_ky,:,:) = zeros([size(k_sort_on,3), size(k_sort_on,4)]);
                    k_sort_off2(this_kx,this_ky,:,:) = zeros([size(k_sort_on,3), size(k_sort_on,4)]);

                    on_spec_k_1(this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,4));
                    on_spec_k_2(this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,4));
                    off_spec_k_1(this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,4));
                    off_spec_k_2(this_kx,this_ky,:) = zeros(1,spec_zfill*size(k_sort_on,4));
                else
                    if all(all_p_thresh == [1 0 0 0 0 0]) % [1] Replace ON1 OR ON2
                        this_comp = [sum(all_p(2:3)), sum(all_p(4:5))];
                        keep = find(this_comp == max(this_comp));
                        
                        if keep == 1 %Replace ON2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 1 0 0 0 0]) % [2] Replace ON1 w/ON2
                        k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                        on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:); 
                    elseif all(all_p_thresh == [1 1 1 0 0 0]) % [3] Replace ON1 w/ON2
                        k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                        on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:); 
                    elseif all(all_p_thresh == [1 1 1 1 0 0]) % [4] Replace ON1 w/ON2, and replace OFF1 with OFF2
                        k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                        on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                        
                        k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                        off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:); 
                    elseif all(all_p_thresh == [1 1 1 1 1 0]) % [5] Replace ON1 or ON2
                        this_comp = [sum(all_p(2:3)), sum(all_p(4:5))];
                        keep = find(this_comp == max(this_comp));
                        
                        if keep == 1 %Replace ON2 with ON1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 1 1 1 1 1]) % [6] Pick the best ON/OFF pair
                        this_comp = all_p(2:5);
                        keep = find(this_comp == max(this_comp));
                        
                        if keep == 1 %Replace ON2 with ON1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 3 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        elseif keep == 4 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 1 0 0 0 0]) % [7]. Replace ON1 or OFF1
                        this_comp = [sum(all_p([1 3])), sum(all_p([4 6]))];
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace OFF1
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 1 1 0 0 0]) % [8] Replace ON1 w/ON2
                        k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                        on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [0 1 1 1 0 0]) % [9] Replace ON1 w/ON2
                        k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                        on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [0 1 1 1 1 0]) % [10] Pick the best ON/OFF pair
                        this_comp = all_p(2:5);
                        keep = find(this_comp == max(this_comp));
                        
                        if keep == 1 %Replace ON2 with ON1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 3 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        elseif keep == 4 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 1 1 1 1 1]) % [11] Pick the best ON/OFF pair
                        this_comp = all_p(2:5);
                        keep = find(this_comp == max(this_comp));
                        
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 3 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        elseif keep == 4 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 0 1 0 0 0]) % [12] Replace ON1 w/ON2
                        k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                        on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [1 0 0 1 0 0]) % [13] Replace ON2 w/ON1
                        k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                        on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [1 0 0 0 1 0]) % [14] Replace ON2 w/ON1
                        k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                        on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [1 0 0 0 0 1]) % [15] Pick the best ON/OFF pair
                        this_comp = all_p(2:5);
                        keep = find(this_comp == max(this_comp));
                        
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 3 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        elseif keep == 4 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 1 0 1 0 0]) % [16] Find max of 3, 5
                        this_comp = all_p([3 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 1 0 0 1 0]) % [17] Choose best ON/OFF pair other than #2,5
                        this_comp = all_p(3:4);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 1 0 0 0 1]) % [18] Replace OFF1
                        k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                        off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);                      
                    elseif all(all_p_thresh == [0 0 1 1 0 0]) % [19] Pick the best ON/OFF pair except #3,4
                        this_comp = all_p([2 5]);
                        keep = find(this_comp == max(this_comp));
                        
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 0 1 0 1 0]) % [20] Replace OFF2
                        k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                        off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [0 0 1 0 0 1]) % [21] Replace OFF2
                        k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                        off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                    elseif all(all_p_thresh == [0 0 0 1 1 0]) % [22] Replace ON2 w/ON1
                        k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                        on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [0 0 0 1 0 1]) % [23] Replace OFF1
                        k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                        off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);   
                    elseif all(all_p_thresh == [0 0 0 0 1 1]) % [24] Replace OFF2
                        k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                        off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                    elseif all(all_p_thresh == [0 0 1 0 0 0]) % [25] Replace ON1 or OFF2
                        this_comp = [sum(all_p([1 2])), sum(all_p([5 6]))];
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace OFF2
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 0 0 1 0 0]) % [26] Replace ON1 or OFF2
                        this_comp = [sum(all_p([1 5])), sum(all_p([2 6]))];
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace OFF1
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 0 0 0 1 0]) % [27] Replace ON1 or OFF2
                        this_comp = [sum(all_p([1 4])), sum(all_p([3 6]))];
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace OFF2
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 0 0 0 0 1]) % [28] Replace ON1 or OFF2
                        this_comp = [sum(all_p([2 4])), sum(all_p([3 5]))];
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace OFF2
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace OFF1
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:); 
                        end
                    elseif all(all_p_thresh == [1 1 0 1 0 0]) % [29] Replace ON1 or OFF2
                        this_comp = all_p([3 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                        end
                        % Replace OFF1
                        k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                        off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:); 
                    elseif all(all_p_thresh == [1 1 0 0 1 0]) % [30] Replace ON1 or ON2, replace OFF1 or OFF2
                        this_comp = all_p([3 4]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            % Replace OFF1
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            % Replace OFF2
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        end
                    elseif all(all_p_thresh == [1 1 0 0 0 1]) % [31] Replace ON2 and OFF1, etc. keep good pair
                        this_comp = all_p(3:5);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        elseif keep == 3 %Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 0 1 1 0 0]) % [32] Replace ON2 and OFF2, replace ON1 and OFF1
                        this_comp = all_p([2 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                        
                    elseif all(all_p_thresh == [1 0 1 0 1 0]) % [33] Replace ON1 or ON2, replace OFF1 or OFF2
                        this_comp = all_p([2 4]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 0 1 0 0 1]) % [34] Keep good pair
                        this_comp = all_p([2 4 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        elseif keep == 3 %Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 0 0 1 1 0]) % [35] Replace ON2s
                        k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                        on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [1 0 0 1 0 1]) % [36] Pick best pair
                        this_comp = all_p([2 3 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 3 %Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 0 0 0 1 1]) % [37] Pick the best ON/OFF pair
                        this_comp = all_p(2:4);
                        keep = find(this_comp == max(this_comp));
                        
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 3 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 1 1 0 1 0]) % [38] Replace ON1, OFF2
                        k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                        on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);

                        k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                        off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [0 1 1 0 0 1]) % [39] Pick best ON, OFF pair
                        this_comp = all_p([4 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 1 0 1 1 0]) % [40] Replace ON2, OFF1
                        k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                        on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);

                        k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                        off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [0 1 0 1 0 1]) % [41] Pick the best pair
                        this_comp = all_p([3 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 1 0 0 1 1]) % [42] Pick the best ON/OFF pair
                        this_comp = all_p(3:4);
                        keep = find(this_comp == max(this_comp));
                        
                        if keep == 1 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 0 1 1 1 0]) % [43] Replace ON2, OFF2
                        k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                        on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);

                        k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                        off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [0 0 1 1 0 1]) % [44] Pick the best ON/OFF pair
                        this_comp = all_p([2 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 0 1 0 1 1]) % [45] Pick the best ON/OFF pair
                        this_comp = all_p([2 4]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [0 0 0 1 1 1]) % [46] Pick the best ON/OFF pair
                        this_comp = all_p([2 3]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 1 1 0 1 0]) % [47] Replace ON1, OFF2
                        k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                        on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);

                        k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                        off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [1 1 1 0 0 1]) % [48] Pick best ON, OFF pair
                        this_comp = all_p([4 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 1 0 1 1 0]) % [49] Replace ON2, OFF1
                        k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                        on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);

                        k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                        off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [1 1 0 1 0 1]) % [50] Pick the best pair
                        this_comp = all_p([3 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        elseif keep == 2 %Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 1 0 0 1 1]) % [51] Replace ON1 or ON2, replace OFF1 or OFF2
                        this_comp = all_p([3 4]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            % Replace OFF1
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            % Replace OFF2
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        end
                    elseif all(all_p_thresh == [1 0 1 1 1 0]) % [52] Replace ON2, OFF2
                        k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                        on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);

                        k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                        off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [1 0 1 1 0 1]) % [53] Pick best pair
                        this_comp = all_p([2 5]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 0 1 0 1 1]) % [54] Pick the best ON/OFF pair
                        this_comp = all_p([2 4]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                        end
                    elseif all(all_p_thresh == [1 0 0 1 1 1]) % [55] Pick the best ON/OFF pair
                        this_comp = all_p([2 3]);
                        keep = find(this_comp == max(this_comp));
                        if keep == 1 %Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:); 
                        elseif keep == 2 %Replace ON1, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                        end
                        
                    elseif all(all_p_thresh == [0 1 1 1 0 1]) % [56] Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [0 1 1 0 1 1]) % [57] Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [0 1 0 1 1 1]) % [58] Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [0 0 1 1 1 1]) % [59] Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [1 0 1 1 1 1]) % [60] Replace ON2, OFF2
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [1 1 0 1 1 1]) % [61] Replace ON2, OFF1
                            k_sort_on2(this_kx,this_ky,:,:) = k_sort_on(this_kx,this_ky,:,:);
                            on_spec_k_2(this_kx,this_ky,:) = on_spec_k_1(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [1 1 1 0 1 1]) % [62] Replace ON1, OFF2
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off2(this_kx,this_ky,:,:) = k_sort_off(this_kx,this_ky,:,:);
                            off_spec_k_2(this_kx,this_ky,:) = off_spec_k_1(this_kx,this_ky,:,:);
                    elseif all(all_p_thresh == [1 1 1 1 0 1]) % [63] Replace ON1, OFF1
                            k_sort_on(this_kx,this_ky,:,:) = k_sort_on2(this_kx,this_ky,:,:);
                            on_spec_k_1(this_kx,this_ky,:) = on_spec_k_2(this_kx,this_ky,:,:);
                            
                            k_sort_off(this_kx,this_ky,:,:) = k_sort_off2(this_kx,this_ky,:,:);
                            off_spec_k_1(this_kx,this_ky,:) = off_spec_k_2(this_kx,this_ky,:,:);
                    end

                end
            end
        end
    end
end