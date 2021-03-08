function [MRSCont] = load_mrsi_data(MRSCont)
%[off_spec_no_lb,on_spec_no_lb] = load_mrsi_data(FileName, PathName,lb, spec_zfill, k_zfill, seq_type, water_filter, k_fft2_wat_ref, k_fft2_wat_ref_no_k_zfill, k_sort_b4_2dfft, wat_file)

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
%%
lb = 5;
spec_zfill =2;
k_zfill = 1;
seq_type = 'MEGA-PRESS';
water_filter = 0;
k_fft2_wat_ref = [];
k_fft2_wat_ref_no_k_zfill = [];
k_sort_b4_2dfft = [];
wat_file = [];
%%
% Close any remaining open figures
close all;
warning('off','all');
fileID = fopen(fullfile(MRSCont.outputFolder, 'LogFile.txt'),'a+');
if MRSCont.flags.hasMM %re_mm adding functionality to load MM data
    if ((length(MRSCont.files_mm) == 1) && (MRSCont.nDatasets>1))   %re_mm seems like specificy one MM file for a batch is also an option to plan to accomodate
        for kk=2:MRSCont.nDatasets %re_mm 
            MRSCont.files_mm{kk} = MRSCont.files_mm{1}; % re_mm allowable to specify one MM file for the whole batch
        end %re_mm 
    end   %re_mm 
    if ((length(MRSCont.files_mm) ~= MRSCont.nDatasets) )   %re_mm 
        msg = 'Number of specified MM files does not match number of specified metabolite files.'; %re_mm 
        fprintf(fileID,msg);
        error(msg);
    end   %re_mm 
end   %re_mm 
if MRSCont.flags.hasRef
    if length(MRSCont.files_ref) ~= MRSCont.nDatasets
        msg = 'Number of specified reference files does not match number of specified metabolite files.'; %re_mm 
        fprintf(fileID,msg);
        error(msg);
    end
end
if MRSCont.flags.hasWater
    if length(MRSCont.files_w) ~= MRSCont.nDatasets
        msg = 'Number of specified water files does not match number of specified metabolite files.'; %re_mm 
        fprintf(fileID,msg);
        error(msg);
    end
end

%% Get the data (loop over all datasets)
refLoadTime = tic;
reverseStr = '';
if MRSCont.flags.isGUI
    progressText = MRSCont.flags.inProgress;
end
fileID = fopen(fullfile(MRSCont.outputFolder, 'LogFile.txt'),'a+');
for kk = 1:MRSCont.nDatasets
msg = sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets);
    fprintf([reverseStr, msg]);
    fprintf(fileID,[reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    if MRSCont.flags.isGUI        
        set(progressText,'String' ,sprintf('Loading raw data from dataset %d out of %d total datasets...\n', kk, MRSCont.nDatasets));
    end
    
    if ((MRSCont.flags.didLoadData == 1 && MRSCont.flags.speedUp && isfield(MRSCont, 'raw') && (kk > length(MRSCont.raw))) || ~isfield(MRSCont.ver, 'Load') || ~strcmp(MRSCont.ver.Load,MRSCont.ver.CheckLoad))
        
        % Read in the raw metabolite data. Since the Philips DATA loader needs
        % to know the number of sub-spectra (e.g. from spectral editing), the
        % type of sequence needs to be differentiated here already.
        if MRSCont.flags.hasStatfile
            statFile = MRSCont.file_stat;
        else
            statFile = [];
        end
        disp('Opening data.')
        filename = MRSCont.files{kk};



        endian = 'l';
        type = 'float';

        disp('Reading scan parameters.')
        fname_scan_params = [filename(1:(end-4)),'list'];
        scan_params = textread(fname_scan_params, '%s');
        % Find the data lines and ignore the noise channels.
        data_lines = find(strcmp(scan_params,'STD'));
        tot_offset_idx = scan_params(data_lines(end) + 20); % last STD + 20 more offsets
        data_lines = data_lines(3:end); % ignore first 4 STDs.
        offset = str2num(scan_params{data_lines(1) + 20});

        tot_offsets = (str2num(tot_offset_idx{1}) - offset)/8192 + 1; % All offsets divided by 8192 bytes.


        disp('Reading data.')
        fp=fopen(filename, 'rb', endian);

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
%                 [k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
%                 on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,...
%                 k_ph_corr_on1, k_ph_corr_on2, k_ph_corr_off1, k_ph_corr_off2,...
%                 on1_replace_track, off1_replace_track, on2_replace_track, ...
%                 off2_replace_track, zero_replace_track, k_space_locs, corr_options]  = motion_correct_mrsi_full_ph(k_sort_on, k_sort_on2, k_sort_off, k_sort_off2,....
%                                                                                                                 on_spec_k_1, on_spec_k_2, off_spec_k_1, off_spec_k_2,...
%                                                                                                                 spec_zfill, seq_type, kx_tot, ky_tot);

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
    
   % Now we need to bring everything into the Osprey struct format
   % Extract information from SPAR files that is not in DATA/LIST
    hasSPAR = 0;
    [~,~,ext] = fileparts(filename);
    % Get the appropriate raw_act.spar
    spar_file = strrep(filename,ext, '_raw_act.spar'); % it's automatically case-insensitive

    if ~isfile(spar_file)
        spar_file = strrep(filename,ext, '_raw_act.SPAR'); % it's automatically case-insensitive
    end
    if ~isfile(spar_file)    
        if ~isempty(statFile) % Has stat csv file
            statCSV = readtable(statFile, 'Delimiter', ',','ReadVariableNames',1); % Load it
            name = statCSV.Properties.VariableNames;
            tr_idx = find(strcmp(name,'tr'));
            if isempty(tr_idx)
                tr_idx = find(strcmp(name,'TR'));
            end
            te_idx = find(strcmp(name,'te'));        
            if isempty(te_idx)
                te_idx = find(strcmp(name,'TE'));
            end
            sw_idx = find(strcmp(name,'sw'));        
            if isempty(sw_idx)
                sw_idx = find(strcmp(name,'SW'));
            end 
            Larmor_idx = find(strcmp(name,'Larmor')); 
            date_idx = find(strcmp(name,'date'));    
            seq_idx = find(strcmp(name,'seq')); 
            if ~isempty(tr_idx) && ~isempty(te_idx) && ~isempty(sw_idx) && ~isempty(Larmor_idx)
                tr = statCSV{kk,tr_idx};
                te = statCSV{kk,te_idx}; 
                sw = statCSV{kk,sw_idx}; 
                Larmor = statCSV{kk,Larmor_idx}; 
                if ~isempty(date_idx)
                    date = statCSV{kk,date_idx}{1};
                else
                    date = '';
                end
                if ~isempty(seq_idx)
                    seq = statCSV{kk,seq_idx}{1};
                else
                    seq = 'PRESS'; % Let's assume it is a PRESS sequence 
                end
            else
                msg = 'You need a SPAR file that has the same name as the .data/.list file with the extension _raw_act.spar. Otherwise you can supply TR, TE, SW (spectralwidth), Larmor (larmor frequency in MHz), and seq (localization) in the stat.csv file.';
                error(msg);    
            end
        else
            msg = 'You need a SPAR file that has the same name as the .data/.list file with the extension _raw_act.spar. Otherwise you can supply TR, TE, SW (spectralwidth), and Larmor (larmor frequency in MHz), and seq (localization) in the stat.csv file.';
            error(msg);          
        end
    else
        hasSPAR = 1;    
    end

    if hasSPAR
        % Open spar file and read parameters
        sparname = fopen(spar_file,'r');
        sparheader = textscan(sparname, '%s');
        sparheader = sparheader{1};
        sparidx=find(ismember(sparheader, 'repetition_time')==1); % TR
        tr = str2double(sparheader{sparidx+2});
        sparidx=find(ismember(sparheader, 'echo_time')==1); % TE
        te = str2double(sparheader{sparidx+2}); % 
        sparidx=find(ismember(sparheader, 'synthesizer_frequency')==1); % F0
        Larmor = str2double(sparheader{sparidx+2})/1e6;
        sparidx=find(ismember(sparheader, 'sample_frequency')==1); % readout bandwidth
        sw = str2double(sparheader{sparidx+2});

        % Read voxel geometry information.
        % THIS IS IN THE ORDER LR-AP-FH
        sparidx=find(ismember(sparheader, 'scan_date')==1); % voxel size 
        date = sparheader{sparidx+2};
        sparidx=find(ismember(sparheader, 'scan_id')==1); % voxel size 
        seq = sparheader{sparidx+2};
        sparidx=find(ismember(sparheader, 'lr_size')==1); % voxel size 
        geometry.size.lr = str2double(sparheader{sparidx+2});
        sparidx=find(ismember(sparheader, 'ap_size')==1);
        geometry.size.ap = str2double(sparheader{sparidx+2});
        sparidx=find(ismember(sparheader, 'cc_size')==1);
        geometry.size.cc = str2double(sparheader{sparidx+2});
        sparidx=find(ismember(sparheader, 'lr_off_center')==1); % voxel center offset
        geometry.pos.lr = str2double(sparheader{sparidx+2});
        sparidx=find(ismember(sparheader, 'ap_off_center')==1);
        geometry.pos.ap = str2double(sparheader{sparidx+2});
        sparidx=find(ismember(sparheader, 'cc_off_center')==1);
        geometry.pos.cc = str2double(sparheader{sparidx+2});
        sparidx=find(ismember(sparheader, 'lr_angulation')==1); % voxel angulation (radians)
        geometry.rot.lr = str2double(sparheader{sparidx+2});
        sparidx=find(ismember(sparheader, 'ap_angulation')==1);
        geometry.rot.ap = str2double(sparheader{sparidx+2});
        sparidx=find(ismember(sparheader, 'cc_angulation')==1);
        geometry.rot.cc = str2double(sparheader{sparidx+2});
        fclose(sparname);
    end

    %read in the data using the philipsDataLoad.m (adapted from PhilipsRead_data.m)
    [data] = loadRawKspace(filename);
    %[FullData,WaterData]=philipsDataLoad(filename);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine dimensions of the acquisition
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Determine number of channels
        n_coils = length(unique(data.chan));

        % Determine number of mixes
        n_mixes = length(unique(data.mix));

        % Determine number of averages per mix
        n_averages = data.kspace_properties.number_of_signal_averages;
        n_av_edit = n_averages(1);
        if n_mixes == 2
            n_av_water = n_averages(2); % if second mix exists, it is water
        end
        % This will be the NSA as specified on the exam card for the water-suppressed mix (mix = 0)
        % and the NSA as specified on the exam card for the water-suppressed mix (mix = 1);

        % Determine number of dynamics per NSA. It is not stored in the dynamics
        % field, but rather in extra attributes. Could be different for different
        % software versions, need to check!
        n_dyns = data.kspace_properties.number_of_extra_attribute_1_values;
        n_dyns_edit = n_dyns(1);
        %if n_mixes == 2
        %    n_dyns_water = n_dyns(2); % if second mix exists, it is water
        %end
        % Since dynamics are set globally on the exam card, this will be the same
        % for both - it will only be the true value for the water-suppressed mix

        % Determine number of data points per scan
        n_points = data.kspace_properties.F_resolution(1);


    if MRSCont.flags.isMEGA
        subspecs = 2;
    end

    if n_mixes == 2
        fids_w=squeeze(fids_water_ph);
    end



    %Find the magnetic field strength:
    Bo=Larmor/42.577;


    %Now create a record of the dimensions of the data array.  
    dims.t=1;
    dims.coils=0;
    dims.averages=0;
    if subspecs>1
        dims.subSpecs=2;
        sz =size(on_spec_no_lb);
        [~,t_dim] = max(sz);
        kx_dim = find(sz==kx_tot);
        ky_dim = find(sz==ky_tot);
        if subspecs == 2
            specs = zeros(sz(t_dim),2,sz(kx_dim),sz(ky_dim));
            specs(:,1,:,:) = permute(on_spec_no_lb, [t_dim kx_dim ky_dim]);
            specs(:,2,:,:) = permute(off_spec_no_lb, [t_dim kx_dim ky_dim]);
        end
        dims.Xvoxels=3;
        dims.Yvoxels=4;
        %Adding MultiVoxelInfo
        out.nXvoxels = kx_tot;
        out.nYvoxels = ky_tot;
        out.nZvoxels = 1;

    else
        dims.subSpecs=0;
    end
    dims.extras=0;

    if n_mixes == 2
        dims_w.t=1;
        dims_w.coils=0;
        dims_w.averages=2;
        dims_w.subSpecs=0;
        dims_w.extras=0;
    end

    
    if mod(size(specs,dims.t),2)==0
    %disp('Length of vector is even.  Doing normal conversion');
    fids=ifft(fftshift(specs,dims.t),[],dims.t);
    else
        %disp('Length of vector is odd.  Doing circshift by 1');
        fids=ifft(circshift(fftshift(specs,dims.t),1),[],dims.t);
    end
    if n_mixes == 2
        if mod(size(specs,dims.t),2)==0
        %disp('Length of vector is even.  Doing normal conversion');
        fids_w=ifft(fftshift(specs_w,dims.t),[],dims.t);
        else
            %disp('Length of vector is odd.  Doing circshift by 1');
            fids_w=ifft(circshift(fftshift(specs_w,dims.t),1),[],dims.t);
        end
    end

    sz=size(fids);
    if n_mixes == 2
        sz_w=size(fids_w);
    end

    %Now get relevant scan parameters:*****************************

    %Get Spectral width and Dwell Time
    spectralwidth=sw;
    dwelltime=1/spectralwidth;

    %Get TxFrq
    txfrq=Larmor*1e6;


    %Find the number of averages.  'averages' will specify the current number
    %of averages in the dataset as it is processed, which may be subject to
    %change.  'rawAverages' will specify the original number of acquired 
    %averages in the dataset, which is unchangeable.
    %FOR WATER SUPPRESSED DATA:
    if dims.subSpecs ~=0
        if dims.averages~=0
            averages=sz(dims.averages);
            rawAverages=averages*sz(dims.subSpecs);
        else
            averages= 1;
            rawAverages=sz(dims.subSpecs);
        end
    else
        if dims.averages~=0
            averages=sz(dims.averages);
            rawAverages=averages;
        else
            averages=1;
            rawAverages=1;
        end
    end

    %FOR WATER UNSUPPRESSED DATA:
    if n_mixes == 2
        if dims_w.subSpecs ~=0
            if dims_w.averages~=0
                averages_w=sz_w(dims_w.averages)*sz(dims_w.subSpecs);
                rawAverages_w=averages_w;
            else
                averages_w=sz_w(dims_w.subSpecs);
                rawAverages_w=1;
            end
        else
            if dims_w.averages~=0
                averages_w=sz_w(dims_w.averages);
                rawAverages_w=averages_w;
            else
                averages_w=1;
                rawAverages_w=1;
            end
        end
    end


    %Find the number of subspecs.  'subspecs' will specify the current number
    %of subspectra in the dataset as it is processed, which may be subject to
    %change.  'rawSubspecs' will specify the original number of acquired 
    %subspectra in the dataset, which is unchangeable.
    %FOR WATER SUPPRESSED DATA:
    if dims.subSpecs ~=0
        subspecs=sz(dims.subSpecs);
        rawSubspecs=subspecs;
    else
        subspecs=1;
        rawSubspecs=subspecs;
    end

    %FOR WATER UNSUPPRESSED DATA:
    if n_mixes == 2
        if dims_w.subSpecs ~=0
            subspecs_w=sz_w(dims.subSpecs);
            rawSubspecs_w=subspecs_w;
        else
            subspecs_w=1;
            rawSubspecs_w=subspecs_w;
        end
    end

    %****************************************************************

    %Calculate t and ppm arrays using the calculated parameters:
    f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
    ppm=f/(Bo*42.577);

    % Philips data assumes the center frequency to be 4.68 ppm:
    centerFreq = 4.68;
    ppm=ppm + centerFreq;

    t=[0:dwelltime:(sz(1)-1)*dwelltime];


    %FOR WATER SUPPRESSED DATA
    %FILLING IN DATA STRUCTURE
    out.fids=fids;
    out.specs=specs;
    out.sz=sz;
    out.ppm=ppm;  
    out.t=t;    
    out.spectralwidth=spectralwidth;
    out.dwelltime=dwelltime;
    out.txfrq=txfrq;
    out.date=date;
    out.dims=dims;
    out.Bo=Bo;
    out.averages=averages;
    out.rawAverages=rawAverages;
    out.subspecs=subspecs;
    out.rawSubspecs=rawSubspecs;
    out.seq=seq;
    out.te=te;
    out.tr=tr;
    out.pointsToLeftshift=0;
    out.centerFreq = centerFreq;
    if hasSPAR
        out.geometry = geometry;
    end

    %FILLING IN THE FLAGS
    out.flags.writtentostruct=1;
    out.flags.gotparams=1;
    out.flags.leftshifted=0;
    out.flags.filtered=0;
    out.flags.zeropadded=0;
    out.flags.freqcorrected=0;
    out.flags.phasecorrected=0;
    if averages == 1
        out.flags.averaged=1;
    else
        out.flags.averaged=0;
    end
    out.flags.addedrcvrs=1;
    out.flags.subtracted=0;
    out.flags.writtentotext=0;
    out.flags.downsampled=0;
    if out.dims.subSpecs==0
        out.flags.isISIS=0;
    else
        out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
    end

    MRSCont.raw{kk} = out;
    %FOR WATER UNSUPPRESSED DATA
    %FILLING IN DATA STRUCTURE
    if n_mixes == 2
        out_w.fids=fids_w;
        out_w.specs=specs_w;
        out_w.sz=sz_w;
        out_w.ppm=ppm;  
        out_w.t=t;    
        out_w.spectralwidth=spectralwidth;
        out_w.dwelltime=dwelltime;
        out_w.txfrq=txfrq;
        out_w.date=date;
        out_w.dims=dims_w;
        out_w.Bo=Bo;
        out_w.averages=averages_w;
        out_w.rawAverages=rawAverages_w;
        out_w.subspecs=subspecs_w;
        out_w.rawSubspecs=rawSubspecs_w;
        out_w.seq='';
        out_w.te=te;
        out_w.tr=tr;
        out_w.pointsToLeftshift=0;
        out_w.centerFreq = centerFreq;
        if hasSPAR
            out_w.geometry = geometry;
        end


        %FILLING IN THE FLAGS
        out_w.flags.writtentostruct=1;
        out_w.flags.gotparams=1;
        out_w.flags.leftshifted=0;
        out_w.flags.filtered=0;
        out_w.flags.zeropadded=0;
        out_w.flags.freqcorrected=0;
        out_w.flags.phasecorrected=0;
        out_w.flags.averaged=0;
        out_w.flags.addedrcvrs=0;
        out_w.flags.subtracted=0;
        out_w.flags.writtentotext=0;
        out_w.flags.downsampled=0;
        if out_w.dims.subSpecs==0
            out_w.flags.isISIS=0;
        else
            out_w.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
        end
        MRSCont.raw_ref{kk} = out_w;
    end        
end

fprintf('... done.\n');
time = toc(refLoadTime);
if MRSCont.flags.isGUI        
    set(progressText,'String' ,sprintf('... done.\n Elapsed time %f seconds',time));
    pause(1);
end
fprintf(fileID,'... done.\n Elapsed time %f seconds\n',time);
fclose(fileID);
% Set flag
MRSCont.flags.coilsCombined     = 1;
MRSCont.runtime.Load = time;

end