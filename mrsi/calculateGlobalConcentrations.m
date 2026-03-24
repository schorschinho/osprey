function [MRSCont] = calculateGlobalConcentrations(MRSCont)

    kk = 1;

    MRSCont.GlobalConc.metabolites = MRSCont.opts.MRSI.GlobalConc.metabolites;

    for mm = 1: length(MRSCont.opts.MRSI.GlobalConc.metabolites)
        metabolite = MRSCont.opts.MRSI.GlobalConc.metabolites{mm};
        quantification = MRSCont.opts.MRSI.GlobalConc.quantities{1};
    
       % Store original dimensions and create full mask
        original_dims = size(squeeze(MRSCont.quantify.(quantification).(metabolite)(:,:,MRSCont.opts.MRSI.GlobalConc.SliceIndices)));
        n_slices = length(MRSCont.opts.MRSI.GlobalConc.SliceIndices);
        
        % Reshape to vectors
        Q_full = reshape(squeeze(MRSCont.quantify.(quantification).(metabolite)(:,:,MRSCont.opts.MRSI.GlobalConc.SliceIndices)),1,[]); 
        WM_full = reshape(squeeze(MRSCont.seg.tissue.fWM(kk,:,:,MRSCont.opts.MRSI.GlobalConc.SliceIndices)),1,[]);
        GM_full = reshape(squeeze(MRSCont.seg.tissue.fGM(kk,:,:,MRSCont.opts.MRSI.GlobalConc.SliceIndices)),1,[]); 
        
        % Create combined mask for all filters
        valid_mask = (WM_full + GM_full) > MRSCont.opts.MRSI.GlobalConc.fGMpfWM;  % Tissue threshold
        valid_mask = valid_mask & (Q_full ~= 0);  % Non-zero values
        Q_mean = nanmean(Q_full(valid_mask));
        Q_std = nanstd(Q_full(valid_mask));
        valid_mask = valid_mask & (Q_full <= Q_mean + 3*Q_std) & (Q_full >= Q_mean - 3*Q_std);  % Outliers
        
        % Extract valid voxels
        Q = Q_full(valid_mask)';
        WM = WM_full(valid_mask)';
        GM = GM_full(valid_mask)';
        
        W = cat(2,GM,WM);
        MRSCont.GlobalConc.(quantification).(metabolite).WM = WM;
        MRSCont.GlobalConc.(quantification).(metabolite).GM = GM;
        MRSCont.GlobalConc.(quantification).(metabolite).Q = Q;
        
        % Solve using linear regression
        if size(W, 1) >= 2 && rank(W) == 2
            MRSCont.GlobalConc.(quantification).(metabolite).C = W \ Q;
            C_GM = MRSCont.GlobalConc.(quantification).(metabolite).C(1);
            C_WM = MRSCont.GlobalConc.(quantification).(metabolite).C(2);
            
            % Calculate predictions for valid voxels
            Q_predicted = W * MRSCont.GlobalConc.(quantification).(metabolite).C;
            residuals = Q - Q_predicted;
            MRSCont.GlobalConc.(quantification).(metabolite).residuals = residuals;
            
            % Create full-size arrays
            Q_predicted_full = nan(size(Q_full));
            residuals_full = nan(size(Q_full));
            
            % Fill in values where valid_mask is true
            Q_predicted_full(valid_mask) = Q_predicted;
            residuals_full(valid_mask) = residuals;
            
            % Reshape to original dimensions
            if n_slices == 1
                Q_predicted_image = reshape(Q_predicted_full, original_dims);
                residuals_image = reshape(residuals_full, original_dims);
                valid_mask_image = reshape(valid_mask, original_dims);
            else
                Q_predicted_image = reshape(Q_predicted_full, [original_dims(1), original_dims(2), n_slices]);
                residuals_image = reshape(residuals_full, [original_dims(1), original_dims(2), n_slices]);
                valid_mask_image = reshape(valid_mask, [original_dims(1), original_dims(2), n_slices]);
            end
            MRSCont.GlobalConc.(quantification).(metabolite).Q_predicted_image = Q_predicted_image;
            MRSCont.GlobalConc.(quantification).(metabolite).residuals_image = residuals_image;
            MRSCont.GlobalConc.(quantification).(metabolite).valid_mask_image = valid_mask_image;
            
            % Calculate R-squared
            SS_res = sum(residuals.^2);
            SS_tot = sum((Q - mean(Q)).^2);
            R2 = 1 - SS_res/SS_tot;
    
            sigma = sqrt(SS_res / (length(Q) - 2));
            WtW = W' * W;
            cov_matrix = sigma^2 * inv(WtW);
    
            % Standard errors are square roots of diagonal elements
            SE_C_GM = sqrt(cov_matrix(1,1));
            SE_C_WM = sqrt(cov_matrix(2,2));
            MRSCont.GlobalConc.(quantification).(metabolite).SE_C_GM = SE_C_GM;
            MRSCont.GlobalConc.(quantification).(metabolite).SE_C_WM = SE_C_WM;
            
            % Calculate confidence intervals (95%)
            alpha = 0.05;
            t_critical = tinv(1 - alpha/2, length(Q) - 2);
            
            CI_C_GM = [C_GM - t_critical * SE_C_GM, C_GM + t_critical * SE_C_GM];
            CI_C_WM = [C_WM - t_critical * SE_C_WM, C_WM + t_critical * SE_C_WM];
            MRSCont.GlobalConc.(quantification).(metabolite).CI_C_GM = CI_C_GM;
            MRSCont.GlobalConc.(quantification).(metabolite).CI_C_WM = CI_C_WM;

        end
    
        residual_min = min(residuals,[],'all');
        Q_max = max(Q,[],'all');
        MRSCont.GlobalConc.(quantification).(metabolite).residual_min = residual_min;
        MRSCont.GlobalConc.(quantification).(metabolite).Q_max = Q_max;
    end

    

end