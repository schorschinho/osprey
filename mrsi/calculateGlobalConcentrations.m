function [MRSCont] = calculateGlobalConcentrations(MRSCont)
%% [MRSCont] = calculateGlobalConcentrations(MRSCont)
%   This function calculates global concentrations using a linear
%   regression appoach as drescribed in Tal A, Kirov II, Grossman RI, Gonen O. 
%   The role of gray and white matter segmentation in quantitative proton MR 
%   spectroscopic imaging. NMR Biomed. 2012;25(12):1392-1400. doi:10.1002/nbm.2812
%
%   USAGE:
%       MRSCont = calculateGlobalConcentrations(MRSCont);
%
%   INPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   OUTPUTS:
%       MRSCont     = Osprey MRS data container.
%
%   AUTHOR:
%       Helge Zöllner (Johns Hopkins University, 2025-10-31)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2025-10-31: First version of the code.
%%  Linear-regression analysis
%   Perform the analysis of metabolite concentration vs tissue fraction to
%   get an estimate of the global 100% WM or GM concentration
    kk = 1; % Currently working for single subject only

    MRSCont.GlobalConc.metabolites = MRSCont.opts.MRSI.GlobalConc.metabolites;  % Get the metabolites of interst

    for mm = 1: length(MRSCont.opts.MRSI.GlobalConc.metabolites)                % Loop over the metabolites of interest
        metabolite = MRSCont.opts.MRSI.GlobalConc.metabolites{mm};              % metabolite name
        quantification = MRSCont.opts.MRSI.GlobalConc.quantities{1};            % quantification name
    
       % Store original dimensions and create full mask
        original_dims = size(squeeze(MRSCont.quantify.(quantification).(metabolite)(:,:,MRSCont.opts.MRSI.GlobalConc.SliceIndices)));   % Get original dims
        n_slices = length(MRSCont.opts.MRSI.GlobalConc.SliceIndices);           % Get number of slices
        
        % Reshape to vectors 
        Q_full = reshape(squeeze(MRSCont.quantify.(quantification).(metabolite)(:,:,MRSCont.opts.MRSI.GlobalConc.SliceIndices)),1,[]);  % Pick the slices as described in the indices for a certain metabolite
        WM_full = reshape(squeeze(MRSCont.seg.tissue.fWM(kk,:,:,MRSCont.opts.MRSI.GlobalConc.SliceIndices)),1,[]); % Pick white matter content of same slices
        GM_full = reshape(squeeze(MRSCont.seg.tissue.fGM(kk,:,:,MRSCont.opts.MRSI.GlobalConc.SliceIndices)),1,[]); % Pick gray matter content of same slices
        
        % Create combined mask for all filters
        valid_mask = (WM_full + GM_full) > MRSCont.opts.MRSI.GlobalConc.fGMpfWM;  % Tissue threshold from container
        valid_mask = valid_mask & (Q_full ~= 0);                                  % Non-zero values
        Q_mean = nanmean(Q_full(valid_mask));                                     % Mean of metabolites
        Q_std = nanstd(Q_full(valid_mask));                                       % SD of metabolites
        valid_mask = valid_mask & (Q_full <= Q_mean + 3*Q_std) & (Q_full >= Q_mean - 3*Q_std);  % Find outliers
        
        % Extract valid voxels
        Q = Q_full(valid_mask)';                                                % Get metabolites after filtering
        WM = WM_full(valid_mask)';                                              % Remove filtered white matter fractions
        GM = GM_full(valid_mask)';                                              % Remove filtered gray matter fractions
        
        W = cat(2,GM,WM);                                                       % Setup matrix
        MRSCont.GlobalConc.(quantification).(metabolite).WM = WM;               % Store white matter fraction
        MRSCont.GlobalConc.(quantification).(metabolite).GM = GM;               % Store gray matter fraction
        MRSCont.GlobalConc.(quantification).(metabolite).Q = Q;                 % Store metabolites
        
        % Solve using linear regression
        if size(W, 1) >= 2 && rank(W) == 2
            MRSCont.GlobalConc.(quantification).(metabolite).C = W \ Q;         % Linear regression between tissue and metabolite
            C_GM = MRSCont.GlobalConc.(quantification).(metabolite).C(1);       % Store gray matter metabolite conc
            C_WM = MRSCont.GlobalConc.(quantification).(metabolite).C(2);       % Store white matter metabolite conc
            
            Q_predicted = W * MRSCont.GlobalConc.(quantification).(metabolite).C;   % Make linear predictions of metabolites
            residuals = Q - Q_predicted;                                            % Residuals 
            MRSCont.GlobalConc.(quantification).(metabolite).residuals = residuals; % Store residuals
            
            Q_predicted_full = nan(size(Q_full));                               % Create predicted metabolite vector with full size
            residuals_full = nan(size(Q_full));                                 % Create residual vector with full size
            
            Q_predicted_full(valid_mask) = Q_predicted;                     % Fill in the predicted metabolites
            residuals_full(valid_mask) = residuals;                         % Fill in the residuals
            
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

            MRSCont.GlobalConc.(quantification).(metabolite).Q_predicted_image = Q_predicted_image; % Store predicted metabolite image
            MRSCont.GlobalConc.(quantification).(metabolite).residuals_image = residuals_image;     % Store resiudal image
            MRSCont.GlobalConc.(quantification).(metabolite).valid_mask_image = valid_mask_image;   % Store mask 
            
            % Do some statisitcs
            SS_res = sum(residuals.^2);                                     % SSQ
            SS_tot = sum((Q - mean(Q)).^2);
            R2 = 1 - SS_res/SS_tot;                                         % R2
    
            sigma = sqrt(SS_res / (length(Q) - 2));
            WtW = W' * W;
            cov_matrix = sigma^2 * inv(WtW);
    
            SE_C_GM = sqrt(cov_matrix(1,1));                                % Gray matter conc standard error
            SE_C_WM = sqrt(cov_matrix(2,2));                                % White matter conc standard error
            MRSCont.GlobalConc.(quantification).(metabolite).SE_C_GM = SE_C_GM;     % Store standard error of gray matter global metabolite conc
            MRSCont.GlobalConc.(quantification).(metabolite).SE_C_WM = SE_C_WM;     % Store standard error of white matter global metabolite conc
            
            % 95% confidence interval
            alpha = 0.05;
            t_critical = tinv(1 - alpha/2, length(Q) - 2);                          
            
            CI_C_GM = [C_GM - t_critical * SE_C_GM, C_GM + t_critical * SE_C_GM];
            CI_C_WM = [C_WM - t_critical * SE_C_WM, C_WM + t_critical * SE_C_WM];
            MRSCont.GlobalConc.(quantification).(metabolite).CI_C_GM = CI_C_GM;
            MRSCont.GlobalConc.(quantification).(metabolite).CI_C_WM = CI_C_WM;

        end
    
        residual_min = min(residuals,[],'all');                             % Residual minimum for plots
        Q_max = max(Q,[],'all');                                            % Metabolite maximum for plots
        MRSCont.GlobalConc.(quantification).(metabolite).residual_min = residual_min; % Store residual minimum
        MRSCont.GlobalConc.(quantification).(metabolite).Q_max = Q_max;               % Store metabolite maximum
    end

    

end