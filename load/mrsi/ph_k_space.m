function [corr_ph, ph_corr_coeff] = ph_k_space(to_phase, all_spec_wat, all_spec_lip, all_spec)
%sob, sob :(
move_spec = find(to_phase);
static_spec = find(to_phase == 0);

% Try all phases from 5 to 355 (avoiding 0 degrees)
%ph = 5:5:355;
ph = 0:180;
move_spec_w_ph_wat = all_spec_wat(find(to_phase),:);
%move_spec_w_ph_lip = all_spec_lip(find(to_phase),:);

all_corr_coef = zeros(length(move_spec),length(static_spec), length(ph));

for st_idx = 1:length(static_spec)
    for ms_idx = 1:length(move_spec)
        for ph_idx = 1:length(ph)
            % Find correlation coefficients with "good" transients
            this_corr_coef = corrcoef(real(squeeze(move_spec_w_ph_wat(ms_idx,:).*exp(1i*pi*ph(ph_idx)/180))),...
                                     real(squeeze(all_spec_wat(static_spec(st_idx),:))));
            this_corr_coef = this_corr_coef(2);
                                 
%             if ((move_spec(ms_idx) == 1 && static_spec(st_idx) == 2) ||...
%                 (move_spec(ms_idx) == 2 && static_spec(st_idx) == 1) ||...
%                 (move_spec(ms_idx) == 3 && static_spec(st_idx) == 4) ||...
%                 (move_spec(ms_idx) == 4 && static_spec(st_idx) == 3))
%                 this_corr_coef_lip = corrcoef(real(squeeze(move_spec_w_ph_lip(ms_idx,:).*exp(1i*pi*ph(ph_idx)/180))),...
%                                      real(squeeze(all_spec_lip(static_spec(st_idx),:))));
%                 this_corr_coef_lip = this_corr_coef_lip(2);
%                                  
%                 this_corr_coef = mean([this_corr_coef this_corr_coef_lip]);
%             end
                    
            all_corr_coef(ms_idx,st_idx,ph_idx) = this_corr_coef;
        end
    end
end

all_corr_coef = squeeze(mean(all_corr_coef, 2));

if length(move_spec) > 1
    [m,p] = max(all_corr_coef,[],2);
else
    [m,p] = max(all_corr_coef);
end

corr_ph = ph(p);
ph_corr_coeff = m';