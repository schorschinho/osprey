function out = osp_plotMeanSpec(MRSCont, which_spec, g, shift,group, xlab, ylab, figTitle,basis,exp)
%% out = osp_plotMeanSpec(MRSCont, which_spec,g, shift, xlab, ylab, figTitle)
%   Creates a figure mean and standard deviation of the spectra. If the
%   chosen spectra was fitted the mean fit, baseline and residue are shown.
%
%   USAGE:
%       out = osp_plotMeanSpec(MRSCont, which,GUI, g, shift,group, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       which    = String for the spectrum to plot (optional)
%                   OPTIONS:    'A' (default)
%                               'B'
%                               'C'
%                               'D'
%                               'diff1'
%                               'diff2'
%                               'sum'
%                               'ref'
%                               'w'
%       g         = group index
%       shift     = shift between groups
%       group     = control flag to decide whether different groups are plotted (removes residual and baseline)
%       xlab      = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab      = label for the y-axis (optional.  Default = '');
%       figTitle  = label for the title of the plot (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Helge Zoellner (Johns Hopkins University, 2019-10-02)
%       hzoelln2@jhmi.edu
%
%   HISTORY:
%       2020-01-24: First version of the code.

% Check that OspreyOverview has been run before
if ~MRSCont.flags.didOverview
    error('Trying to create overview plots, but no overview data has been created. Run OspreyOverview first.')
end

cb = cbrewer('qual', 'Dark2', 12, 'pchip');
temp = cb(3,:);
cb(3,:) = cb(4,:);
cb(4,:) = temp;

%%% 1. PARSE INPUT ARGUMENTS %%%
% Get the fit method and style
fitMethod   = MRSCont.opts.fit.method;
fitStyle    = MRSCont.opts.fit.style;

% Fall back to defaults if not provided
if nargin<10
    exp = 1;
    if nargin<9
        basis = 1;
        if nargin<8
            figTitle = '';
            if nargin<7
            ylab='';
                if nargin<6
                    xlab='Frequency (ppm)';
                    if nargin<5
                        group = 0;
                        if nargin<4
                            shift = 0.1;
                            if nargin<3
                                g = 1;
                                if nargin < 2
                                    which_spec = 'A';
                                    if nargin<1
                                        error('ERROR: no input Osprey container specified.  Aborting!!');
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if strcmp(which_spec,'MM') % re_mm
    which_spec = 'mm'; % re_mm
end % re_mm

% Check version of Osprey - since we have changed the layout of the Overview struct with the implementation of DualVoxel
if isfield(MRSCont.overview.Osprey, 'sort_data')
    sort_data = 'sort_data';
    sort_fit = 'sort_fit';
else
    sort_data = 'sort_data_voxel_1';
    sort_fit = 'sort_models_voxel_1';
end

% Create a fitMethod-specific theme
switch fitMethod
    case 'Osprey'
        colorFit  = MRSCont.colormap.Accent;
    case 'Osprey_gLCM'
        colorFit  = [255/255 140/255 0/255];
    case 'LCModel'
        colorFit  = 'r';
end

%%% 2. EXTRACT DATA TO PLOT %%%
% Extract normalized spectra and fits
if isnumeric(g)
   if ~(isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM)
        shift = shift * (g-1);
        GroupString = ['g_' num2str(g)];
   else
       shift = shift * (g);
        GroupString = ['g_' num2str(g)];
   end
else
    g = 1;
    GroupString = 'GMean';
    shift = 0;
end

which_spec_split = split(which_spec);
if length(which_spec_split) == 3
    fit = which_spec_split{1};
    spec = which_spec_split{2};
    subspec = which_spec_split{3};
else
    fit = [];
    spec = which_spec_split{1};
    subspec = which_spec_split{2};
end

if MRSCont.flags.didFit
    if isempty(fit) %Is data
        ind = find(strcmp(MRSCont.overview.SubSpecNamesStruct.(spec),subspec));
        if ~(strcmp(spec,'ref') || strcmp(spec,'w') || strcmp(spec,'mm_ref'))
            data_mean = MRSCont.overview.Osprey.(sort_data).(GroupString).(['mean_' spec])(:,ind)';
            if size(MRSCont.overview.Osprey.(sort_data).(GroupString).(['sd_' spec]),2)<=1
                data_sd = [];
            else
                data_sd = MRSCont.overview.Osprey.(sort_data).(GroupString).(['sd_' spec])(:,ind)';
            end
        else
            data_mean = MRSCont.overview.Osprey.(sort_data).(GroupString).(['mean_' spec])(:)';
            data_sd = MRSCont.overview.Osprey.(sort_data).(GroupString).(['sd_' spec])(:)';
        end
        ppm = MRSCont.overview.Osprey.(['ppm_data_' spec]);  
    else %Is fit
        switch spec
            case {'metab','mm'}
                try
                    ind = find(strcmp(MRSCont.overview.FitSpecNamesStruct.(spec)(basis,:,exp),subspec)); 
                catch
                    ind = find(strcmp(MRSCont.overview.FitSpecNamesStruct.(spec)(basis,:,1),subspec)); 
                end
                fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_fit_' spec])(basis,:,ind,exp);
                fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_fit_' spec])(basis,:,ind,exp);
                data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' spec])(basis,:,ind,exp);
                data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' spec])(basis,:,ind,exp);
                baseline_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_baseline_' spec])(basis,:,ind,exp);
                baseline_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_baseline_' spec])(basis,:,ind,exp);
                residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' spec])(basis,:,ind,exp);
                residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' spec])(basis,:,ind,exp);
                ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec])(basis,:,ind,exp);
                if strcmp(spec,'mm')
                    MM_clean_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_MM_clean_' spec])(basis,:,ind,exp);
                    MM_clean_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_MM_clean_' spec])(basis,:,ind,exp);
                end  
                if MRSCont.opts.fit.fitMM
                    MM_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_fittMM_' spec])(basis,:,ind,exp);
                    MM_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_fittMM_' spec])(basis,:,ind,exp);
                end
        case {'ref','w'}
                if ~strcmp(MRSCont.opts.fit.method, 'LCModel')
                    fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_fit_' spec])(1,:,1,exp);
                    fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_fit_' spec])(1,:,1,exp);
                    data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' spec])(1,:,1,exp);
                    data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' spec])(1,:,1,exp);
                    residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' spec])(1,:,1,exp);
                    residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' spec])(1,:,1,exp);
                    ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec])(1,:,1,exp);
                else
                    data_mean = MRSCont.overview.Osprey.(sort_data).(GroupString).(['mean_' spec])(:,1)';
                    data_sd = MRSCont.overview.Osprey.(sort_data).(GroupString).(['sd_' spec])(:,1)';
                    ppm = MRSCont.overview.Osprey.(['ppm_data_' spec]);                        
                end
        end
    end    
else
    ind = find(strcmp(MRSCont.overview.SubSpecNamesStruct.(spec),subspec));
    if ~(strcmp(spec,'ref') || strcmp(spec,'w') || strcmp(spec,'mm_ref'))
        data_mean = MRSCont.overview.Osprey.(sort_data).(GroupString).(['mean_' spec])(:,ind)';
        if size(MRSCont.overview.Osprey.(sort_data).(GroupString).(['sd_' spec]),2)<=1
            data_sd = [];
        else
            data_sd = MRSCont.overview.Osprey.(sort_data).(GroupString).(['sd_' spec])(:,ind)';
        end
    else
        data_mean = MRSCont.overview.Osprey.(sort_data).(GroupString).(['mean_' spec])(:)';
        data_sd = MRSCont.overview.Osprey.(sort_data).(GroupString).(['sd_' spec])(:)';
    end
    ppm = MRSCont.overview.Osprey.(['ppm_data_' spec]);            
end

[~,min_ppm_index] = min(ppm);

if ~exist('name', 'var')
    name = spec;
end

if ~exist('ppmRange', 'var')
    if ~strcmp(GroupString,'GMean')
        if (~strcmp(spec,'w') && ~strcmp(spec,'ref') && ~strcmp(spec,'mm_ref'))
            if exist('fit_mean', 'var')
                figTitle = ['mean data \pm SD & mean model: ' spec ' ' subspec];
            else
                figTitle = ['mean \pm  SD: ' name];
            end
            ppmRange = MRSCont.opts.fit.range;
        else
           figTitle = ['mean data \pm SD & mean model: ' spec ' ' subspec];
           ppmRange = MRSCont.opts.fit.rangeWater;
        end
    else
        if (~strcmp(spec,'w') && ~strcmp(spec,'ref') && ~strcmp(spec,'mm_ref'))
            if exist('fit_mean', 'var')
                figTitle = ['Grand mean data \pm SD &  Grand mean model: ' spec ' ' subspec];
            else
                figTitle = ['Grand mean \pm  SD: ' name];
            end
            ppmRange = MRSCont.opts.fit.range;
        else
           figTitle = ['Grand mean data \pm SD & Grand mean model: ' spec ' ' subspec];
           ppmRange = MRSCont.opts.fit.rangeWater;
        end        
    end
end

%Calculate SD shadows
if length(data_sd) > 1
    data_yu = data_mean + data_sd;
    data_yl = data_mean - data_sd;
end
    
if exist('fit_mean', 'var') && length(fit_sd) > 1
    fit_yu = fit_mean + fit_sd;
    fit_yl = fit_mean - fit_sd;
end
if exist('baseline_mean', 'var') && length(baseline_sd) > 1
    baseline_yu = baseline_mean + baseline_sd;
    baseline_yl = baseline_mean - baseline_sd;
end
if exist('residual_mean', 'var') && length(residual_sd) > 1
    residual_yu = residual_mean + residual_sd;
    residual_yl = residual_mean - residual_sd;
end
if exist('MM_clean_sd', 'var') && length(MM_clean_sd) > 1
    MM_clean_yu = MM_clean_mean + MM_clean_sd;
    MM_clean_yl = MM_clean_mean - MM_clean_sd;
end
if exist('MM_mean', 'var') && length(MM_sd) > 1 && ~isnan(MM_mean(1))
    MM_yu = MM_mean + MM_sd;
    MM_yl = MM_mean - MM_sd;
end
if exist('data_yu', 'var')
    maxshift = max(data_yu);
    maxshift_abs = max(abs(data_yu));
else
    maxshift = max(data_mean);
    maxshift_abs = max(abs(data_mean));    
end
shift = maxshift_abs * shift;
%%% 3. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
out = figure('Visible','on');
hold on

%%% 4. PLOT DATA, FIT, RESIDUAL, BASELINE %%%
% Determine a positive stagger to offset data, fit, residual, and
% baseline from the individual metabolite contributions
if MRSCont.flags.isGUI
    if exist('residual_mean', 'var')
        if ~group && length(residual_sd) > 1
            fill([ppm fliplr(ppm)], [(residual_yu+shift+ max(maxshift +  abs(min(residual_mean)))) (fliplr(residual_yl)+shift+ max(maxshift +  abs(min(residual_mean))))], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow res
        end
    end
    if exist('data_yu', 'var')
        fill([ppm fliplr(ppm)], [data_yu+shift fliplr(data_yl)+shift], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow data
    end

    if exist('fit_mean', 'var')
        plot(ppm,fit_mean+shift ,'color', colorFit, 'LineWidth', 1.5); %Fit
        if ~group
            plot(ppm,residual_mean+shift+ max(maxshift +  abs(min(residual_mean))) ,'color', MRSCont.colormap.Foreground, 'LineWidth', 1);  %Residual
        else
            plot(ppm,residual_mean+shift-maxshift_abs*0.3 ,'color', cb(g,:), 'LineWidth', 1);  %Residual
        end
    end

    if exist('baseline_mean', 'var')
        plot(ppm,baseline_mean+shift ,'color', MRSCont.colormap.LightAccent, 'LineWidth', 1); %Baseline
    end

    if exist('MM_mean', 'var') && ~isnan(MM_mean(1))
        plot(ppm,MM_mean+baseline_mean+shift ,'color', colorFit, 'LineWidth', 1); %MM Baseline
    end

    plot(ppm,data_mean+shift ,'color',cb(g,:), 'LineWidth', 2); % Data
    
    if strcmp(which_spec,'mm') % re_mm 
        if exist('MM_clean_yu', 'var')
            fill([ppm fliplr(ppm)], [MM_clean_yu+maxshift_abs*1.2 fliplr(MM_clean_yl)+maxshift_abs*1.2], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow cleaned MM data
        end
        plot(ppm,MM_clean_mean+maxshift_abs*1.2 ,'color',[1 0 0.1], 'LineWidth', 2); % Data cleaned for MM data       
    end


    if exist('fit_mean', 'var')
        if ~group
            plot(ppm, (zeros(1,length(ppm))), 'Color',MRSCont.colormap.Foreground); % Zeroline
            text(ppm(min_ppm_index)-0.05, 0, '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline text
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift)), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Maximum Data
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean)))), 'Color',MRSCont.colormap.Foreground, 'LineStyle','--', 'LineWidth', 0.5); % Zeroline Residue
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Max Residue
            text(ppm(min_ppm_index)-0.05, 0 + max(maxshift +  abs(min(residual_mean))), '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline Residual text
            text(ppm(min_ppm_index)-0.05, (0 +max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))), [num2str(100/max(fit_mean)*abs(max(residual_mean)),'%10.1f') '%'], 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Max Residue text
        end
    end
else
    if exist('fit_mean', 'var')
        if ~group
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift)), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Maximum Data
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean)))), 'Color',MRSCont.colormap.Foreground, 'LineStyle','--', 'LineWidth', 0.5); % Zeroline Residue
            plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Max Residue
            text(ppm(min_ppm_index)-0.05, 0 + max(maxshift +  abs(min(residual_mean))), '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline Residual text
            text(ppm(min_ppm_index)-0.05, (0 +max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))), [num2str(100/max(fit_mean)*abs(max(residual_mean)),'%10.1f') '%'], 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Max Residue text
        end
    end
    if exist('residual_mean', 'var')
        if ~group && length(residual_sd) > 1
            fill([ppm fliplr(ppm)], [(residual_yu+shift+ max(maxshift +  abs(min(residual_mean)))) (fliplr(residual_yl)+shift+ max(maxshift +  abs(min(residual_mean))))], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow res
        end
    end
    if exist('data_yu', 'var')
        fill([ppm fliplr(ppm)], [data_yu+shift fliplr(data_yl)+shift], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow data
    end
    plot(ppm,data_mean+shift ,'color',cb(g,:), 'LineWidth', 2); % Data
    if exist('fit_mean', 'var')
        plot(ppm,fit_mean+shift ,'color', colorFit, 'LineWidth', 1); %Fit
        if ~group
            plot(ppm,residual_mean+shift+ max(maxshift +  abs(min(residual_mean))) ,'color', MRSCont.colormap.Foreground, 'LineWidth', 1);  %Residual
            plot(ppm, (zeros(1,length(ppm))), 'Color',MRSCont.colormap.Foreground); % Zeroline
            text(ppm(min_ppm_index)-0.05, 0, '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline text
        else
            plot(ppm,residual_mean+shift-maxshift_abs*0.3 ,'color', cb(g,:), 'LineWidth', 1);  %Residual
        end
    end
    if exist('baseline_mean', 'var')
            plot(ppm,baseline_mean+shift ,'color', MRSCont.colormap.LightAccent, 'LineWidth', 1); %Baseline
    end
    if exist('MM_mean', 'var')
        plot(ppm,MM_mean+baseline_mean+shift ,'color', colorFit, 'LineWidth', 1); %MM Baseline
    end
    
    if strcmp(which_spec,'mm') % re_mm 
        if exist('MM_clean_yu', 'var')
            fill([ppm fliplr(ppm)], [MM_clean_yu+maxshift_abs*1.2 fliplr(MM_clean_yl)+maxshift_abs*1.2], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow cleaned MM data
        end
        plot(ppm,MM_clean_mean+maxshift_abs*1.2 ,'color',[1 0 0.1], 'LineWidth', 2); % Data cleaned for MM data       
    end

end

%%% 5. DUAL VOXEL PLOTTING %%%%
%Add second voxle to the plot if necessary
if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM

    sort_data = 'sort_data_voxel_2';
    sort_fit = 'sort_fit_voxel_2';

    if MRSCont.flags.isUnEdited
        switch which_spec
            case 'A'
                fit = 'off';
                fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_' spec]);
                fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_' spec]);
                data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' spec]);
                data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' spec]);
                baseline_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_baseline_' spec]);
                baseline_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_baseline_' spec]);
                residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' spec]);
                residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' spec]);
                ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec]);
                if MRSCont.opts.fit.fitMM
                MM_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_fittMM_' spec]);
                MM_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_fittMM_' spec]);
                end
            case {'ref','w'}
                fit = which_spec;
                fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_' which_spec]);
                fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_which_' spec]);
                data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' which_spec]);
                data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' which_spec]);
                residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' which_spec]);
                residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' which_spec]);
                ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' which_spec]);
            case {'mm'}
    %             name = 'mm';
    %             data_mean = MRSCont.overview.Osprey.sort_data.(GroupString).(['mean_' which_spec]);
    %             data_sd = MRSCont.overview.Osprey.sort_data.(GroupString).(['sd_' which_spec]);
    %             ppm = MRSCont.overview.Osprey.(['ppm_data_' which_spec]);

                fit = 'mm';
                fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_' spec]);
                fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_' spec]);
                data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' spec]);
                data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' spec]);
                baseline_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_baseline_' spec]);
                baseline_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_baseline_' spec]);
                residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' spec]);
                residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' spec]);
                ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec]);
                MM_clean_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_MM_clean_' spec]);
                MM_clean_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_MM_clean_' spec]);
                if MRSCont.opts.fit.fitMM
                    MM_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_fittMM_' spec]);
                    MM_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_fittMM_' spec]);
                end
        end
    end
    if MRSCont.flags.isMEGA
        switch which_spec
            case 'A'
                name = 'off';
                if strcmp(fitStyle,'Concatenated')
                    data_mean = MRSCont.overview.Osprey.(sort_data).(GroupString).(['mean_' which_spec]);
                    data_sd = MRSCont.overview.Osprey.(sort_data).(GroupString).(['sd_' which_spec]);
                    ppm = MRSCont.overview.Osprey.(['ppm_data_' which_spec]);
                else
                    fit = 'off';
                    fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_' spec]);
                    fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_' spec]);
                    data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' spec]);
                    data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' spec]);
                    baseline_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_baseline_' spec]);
                    baseline_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_baseline_' spec]);
                    residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' spec]);
                    residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' spec]);
                    ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec]);
                end
            case 'B'
                name = 'on';
                data_mean = MRSCont.overview.Osprey.(sort_data).(GroupString).(['mean_' which_spec]);
                data_sd = MRSCont.overview.Osprey.(sort_data).(GroupString).(['sd_' which_spec]);
                ppm = MRSCont.overview.Osprey.(['ppm_data_' which_spec]);
            case 'diff1'
                if ~strcmp(fitStyle,'Concatenated')
                    fit = which_spec;
                else
                    fit = 'conc';
                end
                fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_' spec]);
                fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_' spec]);
                data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' spec]);
                data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' spec]);
                baseline_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_baseline_' spec]);
                baseline_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_baseline_' spec]);
                residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' spec]);
                residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' spec]);
                ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec]);
            case {'sum'}
                if ~strcmp(fitStyle,'Concatenated')
                    fit = which_spec;
                    data_mean = MRSCont.overview.Osprey.(sort_data).(GroupString).(['mean_' which_spec]);
                    data_sd = MRSCont.overview.Osprey.(sort_data).(GroupString).(['sd_' which_spec]);
                    ppm = MRSCont.overview.Osprey.(['ppm_data_' which_spec]);
                else
                    fit = 'conc';
                    fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_' spec]);
                    fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_' spec]);
                    data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' spec]);
                    data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' spec]);
                    baseline_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_baseline_' spec]);
                    baseline_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_baseline_' spec]);
                    residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' spec]);
                    residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' spec]);
                    ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec]);
                end
            case {'ref','w'}
                fit = which_spec;
                fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_' which_spec]);
                fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_' which_spec]);
                data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' which_spec]);
                data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' which_spec]);
                residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' which_spec]);
                residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' which_spec]);
                ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' which_spec]);
        end
    end
    if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
        switch which_spec
            case {'A','B','C','D'}
                data_mean = MRSCont.overview.Osprey.(sort_data).(GroupString).(['mean_' which_spec]);
                data_sd = MRSCont.overview.Osprey.(sort_data).(GroupString).(['sd_' which_spec]);
                ppm = MRSCont.overview.Osprey.(['ppm_data_' which_spec]);
            case {'diff1','diff2','sum'}
                if ~strcmp(fitStyle,'Concatenated')
                    fit = which_spec;
                    fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_' spec]);
                    fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_' spec]);
                    data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' spec]);
                    data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' spec]);
                    baseline_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_baseline_' spec]);
                    baseline_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_baseline_' spec]);
                    residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' spec]);
                    residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' spec]);
                    ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec]);
                else
                    fit = 'conc';
                    fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_' spec]);
                    fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_' spec]);
                    data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' spec]);
                    data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' spec]);
                    baseline_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_baseline_' spec]);
                    baseline_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_baseline_' spec]);
                    residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' spec]);
                    residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' spec]);
                    ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec]);
                end
            case {'ref','w'}
                fit = which_spec;
                fit_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_' which_spec]);
                fit_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_' which_spec]);
                data_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_data_' which_spec]);
                data_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_data_' which_spec]);
                residual_mean = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['mean_res_' which_spec]);
                residual_sd = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['sd_res_' which_spec]);
                ppm = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' which_spec]);
        end
    end

    if ~exist('name', 'var')
        name = which_spec;
    end

    if ~exist('ppmRange', 'var')
        if ~strcmp(GroupString,'GMean')
            if (~strcmp(which_spec,'w') && ~strcmp(which_spec,'ref'))
                if exist('fit_mean', 'var')
                    figTitle = ['mean data \pm SD & mean fit: ' which_spec];
                else
                    figTitle = ['mean \pm  SD: ' name];
                end
                ppmRange = MRSCont.opts.fit.range;
            else
               figTitle = ['mean data \pm SD & mean fit: ' which_spec];
               ppmRange = MRSCont.opts.fit.rangeWater;
            end
        else
            if (~strcmp(which_spec,'w') && ~strcmp(which_spec,'ref'))
                if exist('fit_mean', 'var')
                    figTitle = ['Grand mean data \pm SD &  Grand mean fit: ' which_spec];
                else
                    figTitle = ['Grand mean \pm  SD: ' name];
                end
                ppmRange = MRSCont.opts.fit.range;
            else
               figTitle = ['Grand mean data \pm SD & Grand mean fit: ' which_spec];
               ppmRange = MRSCont.opts.fit.rangeWater;
            end        
        end
    end

    %Calculate SD shadows
    if length(data_sd) > 1
        data_yu = data_mean + data_sd;
        data_yl = data_mean - data_sd;
    end



    if exist('fit_mean', 'var') && length(fit_sd) > 1
        fit_yu = fit_mean + fit_sd;
        fit_yl = fit_mean - fit_sd;
    end
    if exist('baseline_mean', 'var') && length(baseline_sd) > 1
        baseline_yu = baseline_mean + baseline_sd;
        baseline_yl = baseline_mean - baseline_sd;
    end
    if exist('residual_mean', 'var') && length(residual_sd) > 1
        residual_yu = residual_mean + residual_sd;
        residual_yl = residual_mean - residual_sd;
    end
    if exist('MM_clean_sd', 'var') && length(MM_clean_sd) > 1
        MM_clean_yu = MM_clean_mean + MM_clean_sd;
        MM_clean_yl = MM_clean_mean - MM_clean_sd;
    end
    if exist('MM_mean', 'var') && length(MM_sd) > 1
        MM_yu = MM_mean + MM_sd;
        MM_yl = MM_mean - MM_sd;
    end
    if exist('data_yu', 'var')
        maxshift = max(data_yu);
        maxshift_abs = max(abs(data_yu));
    else
        maxshift = max(data_mean);
        maxshift_abs = max(abs(data_mean));    
    end
    shift = maxshift_abs * shift + 0.15*maxshift_abs;


    % Determine a positive stagger to offset data, fit, residual, and
    % baseline from the individual metabolite contributions
    if MRSCont.flags.isGUI
        if exist('residual_mean', 'var')
            if ~group && length(residual_sd) > 1
                fill([ppm fliplr(ppm)], [(residual_yu+shift+ max(maxshift +  abs(min(residual_mean)))) (fliplr(residual_yl)+shift+ max(maxshift +  abs(min(residual_mean))))], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow res
            end
        end
        if exist('data_yu', 'var')
            fill([ppm fliplr(ppm)], [data_yu+shift fliplr(data_yl)+shift], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow data
        end

        if exist('fit_mean', 'var')
            plot(ppm,fit_mean+shift ,'color', colorFit, 'LineWidth', 1.5); %Fit
            if ~group
                plot(ppm,residual_mean+shift+ max(maxshift +  abs(min(residual_mean))),':' ,'color', MRSCont.colormap.Foreground, 'LineWidth', 1);  %Residual
            else
                plot(ppm,residual_mean+shift-maxshift_abs*0.3,':' ,'color', cb(g,:), 'LineWidth', 1);  %Residual
            end
        end

        if exist('baseline_mean', 'var')
            plot(ppm,baseline_mean+shift,':' ,'color', MRSCont.colormap.LightAccent, 'LineWidth', 1); %Baseline
        end

        if exist('MM_mean', 'var')
            plot(ppm,MM_mean+baseline_mean+shift,':'  ,'color', colorFit, 'LineWidth', 1); %MM Baseline
        end

        plot(ppm,data_mean+shift,':'  ,'color',cb(g,:), 'LineWidth', 2); % Data

        if strcmp(which_spec,'mm') % re_mm 
            if exist('MM_clean_yu', 'var')
                fill([ppm fliplr(ppm)], [MM_clean_yu+maxshift_abs*1.2 fliplr(MM_clean_yl)+maxshift_abs*1.2], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow cleaned MM data
            end
            plot(ppm,MM_clean_mean+maxshift_abs*1.2,':'  ,'color',[1 0 0.1], 'LineWidth', 2); % Data cleaned for MM data       
        end


        if exist('fit_mean', 'var')
            if ~group
                plot(ppm, (zeros(1,length(ppm))), 'Color',MRSCont.colormap.Foreground); % Zeroline
                text(ppm(min_ppm_index)-0.05, 0, '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline text
                plot(ppm, (zeros(1,length(ppm)) + max(maxshift)),':' , 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Maximum Data
                plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean)))),':' , 'Color',MRSCont.colormap.Foreground, 'LineStyle','--', 'LineWidth', 0.5); % Zeroline Residue
                plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))),':' , 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Max Residue
                text(ppm(min_ppm_index)-0.05, 0 + max(maxshift +  abs(min(residual_mean))), '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline Residual text
                text(ppm(min_ppm_index)-0.05, (0 +max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))), [num2str(100/max(fit_mean)*abs(max(residual_mean)),'%10.1f') '%'], 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Max Residue text
            end
        end
    else
        if exist('fit_mean', 'var')
            if ~group
                plot(ppm, (zeros(1,length(ppm)) + max(maxshift)), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Maximum Data
                plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean)))),':' , 'Color',MRSCont.colormap.Foreground, 'LineStyle','--', 'LineWidth', 0.5); % Zeroline Residue
                plot(ppm, (zeros(1,length(ppm)) + max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1); % Max Residue
                text(ppm(min_ppm_index)-0.05, 0 + max(maxshift +  abs(min(residual_mean))), '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline Residual text
                text(ppm(min_ppm_index)-0.05, (0 +max(maxshift +  abs(min(residual_mean))) + abs(max(residual_mean))), [num2str(100/max(fit_mean)*abs(max(residual_mean)),'%10.1f') '%'], 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Max Residue text
            end
        end
        if exist('residual_mean', 'var')
            if ~group && length(residual_sd) > 1
                fill([ppm fliplr(ppm)], [(residual_yu+shift+ max(maxshift +  abs(min(residual_mean)))) (fliplr(residual_yl)+shift+ max(maxshift +  abs(min(residual_mean))))], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow res
            end
        end
        if exist('data_yu', 'var')
            fill([ppm fliplr(ppm)], [data_yu+shift fliplr(data_yl)+shift], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow data
        end
        plot(ppm,data_mean+shift,':'  ,'color',cb(g,:), 'LineWidth', 2); % Data
        if exist('fit_mean', 'var')
            plot(ppm,fit_mean+shift,':'  ,'color', colorFit, 'LineWidth', 1); %Fit
            if ~group
                plot(ppm,residual_mean+shift+ max(maxshift +  abs(min(residual_mean))),':'  ,'color', MRSCont.colormap.Foreground, 'LineWidth', 1);  %Residual
                plot(ppm, (zeros(1,length(ppm))), 'Color',MRSCont.colormap.Foreground); % Zeroline
                text(ppm(min_ppm_index)-0.05, 0, '0', 'FontSize', 10,'Color',MRSCont.colormap.Foreground); %Zeroline text
            else
                plot(ppm,residual_mean+shift-maxshift_abs*0.3,':'  ,'color', cb(g,:), 'LineWidth', 1);  %Residual
            end
        end
        if exist('baseline_mean', 'var')
                plot(ppm,baseline_mean+shift,':'  ,'color', MRSCont.colormap.LightAccent, 'LineWidth', 1); %Baseline
        end
        if exist('MM_mean', 'var')
            plot(ppm,MM_mean+baseline_mean+shift,':'  ,'color', colorFit, 'LineWidth', 1); %MM Baseline
        end

        if strcmp(which_spec,'mm') % re_mm 
            if exist('MM_clean_yu', 'var')
                fill([ppm fliplr(ppm)], [MM_clean_yu+maxshift_abs*1.2 fliplr(MM_clean_yl)+maxshift_abs*1.2], [0 0 0],'FaceAlpha',0.15, 'linestyle', 'none'); %SD Shadow cleaned MM data
            end
            plot(ppm,MM_clean_mean+maxshift_abs*1.2,':'  ,'color',[1 0 0.1], 'LineWidth', 2); % Data cleaned for MM data       
        end

    end    
end

%%% 6. DESIGN FINETUNING %%%
% Adapt common style for all axes
set(gca, 'XDir', 'reverse', 'XLim', [ppmRange(1), ppmRange(end)],'XMinorTick','on');
ticks = get(gca,'XTick');
set(gca, 'XTick', unique(round(ticks)));
set(gca, 'LineWidth', 1, 'TickDir', 'out');
set(gca, 'FontSize', 16);
% If no y caption, remove y axis
if isempty(ylab)
        set(gca, 'YColor', MRSCont.colormap.Background);
        set(gca,'YTickLabel',{});
        set(gca,'YTick',{});
        % Dirtywhite axes, light gray background
        set(gca, 'XColor', MRSCont.colormap.Foreground);
        set(gca, 'Color', MRSCont.colormap.Background);
        set(gcf, 'Color', MRSCont.colormap.Background);
        title(figTitle, 'Interpreter', 'Tex', 'Color', MRSCont.colormap.Foreground);
else
    set(gca, 'YColor', MRSCont.colormap.Foreground);
end

box off;
xlabel(xlab, 'FontSize', 16);
ylabel(ylab, 'FontSize', 16);


%%% 6. ADD OSPREY LOGO %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end

end
