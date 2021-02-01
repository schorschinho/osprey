function out = osp_plotOverviewSpec(MRSCont, which_spec, g, shift, xlab, ylab, figTitle)
%% out = osp_plotOverviewSpec(MRSCont, which_spec, g, shift, xlab, ylab, figTitle)
%   Creates a figure with all spectra or fits.
%
%   USAGE:
%       out = osp_plotOverviewSpec, which,GUI, g, shift, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       which    = String for the spectrum or fit to plot (optional)
%                   OPTIONS:    'A' (default)
%                               'B'
%                               'C'
%                               'D'
%                               'diff1'
%                               'diff2'
%                               'sum'
%                               'mm'
%                               'ref'
%                               'w'
%                               'Fit:off'
%                               'Fit:diff1'
%                               'Fit:diff2'
%                               'Fit:sum'
%                               'Fit:mm'
%                               'Fit:ref'
%                               'Fit:w'
%                               'MM_clean'
%       g         = group index
%       shift     = shift between groups
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
fitStyle    = MRSCont.opts.fit.style;
% Fall back to defaults if not provided
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

 if strcmp(which_spec,'MM') %re_mm
    which_spec = 'mm'; %re_mm
 end %re_mm

% Check version of Osprey - since we have changed the layout of the Overview struct with the implementation of DualVoxel
if isfield(MRSCont.overview.Osprey, 'sort_data')
    sort_data = 'sort_data';
    sort_fit = 'sort_fit';
    mm = 'all_models';
else
    sort_data = 'sort_data_voxel_1';
    sort_fit = 'sort_fit_voxel_1';
    mm = 'all_models_voxel_1';
end
 
%%% 2. EXTRACT DATA TO PLOT %%%
% Extract normalized spectra and fits
if isnumeric(g)
    shift = shift * (g-1);
    GroupString = ['g_' num2str(g)];
else
    g = 1;
    GroupString = 'GMean';
    shift = 0;
end

%Is spectrum?

if (strcmp(which_spec,'A') || strcmp(which_spec,'B') || strcmp(which_spec,'C') || strcmp(which_spec,'D') || strcmp(which_spec,'diff1') || strcmp(which_spec,'diff2') || strcmp(which_spec,'sum') ||strcmp(which_spec,'mm') ||  strcmp(which_spec,'ref') || strcmp(which_spec,'w'))
    data = MRSCont.overview.Osprey.(sort_data).(GroupString).(which_spec);
    if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
        data2 = MRSCont.overview.Osprey.sort_data_voxel_2.(GroupString).(which_spec);
    end
    if nargin<8
        if (~strcmp(which_spec,'w') && ~strcmp(which_spec,'ref'))
            figTitle = ['Individual specs: ' which_spec];
            ppmRange = MRSCont.opts.fit.range;
        else
           figTitle = ['Individual specs: ' which_spec];
           ppmRange = MRSCont.opts.fit.rangeWater;
        end
    end
else % Is fit?
    fitwhich = which_spec(5:end);
    if MRSCont.flags.isUnEdited %Is UnEdited
        switch fitwhich
            case 'off'
                fit = 'A';
                data = MRSCont.overview.Osprey.(sort_fit).(GroupString).([fitwhich '_' fit]);
                    if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        data2 = MRSCont.overview.Osprey.sort_fit_voxel_2.(GroupString).([fitwhich '_' fit]);
                    end
            case {'ref','w','mm'}
                fit = fitwhich;
                data = MRSCont.overview.Osprey.(sort_fit).(GroupString).([fitwhich '_' fit]);
                if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    data2 = MRSCont.overview.Osprey.sort_fit_voxel_2.(GroupString).([fitwhich '_' fit]);
                end

            case {'lean'}
                fit = fitwhich;
                data = MRSCont.overview.Osprey.(mm).mm_mm;
                if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    data2 = MRSCont.overview.Osprey.(mm).mm_mm;
                end
        end
    end
    if MRSCont.flags.isMEGA %Is MEGA
        switch fitwhich
            case {'diff1','sum'}
                if strcmp(fitStyle,'Concatenated') %Is Concatenated?
                    fit = 'conc';
                    data = MRSCont.overview.Osprey.(sort_fit).(GroupString).([fit '_' fitwhich]);
                    if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        data2 = MRSCont.overview.Osprey.sort_fit_voxel_2.(GroupString).([fit '_' fitwhich]);
                    end
                else
                    fit = 'diff1';
                    data = MRSCont.overview.Osprey.(sort_fit).(GroupString).([fitwhich '_' fit]);
                    if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        data2 = MRSCont.overview.Osprey.sort_fit_voxel_2.(GroupString).([fitwhich '_' fit]);
                    end
                end
            case {'off'}
                    fit = 'A';
                    data = MRSCont.overview.Osprey.(sort_fit).(GroupString).([fitwhich '_' fit]);
                    if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        data2 = MRSCont.overview.Osprey.sort_fit_voxel_2.(GroupString).([fitwhich '_' fit]);
                    end
            case {'ref','w','mm'}
                fit = fitwhich;
                data = MRSCont.overview.Osprey.(sort_fit).(GroupString).([fit '_' fitwhich]);
                if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    data2 = MRSCont.overview.Osprey.sort_fit_voxel_2.(GroupString).([fit '_' fitwhich]);
                    end
        end
    end
    if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES) %Is Multiplexed
        switch fitwhich
            case {'diff1','diff2','sum'}
                if strcmp(fitStyle,'Concatenated') %Is Concatenated?
                    fit = 'conc';
                    data = MRSCont.overview.Osprey.(sort_fit).(GroupString).([fit '_' fitwhich]);
                    if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        data2 = MRSCont.overview.Osprey.sort_fit_voxel_2.(GroupString).([fit '_' fitwhich]);
                    end
                else
                    fit = fitwhich;
                    data = MRSCont.overview.Osprey.(sort_fit).(GroupString).([fit '_' fitwhich]);
                    if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                        data2 = MRSCont.overview.Osprey.sort_fit_voxel_2.(GroupString).([fit '_' fitwhich]);
                    end
                end
            case {'ref','w'}
                fit = fitwhich;
                data = MRSCont.overview.Osprey.(sort_fit).(GroupString).([fitwhich '_' fit]);
                if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    data2 = MRSCont.overview.Osprey.sort_fit_voxel_2.(GroupString).([fitwhich '_' fit]);
                end
        end
    end

    if nargin<7
        if (~strcmp(fitwhich,'w') && ~strcmp(fitwhich,'ref'))
            figTitle = ['Individual fits: ' fit ' ' fitwhich];
            ppmRange = MRSCont.opts.fit.range;
        else
           figTitle = ['Individual fits: ' fit ' ' fitwhich];
           ppmRange = MRSCont.opts.fit.rangeWater;
        end
    end
end



%%% 3. PLOT DATA %%%
if length(which_spec)>4
    if ~strcmp(which_spec(1:4),'MM_c') %re_mm
        if ~strcmp(which_spec(1:4),'Fit:')
            maxshift_abs = max(abs(data{1}.specs));
            shift = maxshift_abs * shift;
            out = plot(data{1}.ppm,data{1}.specs+shift ,'color', cb(g,:), 'LineWidth', 1); %data
            hold on;
            for kk = 2 : length(data)
                plot(data{kk}.ppm,data{kk}.specs+shift ,'color', cb(g,:), 'LineWidth', 1); %data
            end
            if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                maxshift_abs = max(abs(data2{1}.specs));
                shift = maxshift_abs * shift + maxshift_abs*0.15;
                out = plot(data2{1}.ppm,data2{1}.specs+shift ,':','color', cb(g,:), 'LineWidth', 2); %data
                hold on;
                for kk = 2 : length(data2)
                    plot(data2{kk}.ppm,data2{kk}.specs+shift ,':','color', cb(g,:), 'LineWidth', 2); %data
                end    
            end
        else
            maxshift_abs = max(abs(data{1}.fit));
            shift = maxshift_abs * shift;
            out = plot(data{1}.ppm,data{1}.fit+shift ,'color', cb(g,:), 'LineWidth', 1); %Fit
            hold on;
            for kk = 2 : length(data)
                plot(data{kk}.ppm,data{kk}.fit+shift ,'color', cb(g,:), 'LineWidth', 1); %Fit
            end
            if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                maxshift_abs = max(abs(data2{1}.fit));
                shift = maxshift_abs * shift+ maxshift_abs*0.15;
                out = plot(data2{1}.ppm,data2{1}.fit+shift,':' ,'color', cb(g,:), 'LineWidth', 2); %Fit
                hold on;
                for kk = 2 : length(data2)
                    plot(data2{kk}.ppm,data2{kk}.fit+shift,':' ,'color', cb(g,:), 'LineWidth', 2); %Fit
                end
            end
        end
    else %re_mm
        maxshift_abs = max(abs(data{1}.MM_clean));
        shift = maxshift_abs * shift;
        out = plot(data{1}.ppm,data{1}.MM_clean+shift ,'color', cb(g,:), 'LineWidth', 1); %data
        hold on;
        for kk = 2 : length(data)
            plot(data{kk}.ppm,data{kk}.MM_clean+shift ,'color', cb(g,:), 'LineWidth', 1); %data
        end
        if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
            maxshift_abs = max(abs(data2{1}.MM_clean));
            shift = maxshift_abs * shift+ maxshift_abs*0.15;
            out = plot(data2{1}.ppm,data2{1}.MM_clean+shift,':' ,'color', cb(g,:), 'LineWidth', 2); %data
            hold on;
            for kk = 2 : length(data2)
                plot(data2{kk}.ppm,data2{kk}.MM_clean+shift,':' ,'color', cb(g,:), 'LineWidth', 2); %data
            end
        end
    end %re_mm
else
        maxshift_abs = max(abs(data{1}.specs));
        shift = maxshift_abs * shift;
        out = plot(data{1}.ppm,data{1}.specs+shift ,'color', cb(g,:), 'LineWidth', 1); %data
        hold on;
        for kk = 2 : length(data)
            plot(data{kk}.ppm,data{kk}.specs+shift ,'color', cb(g,:), 'LineWidth', 1); %data
        end
        if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
            maxshift_abs = max(abs(data2{1}.specs));
            shift = maxshift_abs * shift+ maxshift_abs*0.15;
            out = plot(data2{1}.ppm,data2{1}.specs+shift,':' ,'color', cb(g,:), 'LineWidth', 2); %data
            hold on;
            for kk = 2 : length(data2)
                plot(data2{kk}.ppm,data2{kk}.specs+shift,':' ,'color', cb(g,:), 'LineWidth', 2); %data
            end
        end
end
%%% 4. DESIGN FINETUNING %%%
% Adapt common style for all axes
set(gca, 'XDir', 'reverse', 'XLim', [ppmRange(1), ppmRange(end)], 'XMinorTick', 'On');
set(gca, 'LineWidth', 1, 'TickDir', 'out');
set(gca, 'FontSize', 16);
% If no y caption, remove y axis
if isempty(ylab)
    if ~MRSCont.flags.isGUI
        set(gca, 'YColor', 'w');
        % Black axes, white background
        set(gca, 'XColor', 'k');
        set(gca, 'Color', 'w');
        set(gcf, 'Color', 'w');
        title(figTitle, 'Interpreter', 'Tex');
    else
        set(gca, 'YColor', MRSCont.colormap.Background);
        set(gca,'YTickLabel',{});
        set(gca,'YTick',{});
        % Dirtywhite axes, light gray background
        set(gca, 'XColor', MRSCont.colormap.Foreground);
        set(gca, 'Color', MRSCont.colormap.Background);
        set(gcf, 'Color', MRSCont.colormap.Background);
        title(figTitle, 'Interpreter', 'Tex', 'Color', MRSCont.colormap.Foreground);
    end
else
    set(gca, 'YColor', 'k');
end

box off;
xlabel(xlab, 'FontSize', 16);
ylabel(ylab, 'FontSize', 16);


%%% 5. ADD OSPREY LOGO %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out.Parent.Parent, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end

end
