function out = osp_plotOverviewSpec(MRSCont, which_spec, g, shift, xlab, ylab, figTitle, basis)
%% out = osp_plotOverviewSpec(MRSCont, which_spec, g, shift, xlab, ylab, figTitle, basis)
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
if nargin<8 
    basis =1;
    if nargin<7
        figTitle=[];
        if nargin<6
            ylab='';
            if nargin<5
                xlab='Frequency (ppm)';
                if nargin<4
                    shift = 0.1;
                    if nargin<3
                        g = 1;
                        if nargin < 2
                            which_spec = 'metab A';
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


% Check version of Osprey - since we have changed the layout of the Overview struct with the implementation of DualVoxel
if isfield(MRSCont.overview.Osprey, 'sort_data')
    sort_data = 'sort_data';
    sort_fit = 'sort_fit';
else
    sort_data = 'sort_data_voxel_1';
    sort_fit = 'sort_models_voxel_1';
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

        
if isempty(fit) %Is data
    if ~(strcmp(spec,'ref') || strcmp(spec,'w') || strcmp(spec,'mm_ref'))
        ind = find(strcmp(MRSCont.overview.SubSpecNamesStruct.(spec),subspec));  
        ppmRange = MRSCont.opts.fit.range; 
    else
        ind =1;
        ppmRange = MRSCont.opts.fit.rangeWater;
    end
    data = MRSCont.overview.Osprey.(sort_data).(GroupString).(spec);
    if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
        data2 = MRSCont.overview.Osprey.sort_data_voxel_2.(GroupString).(spec);
    end
    if isempty(figTitle)
        figTitle = ['Individual specs: ' spec ' ' subspec];
    end
else % Is fit
    switch spec
        case {'metab','mm'}
            ind = find(strcmp(MRSCont.overview.FitSpecNamesStruct.(spec)(basis,:),subspec));   
            data = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['fit_' spec]);
                if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                    data2 = MRSCont.overview.Osprey.sort_fit_voxel_2.(GroupString).(spec);
                end
                ppm=MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec]);
                ppmRange = MRSCont.opts.fit.range;
        case {'ref','w','mm_ref'}
            data = MRSCont.overview.Osprey.(sort_fit).(GroupString).(['fit_' spec]);
            if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                data2 = MRSCont.overview.Osprey.sort_fit_voxel_2.(GroupString).(spec);
            end
            basis = 1;
            ind = 1;
            ppm=MRSCont.overview.Osprey.(sort_fit).(GroupString).(['ppm_fit_' spec]);
            ppmRange = MRSCont.opts.fit.rangeWater;
    end        
    if isempty(figTitle)
       figTitle = ['Individual fits: ' spec ' ' fit];
    end
end

out = figure;

%%% 3. PLOT DATA %%%
if ~isempty(fit) %Is fit
            maxshift_abs = max(abs(data(:,:,basis,ind)));
            shift = maxshift_abs * shift;
            
            plot(ppm(basis,:,ind),data(1,:,basis,ind)+shift ,'color', cb(g,:), 'LineWidth', 1); %Fit
            hold on;
            if size(data,1) > 1
                for kk = 2 : size(data,1)
                    plot(ppm(basis,:,ind),data(kk,:,basis,ind)+shift ,'color', cb(g,:), 'LineWidth', 1); %Fit
                end
            end
            if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                maxshift_abs = max(abs(data2{1}.fit));
                shift = maxshift_abs * shift+ maxshift_abs*0.15;
                plot(data2{1}.ppm,data2{1}.fit+shift,':' ,'color', cb(g,:), 'LineWidth', 2); %Fit
                hold on;
                if size(data2,1) > 1
                    for kk = 2 : length(data2)
                        plot(data2{kk}.ppm,data2{kk}.fit+shift,':' ,'color', cb(g,:), 'LineWidth', 2); %Fit
                    end
                end
            end
        else
            maxshift_abs = max(abs(data{1}.specs(:,ind)));
            shift = maxshift_abs * shift;
            plot(data{1}.ppm,data{1}.specs(:,ind)+shift ,'color', cb(g,:), 'LineWidth', 1); %data
            hold on;
            for kk = 2 : length(data)
                plot(data{kk}.ppm,data{kk}.specs(:,ind)+shift ,'color', cb(g,:), 'LineWidth', 1); %data
            end
            if isfield(MRSCont.flags,'isPRIAM')  && MRSCont.flags.isPRIAM
                maxshift_abs = max(abs(data2{1}.specs));
                shift = maxshift_abs * shift + maxshift_abs*0.15;
                plot(data2{1}.ppm,data2{1}.specs+shift ,':','color', cb(g,:), 'LineWidth', 2); %data
                hold on;
                for kk = 2 : length(data2)
                    plot(data2{kk}.ppm,data2{kk}.specs+shift ,':','color', cb(g,:), 'LineWidth', 2); %data
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
