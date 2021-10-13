function out = osp_plotProcessMRSIMatrix(MRSCont, kk, which_spec,coord, ppmmin, ppmmax,lb,add_No_MoCo,which_slice,yLim)
%% out = osp_plotProcessMRSI(MRSCont, kk, which, ppmmin, ppmmax)
%   Creates a figure showing processed data stored in an Osprey data container,
%   ie in the raw fields. This function will display the *processed and
%   averaged* data, i.e. after spectral alignment, averaging, water removal,
%   and other processing steps carried out in OspreyProcess.
%
%   USAGE:
%       out = osp_plotProcess(MRSCont, kk, which, ppmmin, ppmmax, xlab, ylab)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       which    = String for the spectrum to fit (optional)
%                   OPTIONS:    'A' (default)
%                               'B' (for MEGA, HERMES, HERCULES)
%                               'C' (for HERMES, HERCULES)
%                               'D' (for HERMES, HERCULES)
%                               'diff1' (for MEGA, HERMES, HERCULES)
%                               'diff2' (for HERMES, HERCULES)
%                               'sum' (for MEGA, HERMES, HERCULES)
%                               'ref'
%                               'w'
%       VoxelIndex = Index for the Voxel
%       xlab     = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab     = label for the y-axis (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-10-02)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-10-02: First version of the code.

% Check that OspreyProcess has been run before
if ~MRSCont.flags.didProcess
    error('Trying to plot processed data, but data has not been processed yet. Run OspreyProcess first.')
end

%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
if nargin <10
    yLim = [];
    if nargin < 9
        which_slice =1;
        if nargin < 8
        add_No_MoCo = 0;
            if nargin < 7
            lb = 0;
                if nargin<6
                switch which_spec{1}
                case {'A', 'B', 'C', 'D', 'diff1', 'diff2','diff3', 'sum','mm'}
                ppmmax = 5;
                case {'ref', 'w'}
                ppmmax = 2*4.68;
                otherwise
                error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
                end
                    if nargin<5
                    switch which_spec{1}
                    case {'A', 'B', 'C', 'D', 'diff1', 'diff2','diff3', 'sum'}
                    ppmmin = 0.2;
                    case {'ref', 'w','mm'}
                    ppmmin = 0;
                    otherwise
                    error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
                    end
                        if nargin < 4
                        coord = [2 8; 2 8];
                            if nargin < 3
                            which_spec = 'A';
                                if nargin < 2
                                kk = 1;
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
% Set up colormaps
if isfield(MRSCont,'colormap')
    colormap = MRSCont.colormap;
    tintFactor = 0.50;
    colormap.ForegroundTint{1} = [colormap.Foreground(1)+(1-colormap.Foreground(1))*tintFactor...
                               colormap.Foreground(2)+(1-colormap.Foreground(2))*tintFactor...
                               colormap.Foreground(3)+(1-colormap.Foreground(3))*tintFactor ];
                           
    tintFactor = 0.30;
    colormap.ForegroundTint{2} = [colormap.Foreground(1)+(1-colormap.Foreground(1))*tintFactor...
                               colormap.Foreground(2)+(1-colormap.Foreground(2))*tintFactor...
                               colormap.Foreground(3)+(1-colormap.Foreground(3))*tintFactor ];
   
else
    colormap.Background     = [1 1 1];
    colormap.LightAccent    = [110/255 136/255 164/255];
    colormap.Foreground     = [0 0 0];
    colormap.Accent         = [11/255 71/255 111/255];
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract raw and processed spectra in the plot range


if MRSCont.flags.hasRef || MRSCont.flags.hasWater
    if MRSCont.flags.hasRef
        Norm = MRSCont.plot.processed.ref.ContMax;
    else
        Norm = MRSCont.plot.processed.w.ContMax;
    end
    if isempty(yLim)
        yLim = [min(MRSCont.plot.processed.(which_spec{1}).min) max(MRSCont.plot.processed.(which_spec{1}).max*2)]/Norm;
        if strcmp(which_spec{1},'diff1') || strcmp(which_spec{1},'diff2')
            yLim = [min(MRSCont.plot.processed.(which_spec{1}).min)/2 max(MRSCont.plot.processed.(which_spec{1}).max)]/Norm;
        end
    end
else
    Norm = MRSCont.plot.processed.A.ContMax;
    if isempty(yLim)
        yLim = [min(MRSCont.plot.processed.A.min) max(MRSCont.plot.processed.A.max)]/Norm;
    end
end


%%% 3. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized

plot_index = 1;

for y = coord(2,1) : coord(2,2)
    for x = coord(1,1) : coord(1,2)
        out = subaxis(coord(2,2)-coord(2,1)+1,coord(1,2)-coord(1,1)+1,plot_index, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0, 'MarginLeft', 0.075, 'MarginBottom', 0.1);
        hold on
        for ss = 1 : length(which_spec)
            if add_No_MoCo
                procNoMoCoData=op_takeVoxel(MRSCont.processed_no_MoCo.(which_spec{ss}){kk},[x y which_slice]);
                if ~(lb == 0)
                    [procNoMoCoData,~]=op_filter(procNoMoCoData,lb);
                end
                plot(procNoMoCoData.ppm,real(procNoMoCoData.specs)/Norm,'Color',colormap.ForegroundTint{1});
            end
            procData=op_takeVoxel(MRSCont.processed.(which_spec{ss}){kk},[x y which_slice]);
            if ~(lb == 0)
                    [procData,~]=op_filter(procData,lb);
            end
            if ss == 1
                plot(procData.ppm,real(procData.specs)/Norm,'Color',colormap.Foreground);    
            else
                plot(procData.ppm,real(procData.specs)/Norm,'Color',colormap.ForegroundTint{2});   
            end
        end
        
        
        plot_index = plot_index +1;    
        box off;
        set(gca, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLim);
        set(gca, 'LineWidth', 1, 'TickDir', 'in', 'XMinorTick', 'Off','TickLength', [0.03 0.025]);
        set(gca, 'XColor', colormap.Foreground, 'YColor', colormap.Foreground);
        
        if x ~= coord(1,1)
            set(gca,'YTickLabel',{})
        end
        if y ~= coord(2,2)
            set(gca,'XTickLabel',{})
        else
            xlabel('Frequency (ppm)', 'Color', colormap.Foreground)
        end
    end
end

set(gcf, 'Color', MRSCont.colormap.Background);   

end

   