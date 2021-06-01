function out = osp_plotProcessMRSI(MRSCont, kk, which_spec,which_slice, ppmmin, ppmmax,mask, lb,shift_factor,add_text,add_No_MoCo)
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
if nargin < 11
    add_No_MoCo = 0;
    if nargin < 10
    add_text = 0;
        if nargin < 9
            shift_factor = 2;
            if nargin < 8
                lb = 0;
                if nargin < 7
                    mask = 1;
                        if nargin < 6
                        which_slice = 1;
                        if nargin<5
                            switch which_spec
                                case {'A', 'B', 'C', 'D', 'diff1', 'diff2','diff3', 'sum','mm'}
                                    ppmmax = 5;
                                case {'ref', 'w'}
                                    ppmmax = 2*4.68;
                                otherwise
                                    error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
                            end
                            if nargin<4
                                switch which_spec
                                    case {'A', 'B', 'C', 'D', 'diff1', 'diff2','diff3', 'sum'}
                                        ppmmin = 0.2;
                                    case {'ref', 'w','mm'}
                                        ppmmin = 0;
                                    otherwise
                                        error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
                                end
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
end

figTitle = ['MRSI spectra: ' which_spec ' ' num2str(ppmmin) ' to ' num2str(ppmmax) ' ppm'];
% Set up colormaps
if isfield(MRSCont,'colormap')
    colormap = MRSCont.colormap;
    tintFactor = 0.50;
    colormap.ForegroundTint = [colormap.Foreground(1)+(1-colormap.Foreground(1))*tintFactor...
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
XVox = MRSCont.raw{kk}.nXvoxels;
YVox = MRSCont.raw{kk}.nYvoxels;  
procData=op_takeVoxel(MRSCont.processed.(which_spec){kk},[1, 1]);
procData     = op_freqrange(procData,ppmmin,ppmmax);
procDataMarixToPlot = zeros((procData.sz(1)+50) * XVox,YVox);
zeroDataMarixToPlot = procDataMarixToPlot;
if add_No_MoCo
    procNoMoCoDataMarixToPlot = procDataMarixToPlot;
end
ppmLineToPlot = [];
for x = 1 : XVox
    ppmLineToPlot = horzcat(ppmLineToPlot, fliplr(procData.ppm),ones(1,50)*nan);
end


if MRSCont.flags.hasRef || MRSCont.flags.hasWater
    if MRSCont.flags.hasRef
        Norm = MRSCont.plot.processed.ref.ContMax;
    else
        Norm = MRSCont.plot.processed.w.ContMax;
    end        
    yLim = [min(MRSCont.plot.processed.(which_spec).min) max(MRSCont.plot.processed.(which_spec).max)]/Norm;
else
    Norm = MRSCont.plot.processed.A.ContMax;
    yLim = [min(MRSCont.plot.processed.A.min) max(MRSCont.plot.processed.A.max)]/Norm;
end


for y = 1 : YVox    
    procDataLineToPlot = [];
    zeroDataLineToPlot = [];
    for x = 1 : XVox
         procData=op_takeVoxel(MRSCont.processed.(which_spec){kk},[x, y,which_slice]);
         if ~(lb == 0)
            [procData,~]=op_filter(procData,lb);
         end
         procData     = op_freqrange(procData,ppmmin,ppmmax);
         procDataLineToPlot = vertcat(procDataLineToPlot,real(procData.specs)/Norm,ones(50,1)*nan);
         zeroDataLineToPlot = vertcat(zeroDataLineToPlot,zeros(length(procData.specs),1),ones(50,1)*nan);
    end
    procDataMarixToPlot(:,y) = procDataLineToPlot;
    zeroDataMarixToPlot(:,y) = zeroDataLineToPlot;
end

if add_No_MoCo
    for y = 1 : YVox  
    procNoMoCoDataLineToPlot = [];
        for x = 1 : XVox
             procNoMoCoData=op_takeVoxel(MRSCont.processed_no_MoCo.(which_spec){kk},[x, y,which_slice]);
             if ~(lb == 0)
                [procNoMoCoData,~]=op_filter(procNoMoCoData,lb);
             end
             procData     = op_freqrange(procNoMoCoData,ppmmin,ppmmax);
             procNoMoCoDataLineToPlot = vertcat(procNoMoCoDataLineToPlot,real(procData.specs)/Norm,ones(50,1)*nan);
        end
        procNoMoCoDataMarixToPlot(:,y) = procNoMoCoDataLineToPlot;
    end
end



%%% 3. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
if ~MRSCont.flags.isGUI
    out = figure;
else
    out = figure('Visible','off');
end


shift = abs(yLim(1)) + abs(yLim(2));
shift = shift* shift_factor;
for y = 1 : YVox
    for x = 1 : XVox
        if mask
            if MRSCont.mask{kk}(y,x,which_slice)
                plot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),y)+shift*y,'Color',colormap.ForegroundTint);
                hold on
                if add_text
                    text((x-1)*(procData.sz(1)+50)+1, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y-shift/18, num2str(ppmmin), 'Color', colormap.ForegroundTint);
                    text(x*(procData.sz(1)+50)-50, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y-shift/18, num2str(ppmmax), 'Color', colormap.ForegroundTint);
                    text(x*(procData.sz(1)+50)-50-procData.sz(1)/2, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y-shift/18, num2str((ppmmax-ppmmin)/2), 'Color', colormap.ForegroundTint);
                    text(x*(procData.sz(1)+50)-50-(procData.sz(1)/4)*3, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y-shift/18, num2str((ppmmax-ppmmin)/4), 'Color', colormap.ForegroundTint);
                    text(x*(procData.sz(1)+50)-50-(procData.sz(1)/4), zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y-shift/18, num2str(((ppmmax-ppmmin)/4)*3), 'Color', colormap.ForegroundTint);

                    text((x-1)*(procData.sz(1)+50)+1, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y, 'I', 'Color', colormap.ForegroundTint);
                    text(x*(procData.sz(1)+50)-50, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y, 'I', 'Color', colormap.ForegroundTint);
                    text(x*(procData.sz(1)+50)-50-procData.sz(1)/2, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y, 'I', 'Color', colormap.ForegroundTint);
                    text(x*(procData.sz(1)+50)-50-(procData.sz(1)/4)*3, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y, 'I', 'Color', colormap.ForegroundTint);
                    text(x*(procData.sz(1)+50)-50-(procData.sz(1)/4), zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y, 'I', 'Color', colormap.ForegroundTint);
                    text(x*(procData.sz(1)+50)-50-(procData.sz(1)/8), zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y, 'I', 'Color', colormap.ForegroundTint,'FontSize',8);
                    text(x*(procData.sz(1)+50)-50-(procData.sz(1)/8)*3, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y, 'I', 'Color', colormap.ForegroundTint,'FontSize',8);
                    text(x*(procData.sz(1)+50)-50-(procData.sz(1)/8)*5, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y, 'I', 'Color', colormap.ForegroundTint,'FontSize',8);
                    text(x*(procData.sz(1)+50)-50-(procData.sz(1)/8)*7, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y, 'I', 'Color', colormap.ForegroundTint,'FontSize',8);
                end
                 if add_No_MoCo
                    plot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),procNoMoCoDataMarixToPlot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),y)+shift*y,'Color',colormap.ForegroundTint);
                 end
                plot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),procDataMarixToPlot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),y)+shift*y,'Color',colormap.Foreground);                                   
               
                
            else
%                 plot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),procDataMarixToPlot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),y)+shift*y,'Color',colormap.ForegroundTint)
%                 hold on
%                 plot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),y)+shift*y,'Color',colormap.ForegroundTint);              
                
            end
        else
            plot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),procDataMarixToPlot((x-1)*(procData.sz(1)+50)+1:x*(procData.sz(1)+50),y)+shift*y,'Color',colormap.Foreground);
            if add_text
            text((x-1)*(procData.sz(1)+50)+1, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y-shift/15, num2str(ppmmin), 'Color', colormap.ForegroundTint);
                text(x*(procData.sz(1)+50)-50, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y-shift/15, num2str(ppmmax), 'Color', colormap.ForegroundTint);
                text(x*(procData.sz(1)+50)-50-procData.sz(1)/2, zeroDataMarixToPlot((x-1)*(procData.sz(1)+50)+1,y)+shift*y-shift/15, num2str((ppmmax-ppmmin)/2), 'Color', colormap.ForegroundTint);
            end
            hold on
        end
    end
end

hold off



%%% 8. DESIGN FINETUNING %%%
% Adapt common style for all axes

    

gcf = out;
set(gcf, 'Color', MRSCont.colormap.Background);        
box off;
set(gca, 'XDir', 'reverse');


if ~MRSCont.flags.isGUI
    set(gca, 'YColor', 'w');
    % Black axes, white background
    set(gca, 'XColor', 'w');
    set(gca, 'Color', 'w');
    set(gcf, 'Color', 'w');
    title(figTitle, 'Interpreter', 'none');
else
    set(gca, 'YColor', MRSCont.colormap.Background);
    set(gca,'YTickLabel',{});
    set(gca,'YTick',{});
    set(gca,'XTickLabel',{});
    set(gca,'XTick',{});
    set(gca, 'XColor', MRSCont.colormap.Background);
    set(gca, 'Color', MRSCont.colormap.Background);
    set(gcf, 'Color', MRSCont.colormap.Background);
    title(figTitle, 'Interpreter', 'none', 'Color', MRSCont.colormap.Foreground);
end

%%% 9. ADD OSPREY LOGO %%%
% Add to the printout, but not if displayed in the GUI.
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    text(gca, 0, -0.1, [MRSCont.ver.Osp],'Color', colormap.Foreground);
    imshow(I, map);
    axis off;
end
end

   