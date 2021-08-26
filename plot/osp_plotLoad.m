function out = osp_plotLoad(MRSCont, kk, which, VoxelIndex, stag, ppmmin, ppmmax, xlab, ylab, figTitle)
%% out = osp_plotLoad(MRSCont, kk, which, VoxelIndex, stag, ppmmin, ppmmax, xlab, ylab, figTitle)
%   Creates a figure showing raw data stored in an Osprey data container,
%   ie in the raw fields. This function will display the *unprocessed*
%   data, i.e. all averages will be shown prior to spectral alignment,
%   averaging, water removal, and other processing steps carried out in
%   OspreyProcess.
%
%   USAGE:
%       out = osp_plotLoad(MRSCont, kk, which, stag, ppmmin, ppmmax, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset
%       which    = String for the spectrum to fit
%                   OPTIONS:    'mets'
%                                       'mm'
%                               'ref'
%                               'w'
%       VoxelIndex = Index for the Voxel
%       stag     = Numerical value representing the fraction of the maximum
%                   by which individual averages should be plotted
%                   vertically staggered (optional. Default - 0, i.e. no stagger)
%       xlab     = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab     = label for the y-axis (optional.  Default = '');
%       figTitle = label for the title of the plot (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-10-02)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-10-02: First version of the code.

% Check that OspreyLoad has been run before
if ~MRSCont.flags.didLoadData
    error('Trying to plot raw data, but no data has been loaded yet. Run OspreyLoad first.')
end


%%% 1. PARSE INPUT ARGUMENTS %%%
% Fall back to defaults if not provided
if nargin<10
    if ~(isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
        switch which
            case 'mets'
                [~,filen,ext] = fileparts(MRSCont.files{kk});
                figTitle = sprintf(['Load metabolite data plot: ' filen ext '\n']);
            case 'mm'% re_mm
                    [~,filen,ext] = fileparts(MRSCont.files_mm{kk});% re_mm
                    figTitle = sprintf(['Load MM data plot: ' filen ext '\n']);% re_mm
            case 'ref'
                if ~(strcmp(MRSCont.datatype,'P') || strcmp(MRSCont.datatype,'DATA'))
                    [~,filen,ext] = fileparts(MRSCont.files_ref{kk});
                    figTitle = sprintf(['Load water reference data plot: ' filen ext '\n']);
                else
                    [~,filen,ext] = fileparts(MRSCont.files{kk});
                    figTitle = sprintf(['Load interleaved water reference data plot: ' filen ext '\n']);
                end
            case 'w'
                [~,filen,ext] = fileparts(MRSCont.files_w{kk});
                figTitle = sprintf(['Load water data plot: ' filen ext '\n']);
            otherwise
                error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
        end
    else
        if nargin<4
            VoxelIndex = 1; 
        end
        switch which
            case 'mets'
                [~,filen,ext] = fileparts(MRSCont.files{kk});
                figTitle = sprintf(['Load metabolite data plot: ' filen ext '\n Voxel ' num2str(VoxelIndex) ' ']);
            case 'mm'% re_mm
                    [~,filen,ext] = fileparts(MRSCont.files_mm{kk});% re_mm
                    figTitle = sprintf(['Load MM data plot: ' filen ext '\n Voxel ' num2str(VoxelIndex) ' ']);% re_mm
            case 'ref'
                if ~(strcmp(MRSCont.datatype,'P') || strcmp(MRSCont.datatype,'DATA'))
                    [~,filen,ext] = fileparts(MRSCont.files_ref{kk});
                    figTitle = sprintf(['Load water reference data plot: ' filen ext '\n Voxel ' num2str(VoxelIndex) ' ']);
                else
                    [~,filen,ext] = fileparts(MRSCont.files{kk});
                    figTitle = sprintf(['Load interleaved water reference data plot: ' filen ext '\n Voxel ' num2str(VoxelIndex) ' ']);
                end
            case 'w'
                [~,filen,ext] = fileparts(MRSCont.files_w{kk});
                figTitle = sprintf(['Load water data plot: ' filen ext '\n Voxel ' num2str(VoxelIndex) ' ']);
            otherwise
                error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
        end      
    end
    if nargin<9
        ylab='';
        if nargin<8
            xlab='chemical shift (ppm)';
            if nargin<7
                switch which
                    case 'mets'
                        ppmmax = 4.5;
                      case 'mm' %re_mm
                        ppmmax = 4.6;   %re_mm
                    case {'ref', 'w'}
                        ppmmax = 2*4.68;
                    otherwise
                        error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
                end
                if nargin<6
                    switch which
                        case 'mets'
                            ppmmin = 0.2;
                            case 'mm' %re_mm
                        ppmmin = 0.0;   %re_mm
                        case {'ref', 'w'}
                            ppmmin = 0;
                        otherwise
                            error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
                    end
                    if nargin<5
                        stag = 0;
                        if nargin<4
%                             VoxelIndex = [];                        
                            if nargin < 3
                                which = 'mets';
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

%Get y-axis values
maxRef = ones(1,MRSCont.nDatasets);
if isfield(MRSCont,'plot') && (MRSCont.plot.load.match == 1)
    if MRSCont.flags.hasRef
        maxRef = MRSCont.plot.load.ref.max;
    else if MRSCont.flags.hasWater
        maxRef = MRSCont.plot.load.w.max;
        end
    end    
    ymin = min(MRSCont.plot.load.(which).min/maxRef);
    ymax = 1.2*max(MRSCont.plot.load.(which).max/maxRef);
    if ymin < 0
        ymin = 3 * ymin;
    else
        ymin = ymin - 0.5*ymax;
    end
end

%%% 2. EXTRACT DATA TO PLOT %%%
% Extract raw spectra in the plot range
if (isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1)) || (isfield(MRSCont.flags,'isMRSI') && (MRSCont.flags.isMRSI == 1))
        if ~exist('VoxelIndex') && (MRSCont.flags.isPRIAM == 1)
            VoxelIndex = 1;
        elseif ~exist('VoxelIndex') && (MRSCont.flags.isMRSI == 1)
            VoxelIndex = [1 1 1];  
        end
        
    switch which
        case 'mets'
            dataToPlot=op_takeVoxel(MRSCont.raw{kk},VoxelIndex);
            dataToPlot  = op_freqrange(dataToPlot, ppmmin, ppmmax);
        case 'mm' %re_mm
            dataToPlot=op_takeVoxel(MRSCont.raw_mm{kk},VoxelIndex);
            dataToPlot  = op_freqrange(dataToPlot, ppmmin, ppmmax);   %re_mm
        case 'ref'
            dataToPlot=op_takeVoxel(MRSCont.raw_ref{kk},VoxelIndex);
            dataToPlot  = op_freqrange(dataToPlot, ppmmin, ppmmax);
        case 'w'
            dataToPlot=op_takeVoxel(MRSCont.raw_w{kk},VoxelIndex);
            dataToPlot  = op_freqrange(dataToPlot, ppmmin, ppmmax);
        otherwise
            error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
    end
    dataToPlot.specs = dataToPlot.specs/maxRef(kk);
else
    switch which
        case 'mets'
            dataToPlot  = op_freqrange(MRSCont.raw{kk}, ppmmin, ppmmax);
        case 'mm' %re_mm
            dataToPlot  = op_freqrange(MRSCont.raw_mm{kk}, ppmmin, ppmmax);   %re_mm
        case 'ref'
            dataToPlot  = op_freqrange(MRSCont.raw_ref{kk}, ppmmin, ppmmax);
        case 'w'
            dataToPlot  = op_freqrange(MRSCont.raw_w{kk}, ppmmin, ppmmax);
        otherwise
            error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
    end
    dataToPlot.specs = dataToPlot.specs/maxRef(kk);
end


%%% 3. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
out = figure;
% Divide the figure into tiles depending on the number of subspec
% Add the data and plot   
% Staggered plots will be in all black and separated by the mean of the
% maximum across all averages divided by the number of averages
stag = stag*max(abs(mean(max(real(dataToPlot.specs)),2)), abs(mean(min(real(dataToPlot.specs)),2)));

if MRSCont.flags.isUnEdited
    axesHandles.A = gca();
    nAvgs = dataToPlot.averages;
    % Loop over all averages
    for rr = 1:nAvgs
        plot(axesHandles.A, dataToPlot.ppm, dataToPlot.specs(:,rr) + rr*stag, 'k', 'LineWidth', 0.5, 'Color',MRSCont.colormap.Foreground);
        hold on; 
    end
    axesNames = {'A'};
    TitleNames = {'A'};
end
if MRSCont.flags.isMEGA && ~(strcmp(which, 'w') || strcmp(which, 'ref'))
    axesHandles.A  = subplot(2, 1, 1);
    axesHandles.B  = subplot(2, 1, 2);
    nAvgs = dataToPlot.rawAverages/2;
    hold(axesHandles.A, 'on');
    hold(axesHandles.B, 'on');  
    % Loop over all averages
    for rr = 1:nAvgs
        plot(axesHandles.A,dataToPlot.ppm, dataToPlot.specs(:,rr,1) + rr*stag(1), 'k', 'LineWidth', 0.5, 'Color',MRSCont.colormap.Foreground);
        plot(axesHandles.B,dataToPlot.ppm, dataToPlot.specs(:,rr,2) + rr*stag(2), 'k', 'LineWidth', 0.5, 'Color',MRSCont.colormap.Foreground);
    end
    axesNames = {'A','B'};
    TitleNames = {'A','B'};
else if MRSCont.flags.isMEGA && (strcmp(which, 'w') || strcmp(which, 'ref'))
    axesHandles.A = gca();
    nAvgs = dataToPlot.averages;
    % Loop over all averages
    for rr = 1:nAvgs
        plot(axesHandles.A, dataToPlot.ppm, dataToPlot.specs(:,rr) + rr*stag, 'k', 'LineWidth', 0.5, 'Color',MRSCont.colormap.Foreground);
        hold on; 
    end
    axesNames = {'A'};
    if strcmp(which, 'w')
        TitleNames = {'A'};
    else
        TitleNames = {'A & B'};
    end
    end
end
if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
    if ~(strcmp(which, 'w') || strcmp(which, 'ref'))
        figTitle ='';
        axesHandles.A  = subplot(2, 2, 1);
        axesHandles.B  = subplot(2, 2, 2);
        axesHandles.C  = subplot(2, 2, 3);
        axesHandles.D  = subplot(2, 2, 4);
        nAvgs = dataToPlot.averages/4;
        hold(axesHandles.A, 'on');
        hold(axesHandles.B, 'on');  
        hold(axesHandles.C, 'on');
        hold(axesHandles.D, 'on');  
        % Loop over all averages
        for rr = 1:nAvgs
            plot(axesHandles.A,dataToPlot.ppm, real(dataToPlot.specs(:,rr,1) + rr*stag(1)), 'LineWidth', 0.5, 'Color',MRSCont.colormap.Foreground);
            plot(axesHandles.B,dataToPlot.ppm, real(dataToPlot.specs(:,rr,2) + rr*stag(2)), 'LineWidth', 0.5, 'Color',MRSCont.colormap.Foreground);
            plot(axesHandles.C,dataToPlot.ppm, real(dataToPlot.specs(:,rr,3) + rr*stag(3)), 'LineWidth', 0.5, 'Color',MRSCont.colormap.Foreground); 
            plot(axesHandles.D,dataToPlot.ppm, real(dataToPlot.specs(:,rr,4) + rr*stag(4)), 'LineWidth', 0.5, 'Color',MRSCont.colormap.Foreground);
        end
        axesNames = {'A','B','C','D'};
        TitleNames = {'A','B','C','D'};
    else
        axesHandles.A = gca();
        nAvgs = dataToPlot.averages;
        % Loop over all averages
        for rr = 1:nAvgs
            plot(axesHandles.A, dataToPlot.ppm, real(dataToPlot.specs(:,rr) + rr*stag), 'LineWidth', 0.5, 'Color',MRSCont.colormap.Foreground);
            hold on; 
        end
        axesNames = {'A'}; 
        if strcmp(which, 'w')
            TitleNames = {'A'};
        else
            TitleNames = {'A & B & C & D'};
        end
    end
end
hold off;
%%% 4. DESIGN FINETUNING %%%
for ax = 1 : length(axesNames)
    set(axesHandles.(axesNames{ax}), 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'XMinorTick', 'On');
    if isfield(MRSCont,'plot') && (MRSCont.plot.load.match == 1)
        set(axesHandles.(axesNames{ax}), 'YLim', [ymin, ymax]);
    end
    set(axesHandles.(axesNames{ax}), 'LineWidth', 1, 'TickDir', 'out');
    set(axesHandles.(axesNames{ax}), 'FontSize', 16);
    set(axesHandles.(axesNames{ax}), 'Units', 'normalized');
    % If no y caption, remove y axis
    if isempty(ylab)
        set(axesHandles.(axesNames{ax}), 'YColor', MRSCont.colormap.Background);
        set(axesHandles.(axesNames{ax}),'YTickLabel',[])
        set(axesHandles.(axesNames{ax}),'YTick',[])
        set(axesHandles.(axesNames{ax}), 'XColor', MRSCont.colormap.Foreground);
        set(axesHandles.(axesNames{ax}), 'Color', MRSCont.colormap.Background);
        if ax == 1
            title(axesHandles.(axesNames{ax}),[figTitle 'Subspectra ' TitleNames{ax}], 'Interpreter', 'none', 'Color', MRSCont.colormap.Foreground);
        else
            title(axesHandles.(axesNames{ax}),[' Subspectra ' TitleNames{ax}], 'Interpreter', 'none', 'Color', MRSCont.colormap.Foreground);
        end
    else
    set(axesHandles.(axesNames{ax}), 'YColor', 'k');
    end
    % Black axes, white background
    box off;
    xlabel(axesHandles.(axesNames{ax}),xlab, 'FontSize', 16);
    ylabel(axesHandles.(axesNames{ax}),ylab, 'FontSize', 16);
end

set(gcf, 'Color', MRSCont.colormap.Background);
%%% 5. ADD OSPREY LOGO %%%
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end

end