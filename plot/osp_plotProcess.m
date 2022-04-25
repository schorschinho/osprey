function out = osp_plotProcess(MRSCont, kk, which_spec,SubSpectraIndex,ExtraIndex,VoxelIndex, ppmmin, ppmmax)
%% out = osp_plotProcess(MRSCont, kk, which, ppmmin, ppmmax)
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
%       ExtraIndex = Index for the extra dimension
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
if nargin<8
    switch which_spec
        case {'A', 'B', 'C', 'D', 'diff1', 'diff2','diff3', 'sum','A_mm', 'B_mm', 'C_mm', 'D_mm', 'diff1_mm', 'diff2_mm','diff3_mm', 'sum_mm','mm','metab'}
            ppmmax = 4.2;
        case {'ref', 'w','mm_ref'}
            ppmmax = 2*4.68;
        otherwise
            error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
    end
    if nargin<7
        switch which_spec
            case {'A', 'B', 'C', 'D', 'diff1', 'diff2','diff3', 'sum','metab'}
                ppmmin = 0.2;
            case {'ref', 'w','mm','mm_ref','A_mm', 'B_mm', 'C_mm', 'D_mm', 'diff1_mm', 'diff2_mm','diff3_mm', 'sum_mm','mm'}
                ppmmin = 0;
            otherwise
                error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
        end
        if nargin<6
            if MRSCont.flags.isPRIAM
                VoxelIndex = 1; 
            end
            if MRSCont.flags.isMRSI
                VoxelIndex = [1 1]; 
            end
            if nargin < 5
                ExtraIndex = 1;
                if nargin<4
                    SubSpectraIndex = 1;
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

% Set up colormaps
if isfield(MRSCont,'colormap')
    colormap = MRSCont.colormap;
else
    colormap.Background     = [1 1 1];
    colormap.LightAccent    = [110/255 136/255 164/255];
    colormap.Foreground     = [0 0 0];
    colormap.Accent         = [11/255 71/255 111/255];
end

% Pick right entries from the processed struct
switch which_spec
    case {'A', 'B', 'C', 'D','diff1','diff2','diff3','sum','metab'} 
        which_sub_spec = MRSCont.processed.(which_spec){kk}.names{ExtraIndex,SubSpectraIndex};
        which_spec = 'metab';
    case {'A_mm', 'B_mm', 'C_mm', 'D_mm','diff1_mm','diff2_mm','diff3_mm','sum_mm','mm'}
        which_sub_spec = MRSCont.processed.(which_spec){kk}.names{ExtraIndex,SubSpectraIndex};
        which_spec = 'mm';
    otherwise
        which_sub_spec = which_spec;
end  
            
            

%%% 2. EXTRACT DATA TO PLOT %%%
% Extract raw and processed spectra in the plot range
if (isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1)) || (isfield(MRSCont.flags,'isMRSI') && (MRSCont.flags.isMRSI == 1))
        if ~exist('VoxelIndex') && (MRSCont.flags.isPRIAM == 1)
            VoxelIndex = 1;
        elseif ~exist('VoxelIndex') && (MRSCont.flags.isMRSI == 1)
            VoxelIndex = [1 1];  
        end

        switch which_spec
            case {'A', 'B', 'C', 'D'}           
                 raw=op_takeVoxel(MRSCont.raw{kk},VoxelIndex);
                 procDataToPlot=op_takeVoxel(MRSCont.processed.(which_spec){kk},VoxelIndex);

                % Get sub-spectra, depending on whether they are stored as such
                if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                    if raw.subspecs == 4
                        raw_A   = op_takesubspec(raw,procDataToPlot.commuteOrder);                    % Get first subspectrum
                        raw_B   = op_takesubspec(raw,procDataToPlot.commuteOrder);                    % Get second subspectrum
                        raw_C   = op_takesubspec(raw,procDataToPlot.commuteOrder);                    % Get third subspectrum
                        raw_D   = op_takesubspec(raw,procDataToPlot.commuteOrder);                    % Get fourth subspectrum
                    else
                        raw_A   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get first subspectrum
                        raw_B   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get second subspectrum
                        raw_C   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get third subspectrum
                        raw_D   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get fourth subspectrum
                    end
                elseif MRSCont.flags.isMEGA
                    if raw.subspecs == 2
                        raw_A   = op_takesubspec(raw,1);                    % Get first subspectrum
                        raw_B   = op_takesubspec(raw,2);                    % Get second subspectrum
                    else
                        raw_A   = op_takeaverages(raw,1:2:raw.averages);    % Get first subspectrum
                        raw_B   = op_takeaverages(raw,2:2:raw.averages);    % Get second subspectrum
                    end
                    if MRSCont.processed.diff1{kk}.flags.orderswitched
                        temp_spec = raw_A;
                        raw_A = raw_B;
                        raw_B = temp_spec;            
                    end
                elseif MRSCont.flags.isUnEdited
                        raw_A = raw;                                        % Get all averages
                end

                eval(['rawDataToPlot = raw_' which_spec ';']);
                rawDataToScale = raw_A;                                     % This is used to get consistent yLims
            case {'diff1', 'diff2','diff3', 'sum'}
                 rawDataToPlot=op_takeVoxel(MRSCont.raw{kk},VoxelIndex);
                 procDataToPlot=op_takeVoxel(MRSCont.processed.(which_spec){kk},VoxelIndex);
                if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                        temp_spec = cat(3,rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(1)),rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(2)),...
                                        rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(3)),rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(4)));
                        temp_fid = cat(3,rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(1)),rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(2)),...
                                        rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(3)),rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(4)));
                        rawDataToPlot.fids = temp_fid;
                        rawDataToPlot.specs = temp_spec;
                        proc_A   = MRSCont.processed.A{kk};                   % Get first subspectrum
                        proc_B   = MRSCont.processed.B{kk};                  % Get second subspectrum
                        proc_C   = MRSCont.processed.C{kk};                   % Get third subspectrum
                        proc_D   = MRSCont.processed.D{kk};                     % Get fourth subspectrum
                else
                        proc_A   = MRSCont.processed.A{kk};                      % Get first subspectrum
                        proc_B   = MRSCont.processed.B{kk};                      % Get second subspectrum
                        if procDataToPlot.flags.orderswitched
                            temp_spec = rawDataToPlot.specs(:,:,1);
                            rawDataToPlot.specs(:,:,1) = rawDataToPlot.specs(:,:,2);
                            rawDataToPlot.specs(:,:,2) = temp_spec;
                            temp_fids = rawDataToPlot.fids(:,:,1);
                            rawDataToPlot.fids(:,:,1) = rawDataToPlot.fids(:,:,2);
                            rawDataToPlot.fids(:,:,2) = temp_fids;
                        end
                end
                rawDataToScale = rawDataToPlot;                                      % This is used to get consistent yLims
            case 'mm' %re_mm
                rawDataToPlot=op_takeVoxel(MRSCont.raw_mm{kk},VoxelIndex);
                procDataToPlot=op_takeVoxel(MRSCont.processed.(which_spec){kk},VoxelIndex);
                rawDataToScale = rawDataToPlot;               %re_mm                        % This is used to get consistent yLims
            case 'ref'
                rawDataToPlot=op_takeVoxel(MRSCont.raw_ref{kk},VoxelIndex);
                procDataToPlot=op_takeVoxel(MRSCont.processed.(which_spec){kk},VoxelIndex);
                rawDataToScale = rawDataToPlot;                                      % This is used to get consistent yLims
            case 'w'
                rawDataToPlot=op_takeVoxel(MRSCont.raw_w{kk},VoxelIndex);
                procDataToPlot=op_takeVoxel(MRSCont.processed.(which_spec){kk},VoxelIndex);
                rawDataToScale = rawDataToPlot; 
            case 'mm_ref'
                rawDataToPlot=op_takeVoxel(MRSCont.raw_w{kk},VoxelIndex);
                procDataToPlot=op_takeVoxel(MRSCont.processed.(which_spec){kk},VoxelIndex);
                rawDataToScale = rawDataToPlot;    
                 % This is used to get consistent yLims = rawDataToPlot;                                      % This is used to get consistent yLims
            otherwise
                error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
        end

else
    switch which_spec
        case 'metab'
            switch which_sub_spec
                case {'A','B','C','D'} 
                    raw            = MRSCont.raw{ExtraIndex,kk};
                    procDataToPlot = op_takeextra(MRSCont.processed.(which_spec){kk},ExtraIndex);
                    if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES || MRSCont.flags.isMEGA
                        procDataToPlot = op_takesubspec(procDataToPlot,find(strcmp(which_sub_spec,procDataToPlot.names)));
                    end
                    % Get sub-spectra, depending on whether they are stored as such
                    if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                        if raw.subspecs == 4
                            names = {'A','B','C','D'};
                            raw_A   = op_takesubspec(raw,procDataToPlot.commuteOrder(1));                    % Get first subspectrum
                            raw_B   = op_takesubspec(raw,procDataToPlot.commuteOrder(2));                    % Get second subspectrum
                            raw_C   = op_takesubspec(raw,procDataToPlot.commuteOrder(3));                    % Get third subspectrum
                            raw_D   = op_takesubspec(raw,procDataToPlot.commuteOrder(4));                    % Get fourth subspectrum
                        else
                            raw_A   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get first subspectrum
                            raw_B   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get second subspectrum
                            raw_C   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get third subspectrum
                            raw_D   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get fourth subspectrum
                        end
                    elseif MRSCont.flags.isMEGA
                        if raw.subspecs == 2
                            raw_A   = op_takesubspec(raw,1);                    % Get first subspectrum
                            raw_B   = op_takesubspec(raw,2);                    % Get second subspectrum
                        else
                            raw_A   = op_takeaverages(raw,1:2:raw.averages);    % Get first subspectrum
                            raw_B   = op_takeaverages(raw,2:2:raw.averages);    % Get second subspectrum
                        end
                        if MRSCont.processed.metab{kk}.flags.orderswitched
                            temp_spec = raw_A;
                            raw_A = raw_B;
                            raw_B = temp_spec;            
                        end
                    elseif MRSCont.flags.isUnEdited
                            raw_A = raw;                                        % Get all averages
                    end

                    eval(['rawDataToPlot = raw_' which_sub_spec ';']);
                    rawDataToScale = raw_A;                                     % This is used to get consistent yLims
                case {'diff1', 'diff2','diff3', 'sum'}
                    rawDataToPlot  = MRSCont.raw{ExtraIndex,kk};
                    procDataToPlot = op_takeextra(MRSCont.processed.(which_spec){kk},ExtraIndex);
                    procDataToPlot = op_takesubspec(procDataToPlot,find(strcmp(which_sub_spec,procDataToPlot.names)));
                    if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                            temp_spec = cat(3,rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(1)),rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(2)),...
                                            rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(3)),rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(4)));
                            temp_fid = cat(3,rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(1)),rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(2)),...
                                            rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(3)),rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(4)));
                            rawDataToPlot.fids = temp_fid;
                            rawDataToPlot.specs = temp_spec;
                            proc_A   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'A')));                   % Get first subspectrum
                            proc_B   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'B')));                  % Get second subspectrum
                            proc_C   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'C')));                   % Get third subspectrum
                            proc_D   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'D')));                     % Get fourth subspectrum
                    else
                            proc_A   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'A')));                      % Get first subspectrum
                            proc_B   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'B')));                      % Get second subspectrum
                            if procDataToPlot.flags.orderswitched
                                temp_spec = rawDataToPlot.specs(:,:,1);
                                rawDataToPlot.specs(:,:,1) = rawDataToPlot.specs(:,:,2);
                                rawDataToPlot.specs(:,:,2) = temp_spec;
                                temp_fids = rawDataToPlot.fids(:,:,1);
                                rawDataToPlot.fids(:,:,1) = rawDataToPlot.fids(:,:,2);
                                rawDataToPlot.fids(:,:,2) = temp_fids;
                            end
                    end
                    rawDataToScale = rawDataToPlot;                                      % This is used to get consistent yLims
            end
        case 'mm' %re_mm            
            switch which_sub_spec
                case {'A','B','C','D'} 
                    raw  = MRSCont.raw_mm{ExtraIndex,kk}; %re_mm
                    procDataToPlot = op_takeextra(MRSCont.processed.mm{kk},ExtraIndex); %re_mm
                    if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES || MRSCont.flags.isMEGA
                        procDataToPlot = op_takesubspec(procDataToPlot,find(strcmp(which_sub_spec,procDataToPlot.names)));
                    end
                    % Get sub-spectra, depending on whether they are stored as such
                    if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                        if raw.subspecs == 4
                            raw_A   = op_takesubspec(raw,procDataToPlot.commuteOrder);                    % Get first subspectrum
                            raw_B   = op_takesubspec(raw,procDataToPlot.commuteOrder);                    % Get second subspectrum
                            raw_C   = op_takesubspec(raw,procDataToPlot.commuteOrder);                    % Get third subspectrum
                            raw_D   = op_takesubspec(raw,procDataToPlot.commuteOrder);                    % Get fourth subspectrum
                        else
                            raw_A   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get first subspectrum
                            raw_B   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get second subspectrum
                            raw_C   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get third subspectrum
                            raw_D   = op_takeaverages(raw,procDataToPlot.commuteOrder:4:raw.averages);    % Get fourth subspectrum
                        end
                    elseif MRSCont.flags.isMEGA
                        if raw.subspecs == 2
                            raw_A   = op_takesubspec(raw,1);                    % Get first subspectrum
                            raw_B   = op_takesubspec(raw,2);                    % Get second subspectrum
                        else
                            raw_A   = op_takeaverages(raw,1:2:raw.averages);    % Get first subspectrum
                            raw_B   = op_takeaverages(raw,2:2:raw.averages);    % Get second subspectrum
                        end
                        if MRSCont.processed.metab{kk}.flags.orderswitched
                            temp_spec = raw_A;
                            raw_A = raw_B;
                            raw_B = temp_spec;            
                        end
                    elseif MRSCont.flags.isUnEdited
                            raw_A = raw;                                        % Get all averages
                    end

                    eval(['rawDataToPlot = raw_' which_sub_spec ';']);
                    rawDataToScale = raw_A;                                     % This is used to get consistent yLims
                case {'diff1', 'diff2','diff3', 'sum'}
                    rawDataToPlot  = MRSCont.raw_mm{ExtraIndex,kk}; %re_mm
                    procDataToPlot = op_takeextra(MRSCont.processed.mm{kk},ExtraIndex); %re_mm
                    procDataToPlot = op_takesubspec(procDataToPlot,find(strcmp(which_sub_spec,procDataToPlot.names)));
                    if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                            temp_spec = cat(3,rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(1)),rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(2)),...
                                            rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(3)),rawDataToPlot.specs(:,:,procDataToPlot.commuteOrder(4)));
                            temp_fid = cat(3,rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(1)),rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(2)),...
                                            rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(3)),rawDataToPlot.fids(:,:,procDataToPlot.commuteOrder(4)));
                            rawDataToPlot.fids = temp_fid;
                            rawDataToPlot.specs = temp_spec;
                            proc_A   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'A')));                   % Get first subspectrum
                            proc_B   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'B')));                  % Get second subspectrum
                            proc_C   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'C')));                   % Get third subspectrum
                            proc_D   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'D')));                     % Get fourth subspectrum
                    else
                            proc_A   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'A')));                      % Get first subspectrum
                            proc_B   = op_takesubspec(MRSCont.processed.(which_spec){ExtraIndex,kk},find(strcmp(which_sub_spec,'B')));                      % Get second subspectrum
                            if procDataToPlot.flags.orderswitched
                                temp_spec = rawDataToPlot.specs(:,:,1);
                                rawDataToPlot.specs(:,:,1) = rawDataToPlot.specs(:,:,2);
                                rawDataToPlot.specs(:,:,2) = temp_spec;
                                temp_fids = rawDataToPlot.fids(:,:,1);
                                rawDataToPlot.fids(:,:,1) = rawDataToPlot.fids(:,:,2);
                                rawDataToPlot.fids(:,:,2) = temp_fids;
                            end
                    end
                    rawDataToScale = rawDataToPlot;                                      % This is used to get consistent yLims
            end            
        case 'ref'
            rawDataToPlot  = MRSCont.raw_ref{ExtraIndex,kk};
            procDataToPlot = op_takeextra(MRSCont.processed.ref{kk},ExtraIndex);
            rawDataToScale = rawDataToPlot;                                      % This is used to get consistent yLims
        case 'w'
            rawDataToPlot  = MRSCont.raw_w{ExtraIndex,kk};
            procDataToPlot = op_takeextra(MRSCont.processed.w{kk},ExtraIndex);
            rawDataToScale = rawDataToPlot; 
       case 'mm_ref'
            rawDataToPlot  = MRSCont.raw_mm_ref{ExtraIndex,kk};
            procDataToPlot = op_takeextra(MRSCont.processed.mm_ref{kk},ExtraIndex);
            rawDataToScale = rawDataToPlot;      
             % This is used to get consistent yLims = rawDataToPlot;                                      % This is used to get consistent yLims
        otherwise
            error('Input for variable ''which'' not recognized. Needs to be ''mets'' (metabolite data), ''ref'' (reference data), or ''w'' (short-TE water data).');
    end
end
%%% 3. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
if ~MRSCont.flags.isGUI
    out = figure;
else
    out = figure('Visible','off');
end

% Divide the figure into six tiles, create four axes
ax_raw      = subplot(2, 2, 1);
ax_aligned  = subplot(2, 2, 3);
ax_proc     = subplot(2, 2, 4);
ax_drift    = subplot(2, 2, 2);


%%% 4. PLOT RAW UNALIGNED %%%
% Generate global yLimits
applyDataToScale = rawDataToScale;
t = rawDataToScale.t;
switch which_sub_spec
    case {'A', 'B', 'C', 'D'} 
        if (isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1)) 
            fs = procDataToPlot.specReg{VoxelIndex}.fs;
            phs = procDataToPlot.specReg{VoxelIndex}.phs;
        elseif(isfield(MRSCont.flags,'isMRSI') && (MRSCont.flags.isMRSI == 1))
            fs = procDataToPlot.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
            phs = procDataToPlot.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
        else
            fs = procDataToPlot.specReg{1}.fs(:,SubSpectraIndex);
            phs = procDataToPlot.specReg{1}.phs(:,SubSpectraIndex);           
        end
    case {'diff1', 'diff2','diff3', 'sum'}
        if (isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1)) 
            if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                
                fs{1} = proc_A.specReg{VoxelIndex}.fs;
                phs{1} = proc_A.specReg{VoxelIndex}.phs;
                fs{2} = proc_B.specReg{VoxelIndex}.fs;
                phs{2} = proc_B.specReg{VoxelIndex}.phs;
                fs{3} = proc_C.specReg{VoxelIndex}.fs;
                phs{3} = proc_C.specReg{VoxelIndex}.phs;
                fs{4} = proc_D.specReg{VoxelIndex}.fs;
                phs{4} = proc_D.specReg{VoxelIndex}.phs;
            else
                    fs{1} = proc_A.specReg{VoxelIndex}.fs;
                    phs{1} = proc_A.specReg{VoxelIndex}.phs;
                    fs{2} = proc_B.specReg{VoxelIndex}.fs;
                    phs{2} = proc_B.specReg{VoxelIndex}.phs;
            end
        elseif(isfield(MRSCont.flags,'isMRSI') && (MRSCont.flags.isMRSI == 1))
            if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                
                fs{1} = proc_A.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                phs{1} = proc_A.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
                fs{2} = proc_B.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                phs{2} = proc_B.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
                fs{3} = proc_C.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                phs{3} = proc_C.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
                fs{4} = proc_D.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                phs{4} = proc_D.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
            else
                    fs{1} = proc_A.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                    phs{1} = proc_A.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
                    fs{2} = proc_B.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                    phs{2} = proc_B.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
            end
        else
            if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                fs{1} = proc_A.specReg{1}.fs;
                phs{1} = proc_A.specReg{1}.phs;
                fs{2} = proc_B.specReg{1}.fs;
                phs{2} = proc_B.specReg{1}.phs;
                fs{3} = proc_C.specReg{1}.fs;
                phs{3} = proc_C.specReg{1}.phs;
                fs{4} = proc_D.specReg{1}.fs;
                phs{4} = proc_D.specReg{1}.phs;
            else
                fs{1} = proc_A.specReg{1}.fs;
                phs{1} = proc_A.specReg{1}.phs;
                fs{2} = proc_B.specReg{1}.fs;
                phs{2} = proc_B.specReg{1}.phs;
              end
        end
end

if (isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1)) 
    if isfield(MRSCont.QM{VoxelIndex}.freqShift, which_sub_spec)
        switch which_sub_spec
            case {'A', 'B', 'C', 'D'} 
                refShift = -repmat(MRSCont.QM{VoxelIndex}.freqShift.(which_sub_spec)(kk), size(fs));
                fs = fs - refShift;
                for jj = 1:size(applyDataToScale.fids,2)
                    applyDataToScale.fids(:,jj) = applyDataToScale.fids(:,jj) .* ...
                        exp(1i*fs(jj)*2*pi*t') * exp(1i*pi/180*phs(jj));
                end
            case {'diff1', 'diff2','diff3', 'sum'}
                refShift = -repmat(MRSCont.QM{VoxelIndex}.freqShift.(which_sub_spec)(kk), size(fs{1}));
                for ss = 1 : length(fs)
                    fs{ss} = fs{ss} - refShift;
                    for jj = 1:size(applyDataToScale.fids,2)
                        applyDataToScale.fids(:,jj,ss) = applyDataToScale.fids(:,jj,ss) .* ...
                            exp(1i*fs{ss}(jj)*2*pi*t') * exp(1i*pi/180*phs{ss}(jj));
                    end
                end
        end
    end
elseif (isfield(MRSCont.flags,'isMRSI') && (MRSCont.flags.isMRSI == 1))
    if isfield(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift, which_sub_spec)
        switch which_sub_spec
            case {'A', 'B', 'C', 'D'} 
                refShift = -repmat(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift.(which_sub_spec)(kk), size(fs));
                fs = fs - refShift;
                for jj = 1:size(applyDataToScale.fids,2)
                    applyDataToScale.fids(:,jj) = applyDataToScale.fids(:,jj) .* ...
                        exp(1i*fs(jj)*2*pi*t') * exp(1i*pi/180*phs(jj));
                end
            case {'diff1', 'diff2','diff3', 'sum'}
                refShift = -repmat(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift.(which_sub_spec)(kk), size(fs{1}));
                for ss = 1 : length(fs)
                    fs{ss} = fs{ss} - refShift;
                    for jj = 1:size(applyDataToScale.fids,2)
                        applyDataToScale.fids(:,jj,ss) = applyDataToScale.fids(:,jj,ss) .* ...
                            exp(1i*fs{ss}(jj)*2*pi*t') * exp(1i*pi/180*phs{ss}(jj));
                    end
                end
        end
    end
else
    if isfield(MRSCont.QM.freqShift, which_sub_spec)
        switch which_sub_spec
            case {'A', 'B', 'C', 'D'} 
                refShift = -repmat(MRSCont.QM.freqShift.(which_sub_spec)(kk), size(fs));
                fs = fs - refShift;
                for jj = 1:size(applyDataToScale.fids,2)
                    applyDataToScale.fids(:,jj) = applyDataToScale.fids(:,jj) .* ...
                        exp(1i*fs(jj)*2*pi*t') * exp(1i*pi/180*phs(jj));
                end
            case {'diff1', 'diff2','diff3', 'sum'}
                refShift = -repmat(MRSCont.QM.freqShift.(which_sub_spec)(kk), size(fs{1}));
                for ss = 1 : length(fs)
                    fs{ss} = fs{ss} - refShift;
                    for jj = 1:size(applyDataToScale.fids,2)
                        if size(applyDataToScale.fids) == 3
                            applyDataToScale.fids(:,jj,ss) = applyDataToScale.fids(:,jj,ss) .* ...
                                exp(1i*fs{ss}(jj)*2*pi*t') * exp(1i*pi/180*phs{ss}(jj));
                        else
                            applyDataToScale.fids(:,ss) = applyDataToScale.fids(:,ss) .* ...
                                exp(1i*fs{ss}(1)*2*pi*t') * exp(1i*pi/180*phs{ss}(1));
                        end
                    end
                end
        end
    end    
end


applyDataToScale.specs = fftshift(fft(applyDataToScale.fids,[],rawDataToScale.dims.t),rawDataToScale.dims.t);


plotRangeScale = op_freqrange(applyDataToScale, ppmmin, ppmmax);
yLims= [ min(min(real(plotRangeScale.specs(:,:)))) max(max(real(plotRangeScale.specs(:,:))))];
yLimsAbs = (abs(yLims(1)) +  abs(yLims(2)));
if strcmp(which_sub_spec, 'diff1') || strcmp(which_sub_spec, 'diff2') || strcmp(which_sub_spec, 'diff3') || strcmp(which_sub_spec, 'sum')
    if MRSCont.flags.isMEGA
        yLims = [yLims(1) - (yLimsAbs*0.1) (2*yLims(2)) + (yLimsAbs*0.1)];
    else
        yLims = [yLims(1) - (yLimsAbs*0.1) (3*yLims(2)) + (yLimsAbs*0.1)];
    end
else
    yLims = [yLims(1) - (yLimsAbs*0.1) yLims(2) + (yLimsAbs*0.1)];
end

% Add the data and plot
hold(ax_raw, 'on');    
% Loop over all averages
try
    nAvgsRaw = rawDataToPlot.sz(rawDataToPlot.dims.averages);
catch % This is a wild guess in case no averages dimension is stored 
    nAvgsRaw = rawDataToPlot.sz(2);
end
if MRSCont.flags.isUnEdited
    for rr = 1:nAvgsRaw
        plot(ax_raw, rawDataToPlot.ppm, real(rawDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
    end
    set(ax_raw, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
    y = yLims;
end

if MRSCont.flags.isMEGA 
  if ~strcmp(which_sub_spec, 'w') && ~strcmp(which_sub_spec, 'ref')&& ~strcmp(which_sub_spec, 'mm_ref') && ~strcmp(which_sub_spec, 'A') && ~strcmp(which_sub_spec, 'B')
    stag = [0,0.5] .* yLimsAbs;
    stagText = stag + (0.25.* yLimsAbs);
    for rr = 1:nAvgsRaw
        plot(ax_raw, rawDataToPlot.ppm, real(rawDataToPlot.specs(:,rr,1)), 'LineWidth', 0.5, 'Color', colormap.LightAccent);
        plot(ax_raw, rawDataToPlot.ppm, real(rawDataToPlot.specs(:,rr,2) + stag(2)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
    end

    plotRange = op_freqrange(rawDataToPlot, ppmmin, ppmmax);
    if size(rawDataToPlot.fids) == 3
        yLims = [mean(min(real(plotRange.specs(:,:,1)))) (mean(max(real(plotRange.specs(:,:,2))))+stag(2))].*1.5;
    else
        yLims = [mean(min(real(plotRange.specs(:,1)))) (mean(max(real(plotRange.specs(:,2))))+stag(2))].*1.5;
    end
    text(ax_raw, ppmmin+0.3, stagText(1), 'off', 'Color', colormap.LightAccent);
    text(ax_raw, ppmmin+0.3, stagText(2) , 'on', 'Color', colormap.Foreground); 
    set(ax_raw, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);  
    y = yLims;
  else
        for rr = 1:nAvgsRaw
            plot(ax_raw, rawDataToPlot.ppm, real(rawDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
        end
        set(ax_raw, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
        y = yLims;
  end
end

if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
    if ~strcmp(which_sub_spec, 'w') && ~strcmp(which_sub_spec, 'ref')&& ~strcmp(which_sub_spec, 'mm_ref') && ~strcmp(which_sub_spec, 'A') && ~strcmp(which_sub_spec, 'B') && ~strcmp(which_sub_spec, 'C') && ~strcmp(which_sub_spec, 'D')
        stag = [0,0.5,1,1.5] .* yLimsAbs;
        stagText = stag + (0.25.* yLimsAbs);
        for rr = 1:nAvgsRaw
            plot(ax_raw, rawDataToPlot.ppm, real(rawDataToPlot.specs(:,rr,1)), 'LineWidth', 0.5, 'Color', colormap.LightAccent);
            plot(ax_raw, rawDataToPlot.ppm, real(rawDataToPlot.specs(:,rr,2) + stag(2)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
            plot(ax_raw, rawDataToPlot.ppm, real(rawDataToPlot.specs(:,rr,3) + stag(3)), 'LineWidth', 0.5, 'Color', colormap.LightAccent);
            plot(ax_raw, rawDataToPlot.ppm, real(rawDataToPlot.specs(:,rr,4) + stag(4)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
        end
        text(ax_raw, ppmmin+0.3, stagText(1), 'A', 'Color', colormap.LightAccent);
        text(ax_raw, ppmmin+0.3, stagText(2), 'B', 'Color', colormap.Foreground);
        text(ax_raw, ppmmin+0.3, stagText(3), 'C', 'Color', colormap.LightAccent);
        text(ax_raw, ppmmin+0.3, stagText(4), 'D', 'Color', colormap.Foreground);  
        set(ax_raw, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
        y = yLims;       
    else
        for rr = 1:nAvgsRaw
            plot(ax_raw, rawDataToPlot.ppm, real(rawDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
        end
        set(ax_raw, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
        y = yLims;
    end    
end

if ~(strcmp(which_sub_spec,'w') || strcmp(which_sub_spec,'ref')|| strcmp(which_sub_spec,'mm_ref'))
    plot(ax_raw, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground, 'LineWidth', 0.5);
    plot(ax_raw, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground, 'LineWidth', 0.5);
    if ~strcmp(which_sub_spec, 'mm')
        plot(ax_raw, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5); 
    else
        plot(ax_raw, [3.9 3.9], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5);
    end        
end
hold(ax_raw, 'off');
title(ax_raw, 'Pre-alignment', 'Color', colormap.Foreground);
xlabel(ax_raw, 'chemical shift (ppm)', 'Color', colormap.Foreground)
if MRSCont.flags.isGUI
    set(ax_raw, 'YColor', colormap.Background);
    set(ax_raw,'YTickLabel',{})
    set(ax_raw,'YTick',{})
end


%%% 5. PLOT RAW ALIGNED %%%
% Apply stored corrections to calculate the spectra to display
applyDataToPlot = rawDataToPlot;
t = rawDataToPlot.t;
switch which_sub_spec
    case {'A', 'B', 'C', 'D'} 
        if (isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1)) 
            fs = procDataToPlot.specReg{VoxelIndex}.fs;
            phs = procDataToPlot.specReg{VoxelIndex}.phs;
            weights = MRSCont.processed.(which_spec){kk}.specReg{VoxelIndex(1)}.weights;
        elseif (isfield(MRSCont.flags,'isMRSI') && (MRSCont.flags.isMRSI == 1))
            fs = procDataToPlot.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
            phs = procDataToPlot.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
            weights = MRSCont.processed.(which_spec){kk}.specReg{VoxelIndex(1), VoxelIndex(2)}.weights;
        else
            fs = procDataToPlot.specReg{1}.fs;
            phs = procDataToPlot.specReg{1}.phs;
            weights = MRSCont.processed.(which_spec){kk}.specReg{1}.weights{find(strcmp(which_sub_spec,{'A', 'B', 'C', 'D'}))};            
        end
    case {'diff1', 'diff2', 'diff3', 'sum'}
        if (isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1))
            if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                    fs{1} = proc_A.specReg{VoxelIndex}.fs;
                    phs{1} = proc_A.specReg{VoxelIndex}.phs;
                    fs{2} = proc_B.specReg{VoxelIndex}.fs;
                    phs{2} = proc_B.specReg{VoxelIndex}.phs;
                    fs{3} = proc_C.specReg{VoxelIndex}.fs;
                    phs{3} = proc_C.specReg{VoxelIndex}.phs;
                    fs{4} = proc_D.specReg{VoxelIndex}.fs;
                    phs{4} = proc_D.specReg{VoxelIndex}.phs;
                    weights = vertcat(MRSCont.processed.(which_spec){kk}.specReg{VoxelIndex}.weights{:});
            else
                    fs{1} = proc_A.specReg{VoxelIndex}.fs;
                    phs{1} = proc_A.specReg{VoxelIndex}.phs;
                    fs{2} = proc_B.specReg{VoxelIndex}.fs;
                    phs{2} = proc_B.specReg{VoxelIndex}.phs;
                    weights = vertcat(MRSCont.processed.(which_spec){kk}.specReg{VoxelIndex}.weights{:});
            end
        elseif(isfield(MRSCont.flags,'isMRSI') && (MRSCont.flags.isMRSI == 1))
            if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                    fs{1} = proc_A.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                    phs{1} = proc_A.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
                    fs{2} = proc_B.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                    phs{2} = proc_B.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
                    fs{3} = proc_C.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                    phs{3} = proc_C.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
                    fs{4} = proc_D.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                    phs{4} = proc_D.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
                    weights = vertcat(RSCont.processed.(which_spec){kk}.specReg{VoxelIndex(1), VoxelIndex(2)}.weights{:});
            else
                    fs{1} = proc_A.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                    phs{1} = proc_A.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
                    fs{2} = proc_B.specReg{VoxelIndex(1), VoxelIndex(2)}.fs;
                    phs{2} = proc_B.specReg{VoxelIndex(1), VoxelIndex(2)}.phs;
                    weights = vertcat(MRSCont.processed.(which_spec){kk}.specReg{VoxelIndex(1), VoxelIndex(2)}.weights{:});
            end
        else
            if MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                    fs{1} = proc_A.specReg{1}.fs;
                    phs{1} = proc_A.specReg{1}.phs;
                    fs{2} = proc_B.specReg{1}.fs;
                    phs{2} = proc_B.specReg{1}.phs;
                    fs{3} = proc_C.specReg{1}.fs;
                    phs{3} = proc_C.specReg{1}.phs;
                    fs{4} = proc_D.specReg{1}.fs;
                    phs{4} = proc_D.specReg{1}.phs;
                    weights = vertcat(MRSCont.processed.(which_spec){kk}.specReg{1}.weights{:});
            else
                    fs{1} = proc_A.specReg{1}.fs;
                    phs{1} = proc_A.specReg{1}.phs;
                    fs{2} = proc_B.specReg{1}.fs;
                    phs{2} = proc_B.specReg{1}.phs;
                    weights = vertcat(MRSCont.processed.(which_spec){kk}.specReg{1}.weights{:});
            end
        end
end

if (isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1)) 
    if isfield(MRSCont.QM{VoxelIndex(1)}.freqShift, which_sub_spec)
        switch which_sub_spec
            case {'A', 'B', 'C', 'D'} 
                refShift = -repmat(MRSCont.QM{VoxelIndex(1)}.freqShift.(which_spec)(kk), size(fs));
                fs = fs - refShift;
                for jj = 1:size(applyDataToPlot.fids,2)
                    applyDataToPlot.fids(:,jj) = applyDataToPlot.fids(:,jj) .* ...
                        exp(1i*fs(jj)*2*pi*t') * exp(1i*pi/180*phs(jj));
                end
            case {'diff1', 'diff2', 'sum'}
                refShift = -repmat(MRSCont.QM{VoxelIndex(1)}.freqShift.(which_spec)(kk), size(fs{1}));
                for ss = 1 : length(fs)
                    fs{ss} = fs{ss} - refShift;
                    for jj = 1:size(applyDataToPlot.fids,2)
                        applyDataToPlot.fids(:,jj,ss) = applyDataToPlot.fids(:,jj,ss) .* ...
                            exp(1i*fs{ss}(jj)*2*pi*t') * exp(1i*pi/180*phs{ss}(jj));
                    end
                end
        end
    end
elseif (isfield(MRSCont.flags,'isMRSI') && (MRSCont.flags.isMRSI == 1))
    if isfield(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift, which_sub_spec)
        switch which_sub_spec
            case {'A', 'B', 'C', 'D'} 
                refShift = -repmat(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift.(which_spec)(kk), size(fs));
                fs = fs - refShift;
                for jj = 1:size(applyDataToPlot.fids,2)
                    applyDataToPlot.fids(:,jj) = applyDataToPlot.fids(:,jj) .* ...
                        exp(1i*fs(jj)*2*pi*t') * exp(1i*pi/180*phs(jj));
                end
            case {'diff1', 'diff2', 'diff3', 'sum'}
                refShift = -repmat(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift.(which_spec)(kk), size(fs{1}));
                for ss = 1 : length(fs)
                    fs{ss} = fs{ss} - refShift;
                    for jj = 1:size(applyDataToPlot.fids,2)
                        applyDataToPlot.fids(:,jj,ss) = applyDataToPlot.fids(:,jj,ss) .* ...
                            exp(1i*fs{ss}(jj)*2*pi*t') * exp(1i*pi/180*phs{ss}(jj));
                    end
                end
        end
    end
else
    if isfield(MRSCont.QM.freqShift, which_sub_spec)
        switch which_sub_spec
            case {'A', 'B', 'C', 'D'} 
                refShift = -repmat(MRSCont.QM.freqShift.(which_sub_spec)(kk), size(fs));
                fs = fs - refShift;
                for jj = 1:size(applyDataToPlot.fids,2)
                    applyDataToPlot.fids(:,jj) = applyDataToPlot.fids(:,jj) .* ...
                        exp(1i*fs(jj)*2*pi*t') * exp(1i*pi/180*phs(jj));
                end
            case {'diff1', 'diff2', 'diff3', 'sum'}
                refShift = -repmat(MRSCont.QM.freqShift.(which_sub_spec)(kk), size(fs{1}));
                for ss = 1 : length(fs)
                    fs{ss} = fs{ss} - refShift;
                    for jj = 1:size(applyDataToPlot.fids,2)
                        if size(applyDataToPlot.fids) == 3
                            applyDataToPlot.fids(:,jj,ss) = applyDataToPlot.fids(:,jj,ss) .* ...
                                exp(1i*fs{ss}(jj)*2*pi*t') * exp(1i*pi/180*phs{ss}(jj));
                        else
                            applyDataToPlot.fids(:,ss) = applyDataToPlot.fids(:,ss) .* ...
                                exp(1i*fs{ss}(1)*2*pi*t') * exp(1i*pi/180*phs{ss}(1));
                        end
                    end
                end
        end
    end
end
applyDataToPlot.specs = fftshift(fft(applyDataToPlot.fids,[],rawDataToPlot.dims.t),rawDataToPlot.dims.t);

hold(ax_aligned, 'on');    
% Loop over all averages
if MRSCont.flags.isUnEdited
    for rr = 1:nAvgsRaw
        plot(ax_aligned, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
    end
    set(ax_aligned, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
end

if MRSCont.flags.isMEGA
    if ~strcmp(which_sub_spec, 'w') && ~strcmp(which_sub_spec, 'ref') && ~strcmp(which_sub_spec, 'mm_ref') && ~strcmp(which_sub_spec, 'A') && ~strcmp(which_sub_spec, 'B')
        stag = [0,0.5,1,1.5] .* yLimsAbs;
        stagText = stag + (0.25.* yLimsAbs);
        for rr = 1:nAvgsRaw
            plot(ax_aligned, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,1)), 'LineWidth', 0.5, 'Color', colormap.LightAccent);
            plot(ax_aligned, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,2) + stag(2)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
        end
        text(ax_aligned, ppmmin+0.3, stagText(1), 'off', 'Color', colormap.LightAccent);
        text(ax_aligned, ppmmin+0.3, stagText(2), 'on', 'Color', colormap.Foreground);
        set(ax_aligned, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
    else
        for rr = 1:nAvgsRaw
            plot(ax_aligned, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
        end
        set(ax_aligned, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
    end
end


if (MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES)
    if ~strcmp(which_sub_spec, 'w') && ~strcmp(which_sub_spec, 'ref') && ~strcmp(which_sub_spec, 'A') && ~strcmp(which_sub_spec, 'B') && ~strcmp(which_sub_spec, 'C') && ~strcmp(which_sub_spec, 'D')
        stag = [0,0.5,1,1.5] .* yLimsAbs;
        stagText = stag + (0.25.* yLimsAbs);
        for rr = 1:nAvgsRaw
            plot(ax_aligned, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,1)), 'LineWidth', 0.5, 'Color', colormap.LightAccent);
            plot(ax_aligned, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,2) + stag(2)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
            plot(ax_aligned, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,3) + stag(3)), 'LineWidth', 0.5, 'Color', colormap.LightAccent);
            plot(ax_aligned, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr,4) + stag(4)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
        end
        text(ax_aligned, ppmmin+0.3, stagText(1), 'A', 'Color', colormap.LightAccent);
        text(ax_aligned, ppmmin+0.3, stagText(2), 'B', 'Color', colormap.Foreground);
        text(ax_aligned, ppmmin+0.3, stagText(3), 'C', 'Color', colormap.LightAccent);
        text(ax_aligned, ppmmin+0.3, stagText(4), 'D', 'Color', colormap.Foreground); 
        set(ax_aligned, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);           
    else
        for rr = 1:nAvgsRaw
            plot(ax_aligned, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormap.Foreground);
            if ~strcmp(which_sub_spec, 'w') && ~strcmp(which_sub_spec, 'ref')
                plot(ax_aligned, applyDataToPlot.ppm, real(applyDataToPlot.specs(:,rr)), 'LineWidth', 0.5, 'Color', colormap.Foreground);           
            end
        end
        set(ax_aligned, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', yLims);
    end
end

if ~(strcmp(which_sub_spec,'w') || strcmp(which_sub_spec,'ref')|| strcmp(which_sub_spec,'mm_ref'))
    plot(ax_aligned, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground, 'LineWidth', 0.5);
    plot(ax_aligned, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground, 'LineWidth', 0.5);
    if ~strcmp(which_sub_spec, 'mm')
        plot(ax_aligned, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5); 
    else
        plot(ax_aligned, [3.9 3.9], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5);
    end    
end
hold(ax_aligned, 'off');
title(ax_aligned, 'Post-alignment', 'Color', colormap.Foreground);
xlabel(ax_aligned, 'chemical shift (ppm)', 'Color', colormap.Foreground)
if MRSCont.flags.isGUI
    set(ax_aligned, 'YColor', colormap.Background);
    set(ax_aligned,'YTickLabel',{})
    set(ax_aligned,'YTick',{})
end


%%% 6. PLOT PROCESSED %%%
% Add the data and plot
hold(ax_proc, 'on');
if isfield(MRSCont,'plot') && isfield(MRSCont.plot, 'processed') && (MRSCont.plot.processed.match == 1)
    if MRSCont.flags.hasRef || MRSCont.flags.hasWater
        if MRSCont.flags.hasRef
            plot(ax_proc, procDataToPlot.ppm, real(procDataToPlot.specs)/MRSCont.plot.processed.ref.max(kk), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1.5);
            maxRef = MRSCont.plot.processed.ref.max;
        else
            plot(ax_proc, procDataToPlot.ppm, real(procDataToPlot.specs)/MRSCont.plot.processed.w.max(kk), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1.5);
            maxRef = MRSCont.plot.processed.w.max;
        end
        if strcmp(which_sub_spec, 'diff1') || strcmp(which_sub_spec, 'diff2') || strcmp(which_sub_spec, 'sum')
            if strcmp(which_sub_spec, 'sum')
                y = [-0.2*max(MRSCont.plot.processed.(which_spec).max/maxRef), 1.2*max(MRSCont.plot.processed.(which_spec).max/maxRef)];
            else
                if strcmp(which_sub_spec, 'diff1')
                    switch MRSCont.opts.editTarget{1}
                        case 'GABA'
                            y = [-0.2*max(MRSCont.plot.processed.(which_spec).max/maxRef), 1.2*max(MRSCont.plot.processed.(which_spec).max/maxRef)];
                        case 'GSH'
                            y = [1.1*min(MRSCont.plot.processed.(which_spec).min/maxRef), 1.2*max(MRSCont.plot.processed.(which_spec).max/maxRef)];
                        case 'Lac'
                            y = [1.1*min(MRSCont.plot.processed.(which_spec).min/maxRef), 1.2*max(MRSCont.plot.processed.(which_spec).max/maxRef)];
                        otherwise
                            y = [1.1*min(MRSCont.plot.processed.(which_spec).min/maxRef), 1.2*max(MRSCont.plot.processed.(which_spec).max/maxRef)];
                    end
                else
                    switch MRSCont.opts.editTarget{2}
                        case 'GABA'
                            y = [-0.2*max(MRSCont.plot.processed.(which_spec).max/maxRef), 1.2*max(MRSCont.plot.processed.(which_spec).max/maxRef)];
                        case 'GSH'
                            y = [1.1*min(MRSCont.plot.processed.(which_spec).min/maxRef), 1.2*max(MRSCont.plot.processed.(which_spec).max/maxRef)];
                        case 'Lac'
                            y = [1.1*min(MRSCont.plot.processed.(which_spec).min/maxRef), 1.2*max(MRSCont.plot.processed.(which_spec).max/maxRef)];
                        otherwise
                            y = [1.1*min(MRSCont.plot.processed.(which_spec).min/maxRef), 1.2*max(MRSCont.plot.processed.(which_spec).max/maxRef)];
                    end
                end
                            
            end
        else
            y = [-0.2*max(MRSCont.plot.processed.(which_spec).max/maxRef), 1.2*max(MRSCont.plot.processed.(which_spec).max/maxRef)]; 
        end
    else        
        plot(ax_proc, procDataToPlot.ppm, real(procDataToPlot.specs(:,1)), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1.5);
        if strcmp(which_sub_spec, 'diff1') || strcmp(which_sub_spec, 'diff2') || strcmp(which_sub_spec, 'diff3') || strcmp(which_sub_spec, 'sum')
            y = [min(MRSCont.plot.processed.(which_spec).min) max(MRSCont.plot.processed.(which_spec).max)];
        else
            if MRSCont.flags.isMEGA || MRSCont.flags.isHERMES || MRSCont.flags.isHERCULES
                y = [min(MRSCont.plot.processed.(which_spec).min) max(MRSCont.plot.processed.(which_spec).max)];
            else
                y = [min(MRSCont.plot.processed.(which_spec).min) max(MRSCont.plot.processed.(which_spec).max)];
            end        
        end        
    end
else
    plot(ax_proc, procDataToPlot.ppm, real(procDataToPlot.specs(:,1))/max(real(procDataToPlot.specs(procDataToPlot.ppm>ppmmin&procDataToPlot.ppm<ppmmax))), 'Color',MRSCont.colormap.Foreground, 'LineWidth', 1.5);
    if strcmp(which_sub_spec,'diff2')
        y = [-1.2, 1.2];
    else if strcmp(which_sub_spec,'diff1')
            if strcmp(which_spec,'metab')
                y = [-2, 1.2];   
            else
                y = [-1.5, 1.2];   
            end
        else
            if strcmp(which_spec,'metab')
                y = [-0.2, 1.2];
            else
                y = [-1.5, 1.2];   
            end

        end
    end
end
set(ax_proc, 'XDir', 'reverse', 'XLim', [ppmmin, ppmmax], 'YLim', y);
if ~(strcmp(which_sub_spec,'w') || strcmp(which_sub_spec,'ref')|| strcmp(which_sub_spec,'MM_ref'))
    plot(ax_proc, [2.008 2.008], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5);
    plot(ax_proc, [3.027 3.027], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5);
    if ~strcmp(which_sub_spec, 'mm')
        plot(ax_proc, [3.200 3.200], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5); 
    else
        plot(ax_proc, [3.9 3.9], [y(1)-y(2) y(2)],'LineStyle', ':', 'Color', colormap.Foreground,  'LineWidth', 0.5);
    end
end
hold(ax_proc, 'off');
if MRSCont.plot.processed.match == 0
    title(ax_proc, 'Aligned and averaged', 'Color', colormap.Foreground);
else
    title(ax_proc, 'Aligned, averaged, and scaled', 'Color', colormap.Foreground);
end
xlabel(ax_proc, 'chemical shift (ppm)', 'Color', colormap.Foreground)
if MRSCont.flags.isGUI
    set(ax_proc, 'YColor', colormap.Background);
    set(ax_proc,'YTickLabel',{})
    set(ax_proc,'YTick',{})
end

%%% 7. GENERATE DRIFT PLOT %%%
if (isfield(MRSCont.flags,'isPRIAM') && (MRSCont.flags.isPRIAM == 1)) 
    if isfield(MRSCont.QM{VoxelIndex}.drift.pre, which_sub_spec)
        if length(MRSCont.QM{VoxelIndex}.drift.pre.(which_sub_spec){kk}) > 1
            crDriftPre = MRSCont.QM{VoxelIndex}.drift.pre.(which_sub_spec){kk} + MRSCont.QM{VoxelIndex}.freqShift.(which_sub_spec)(kk)/applyDataToPlot.txfrq*1e6;
            crDriftPost = MRSCont.QM{VoxelIndex}.drift.post.(which_sub_spec){kk} + MRSCont.QM{VoxelIndex}.freqShift.(which_sub_spec)(kk)/applyDataToPlot.txfrq*1e6;
            hold(ax_drift, 'on');
            colors = ones(length(crDriftPre),1).*colormap.Foreground;
            for dots = 1 : length(crDriftPre)
                colors(dots,1) = colors(dots,1) + (1 - colors(dots,1)) * (1-weights(dots));
                colors(dots,2) = colors(dots,2) + (1 - colors(dots,2)) * (1-weights(dots));
                colors(dots,3) = colors(dots,3) + (1 - colors(dots,3)) * (1-weights(dots));
            end
            scatter(ax_drift, [1:length(crDriftPre)],crDriftPre'-(MRSCont.QM{VoxelIndex}.freqShift.(which_sub_spec)(kk)/procDataToPlot.txfrq*1e6),36,ones(length(crDriftPre),1).*colormap.LightAccent);
            scatter(ax_drift, [1:length(crDriftPost)],crDriftPost'-(MRSCont.QM{VoxelIndex}.freqShift.(which_sub_spec)(kk)/procDataToPlot.txfrq*1e6),36,colors,'filled','MarkerEdgeColor',colormap.Foreground);    

            text(ax_drift, length(crDriftPre)*1.05, crDriftPre(end)-(MRSCont.QM{VoxelIndex}.freqShift.(which_sub_spec)(kk)/procDataToPlot.txfrq*1e6), 'Pre', 'Color', colormap.LightAccent);
            text(ax_drift, length(crDriftPost)*1.05, crDriftPost(end)-(MRSCont.QM{VoxelIndex}.freqShift.(which_sub_spec)(kk)/procDataToPlot.txfrq*1e6), 'Post', 'Color', colormap.Foreground);
            set(ax_drift, 'YLim', [3.028-0.1 3.028+0.1]);
            yticks([3.028-0.08 3.028-0.04 3.028 3.028+0.04 3.028+0.08]);
            yticklabels({'2.94' '2.98' '3.02' '3.06' '3.10'});
            x = xlim;
            plot(ax_drift, [x(1) x(2)], [3.028 3.028],'LineStyle', ':', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028-0.04 3.028-0.04],'LineStyle', '--', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028+0.04 3.028+0.04],'LineStyle', '--', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            hold(ax_drift, 'off');
        else 
            x = xlim;
            y = yLims;
            text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color', colormap.Foreground);
        end
    else
        x = xlim;
        y = yLims;
        text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color', colormap.Foreground);
    end
elseif (isfield(MRSCont.flags,'isMRSI') && (MRSCont.flags.isMRSI == 1))
    if isfield(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.drift.pre, which_sub_spec)
        if length(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.drift.pre.(which_sub_spec){kk}) > 1
            crDriftPre = MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.drift.pre.(which_sub_spec){kk} + MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift.(which_sub_spec)(kk)/applyDataToPlot.txfrq*1e6;
            crDriftPost = MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.drift.post.(which_sub_spec){kk} + MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift.(which_sub_spec)(kk)/applyDataToPlot.txfrq*1e6;
            hold(ax_drift, 'on');
            colors = ones(length(crDriftPre),1).*colormap.Foreground;
            for dots = 1 : length(crDriftPre)
                colors(dots,1) = colors(dots,1) + (1 - colors(dots,1)) * (1-weights(dots));
                colors(dots,2) = colors(dots,2) + (1 - colors(dots,2)) * (1-weights(dots));
                colors(dots,3) = colors(dots,3) + (1 - colors(dots,3)) * (1-weights(dots));
            end
            scatter(ax_drift, [1:length(crDriftPre)],crDriftPre'-(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift.(which_sub_spec)(kk)/procDataToPlot.txfrq*1e6),36,ones(length(crDriftPre),1).*colormap.LightAccent);
            scatter(ax_drift, [1:length(crDriftPost)],crDriftPost'-(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift.(which_sub_spec)(kk)/procDataToPlot.txfrq*1e6),36,colors,'filled','MarkerEdgeColor',colormap.Foreground);    

            text(ax_drift, length(crDriftPre)*1.05, crDriftPre(end)-(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift.(which_sub_spec)(kk)/procDataToPlot.txfrq*1e6), 'Pre', 'Color', colormap.LightAccent);
            text(ax_drift, length(crDriftPost)*1.05, crDriftPost(end)-(MRSCont.QM{VoxelIndex(1), VoxelIndex(2)}.freqShift.(which_sub_spec)(kk)/procDataToPlot.txfrq*1e6), 'Post', 'Color', colormap.Foreground);
            set(ax_drift, 'YLim', [3.028-0.1 3.028+0.1]);
            yticks([3.028-0.08 3.028-0.04 3.028 3.028+0.04 3.028+0.08]);
            yticklabels({'2.94' '2.98' '3.02' '3.06' '3.10'});
            x = xlim;
            plot(ax_drift, [x(1) x(2)], [3.028 3.028],'LineStyle', ':', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028-0.04 3.028-0.04],'LineStyle', '--', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028+0.04 3.028+0.04],'LineStyle', '--', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            hold(ax_drift, 'off');
        else 
            x = xlim;
            y = yLims;
            text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color', colormap.Foreground);
        end
    else
        x = xlim;
        y = yLims;
        text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color', colormap.Foreground);
    end
else
    if isfield(MRSCont.QM.drift.pre, which_sub_spec) && ~strcmp(which_spec,'mm')
        if length(MRSCont.QM.drift.pre.(which_sub_spec){kk}) > 1
            crDriftPre = MRSCont.QM.drift.pre.(which_sub_spec){kk} + MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/applyDataToPlot.txfrq*1e6;
            crDriftPost = MRSCont.QM.drift.post.(which_sub_spec){kk} + MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/applyDataToPlot.txfrq*1e6;
            hold(ax_drift, 'on');
            colors = ones(length(crDriftPre),1).*colormap.Foreground;
            for dots = 1 : length(crDriftPre)
                colors(dots,1) = colors(dots,1) + (1 - colors(dots,1)) * (1-weights(dots,1));
                colors(dots,2) = colors(dots,2) + (1 - colors(dots,2)) * (1-weights(dots,1));
                colors(dots,3) = colors(dots,3) + (1 - colors(dots,3)) * (1-weights(dots,1));
            end
            scatter(ax_drift, [1:length(crDriftPre)],crDriftPre'-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6),36,ones(length(crDriftPre),1).*colormap.LightAccent);
            scatter(ax_drift, [1:length(crDriftPost)],crDriftPost'-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6),36,colors,'filled','MarkerEdgeColor',colormap.Foreground);    

            text(ax_drift, length(crDriftPre)*1.05, crDriftPre(end)-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6), 'Pre', 'Color', colormap.LightAccent);
            text(ax_drift, length(crDriftPost)*1.05, crDriftPost(end)-(MRSCont.QM.freqShift.(which_spec)(1,kk,SubSpectraIndex)/procDataToPlot.txfrq*1e6), 'Post', 'Color', colormap.Foreground);
            set(ax_drift, 'YLim', [3.028-0.1 3.028+0.1]);
            yticks([3.028-0.08 3.028-0.04 3.028 3.028+0.04 3.028+0.08]);
            yticklabels({'2.94' '2.98' '3.02' '3.06' '3.10'});
            x = xlim;
            plot(ax_drift, [x(1) x(2)], [3.028 3.028],'LineStyle', ':', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028-0.04 3.028-0.04],'LineStyle', '--', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            plot(ax_drift, [x(1) x(2)], [3.028+0.04 3.028+0.04],'LineStyle', '--', 'Color', colormap.Foreground, 'LineWidth', 0.5);
            hold(ax_drift, 'off');
        else 
            x = xlim;
            y = yLims;
            text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color', colormap.Foreground);
        end
    else
        x = xlim;
        y = yLims;
        text(ax_drift, x(2)/6, y(2)/2, 'No drift data available','Color', colormap.Foreground);
        set(ax_drift, 'YLim', yLims);
    end
end
    xlabel(ax_drift, 'Averages', 'Color', colormap.Foreground);
    ylabel(ax_drift, 'Cr frequency (ppm)', 'Color', colormap.Foreground);
    title(ax_drift, 'chemical shift drift', 'Color', colormap.Foreground);        


%%% 8. DESIGN FINETUNING %%%
% Adapt common style for all axes
axs = {ax_raw, ax_aligned, ax_proc, ax_drift};
for ll = 1:length(axs)
    gca = axs{ll};
    set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XMinorTick', 'On');
    set(gca, 'FontSize', 16);

    % Black axes, white background
    if ~MRSCont.flags.isGUI
        set(gca, 'XColor', 'k');
        set(gca, 'Color', 'w');
        % If no y caption, remove y axis
        if isempty(gca.YLabel.String)
            set(gca, 'YColor', 'w');
        else
            set(gca, 'YColor', 'k');
        end
    else
        set(gca, 'XColor', colormap.Foreground);
        set(gca, 'Color', colormap.Background);
        % If no y caption, remove y axis
        if isempty(gca.YLabel.String)
            set(gca, 'YColor', colormap.Background);
        else
            set(gca, 'YColor', colormap.Foreground);
        end        
    end

end

gcf = out;
set(gcf, 'Color', MRSCont.colormap.Background);        
box off;


%%% 9. ADD OSPREY LOGO %%%
% Add to the printout, but not if displayed in the GUI.
if ~MRSCont.flags.isGUI
    [I, map] = imread('osprey.gif','gif');
    axes(out, 'Position', [0, 0.85, 0.15, 0.15*11.63/14.22]);
    imshow(I, map);
    axis off;
end
end

   