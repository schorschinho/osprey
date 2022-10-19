function out = osp_plotMRSImap(MRSCont, kk,slice, nominator_spec,denominator_spec, nominator, denominator, upsample,figTitle,mask)
%% out = osp_plotFit(MRSCont, kk, which, stagFlag, xlab, ylab, figTitle)
%   Creates a figure showing data stored in an Osprey data container, as
%   well as the fit to it, the baseline, the residual, and contributions
%   from the individual metabolites.
%
%   USAGE:
%       out = osp_plotFit(MRSCont, kk, which, GUI, conc, stagFlag, xlab, ylab, figTitle)
%
%   OUTPUTS:
%       out     = MATLAB figure handle
%
%   OUTPUTS:
%       MRSCont  = Osprey data container.
%       kk       = Index for the kk-th dataset (optional. Default = 1)
%       which    = String for the spectrum to plot (optional)
%                   OPTIONS:    'off' (default)
%                               'diff1'
%                               'diff2'
%                               'sum'
%                               'ref'
%                               'w'
%                                 'mm' re_mm
%
%       xlab      = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
%       ylab      = label for the y-axis (optional.  Default = '');
%       figTitle  = label for the title of the plot (optional.  Default = '');
%
%   AUTHOR:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2019-10-02)
%       goeltzs1@jhmi.edu
%
%   HISTORY:
%       2019-10-02: First version of the code.

% Check that OspreyFit has been run before
if ~MRSCont.flags.didFit
    error('Trying to plot fitted data, but fit has not been performed yet. Run OspreyFit first.')
end


%%% 1. PARSE INPUT ARGUMENTS %%%
% Get the fit method and style
fitMethod   = MRSCont.opts.fit.method;
fitStyle    = MRSCont.opts.fit.style;
% Fall back to defaults if not provided
if nargin < 10
mask = 0;
    if nargin < 9
        figTitle = '';
        if nargin < 8
            upsample = 1;
            if nargin<7
                denominator = [];
                if nargin<6
                    nominator = {'Cr', 'PCr'}; 
                    if nargin<5
                        denominator_spec = 'none'; 
                        if nargin < 4
                            nominator_spec = 'off';
                            if nargin < 3
                                slice = 1;
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
if nargin<10    
    [~,filen,ext] = fileparts(MRSCont.files{kk});
    nominator_name = '';
    denominator_name = '';
    for t = 1 : length(nominator)
        nominator_name = [nominator_name nominator{t}];
        if t < length(nominator)
            nominator_name = [nominator_name '+'];
        end
    end
    if ~isempty(denominator)
        for t = 1 : length(denominator)
            denominator_name = [denominator_name denominator{t}];
            if t < length(denominator)
            denominator_name = [denominator_name '+'];
            end
        end
        figTitle = sprintf([filen ext '\n' fitMethod ' ' fitStyle ': ' nominator_name '/' denominator_name  ]); 
    else
        figTitle = sprintf([filen ext '\n' fitMethod ' ' fitStyle ': ' nominator_name ' raw amplitudes']); 
    end    
end


%%% 2. EXTRACT DATA TO PLOT %%%
% Extract processed spectra and fit parameters
nominator_map = zeros(size(MRSCont.quantify.amplMets{kk}.(nominator_spec).(nominator{1})));

if  (MRSCont.flags.isMRSI == 1)
    for nom = 1 : length(nominator)
        nominator_map = nominator_map + MRSCont.quantify.amplMets{kk}.(nominator_spec).(nominator{nom});
    end
    if (~isempty(denominator) && ~isempty(denominator_spec)) && ~strcmp(denominator_spec, 'w')
        denominator_map = zeros(size(MRSCont.quantify.amplMets{kk}.(denominator_spec).(denominator{1})));
        for denom = 1 : length(denominator)
            denominator_map = denominator_map + MRSCont.quantify.amplMets{kk}.(denominator_spec).(denominator{denom});
        end
        if (sum(size(nominator_map) == size(denominator_map)) == 0)
            sz_denominator_map = size(denominator_map);
            sz_nominator_map = size(nominator_map); 
            ratio = sz_denominator_map./sz_nominator_map;
            scale = prod(ratio);
            denominator_map = imresize(denominator_map,sz_nominator_map);
            denominator_map = denominator_map/scale;
        end
        map = nominator_map ./ denominator_map;
        map(denominator_map == 0) = 0;
    else
        map = nominator_map;
    end    
end



%%% 4. SET UP FIGURE LAYOUT %%%
% Generate a new figure and keep the handle memorized
canvasSize  = get(0,'defaultfigureposition');
if ~MRSCont.flags.isGUI
    out = figure('Position', canvasSize);
else
    out = figure('Position', canvasSize,'Visible','off');
end

sz_map = size(map);

map = map(:,:,slice);
map = rot90(map);



if strcmp(denominator_spec, 'w')
    denominator_map = MRSCont.quantify.SH2O_factor{kk};
    denominator_map = denominator_map(:,:,slice);
    denominator_map = rot90(denominator_map);
    [T1_GM, T1_WM, T2_GM, T2_WM] = lookUpRelaxTimes(nominator);
    T1_Metab = mean([T1_GM T1_WM]);
    T2_Metab = mean([T2_GM T2_WM]);
    if strcmp(nominator_spec, 'off')
        metsTR  = MRSCont.processed.A{kk}.tr * 1e-3;
        metsTE  = MRSCont.processed.A{kk}.te * 1e-3;
    end
    if strcmp(nominator_spec, 'diff1')
        metsTR  = MRSCont.processed.diff1{kk}.tr * 1e-3;
        metsTE  = MRSCont.processed.diff1{kk}.te * 1e-3;
    end
    T1_Factor = 1./ (1-exp(-metsTR./T1_Metab));
    T2_Factor = 1 ./ exp(-metsTE./T2_Metab);
    map = map ./ denominator_map .* T1_Factor .* T2_Factor;
end


if mask
    mask = MRSCont.mask{kk};
    mask = mask(:,:,slice);
    mask(mask > 1) = 1;
    mask = rot90(mask);
    if ~(sum(size(mask) == size(map)) == 0)
        map = mask .* map;
    end
end

map(isnan(map)) = 0;
map(isinf(map)) = 0;

if upsample > 1
    map = imresize(map,sz_map(1:2) .* upsample);
end


map_mean = mean(map(mask==1),'all');
map_std = std(map(mask==1));


colormap = viridis(100);

% heatmap(map,'Colormap',colormap);

map(map > map_mean +2.5* map_std) = 0;
map_mean = mean(map(mask==1),'all');

heatmap(fliplr(map),'Colormap',colormap);
% heatmap(map,'Colormap',gray);

% caxis(out.Children,[0 map_mean]);
% colorbar off

%%% 7. DESIGN FINETUNING %%%
% Adapt common style for all axes

set(gca, 'FontSize', 16);
ax = gca;
% 
iniXLabel = get(ax,'XDisplayLabels');
iniYLabel = get(ax,'YDisplayLabels');

for l = 1 : length(iniXLabel)
        iniXLabel{l,1} = '';
end

for l = 1 : length(iniYLabel)
        iniYLabel{l,1} = '';
end
set(ax, 'XDisplayLabels',iniXLabel)
set(ax, 'YDisplayLabels',iniYLabel)
set(ax, 'FontSize',12)

if ~MRSCont.flags.isGUI
    % Black axes, white background
    title(figTitle);
else
    title(figTitle);
end

set(gcf, 'Color', MRSCont.colormap.Background);   

set(gcf,'Units','Normalized');
set(gcf,'Position',[1.1250 0.1411 0.66 1]);
%%% 8. ADD OSPREY LOGO AND TIGHTEN FIGURE %%%
% if ~MRSCont.flags.isGUI
%     [I, map] = imread('osprey.gif','gif');
%     axes(out, 'Position', [0, 0.9, 0.1, 0.1*11.63/14.22]);
%     imshow(I, map);
%     axis off;
% end


%% Generate nifti map
vol_mask = MRSCont.coreg.vol_mask{kk}{slice};
% maskFile = vol_mask.fname;
% index_maskFileOut            = vol_mask.fname;  
% vol_image = spm_vol(vol_mask.fname);

mask_res = MRSCont.coreg.MRSI_res_mask{kk}{slice};
metab_res = mask_res;
for z = 1 : size(mask_res,3)
     y_ind = size(map,1);
    for y = 1 : size(mask_res,2)
        for x = 1 : size(mask_res,1)
            if mask_res(x,y,z) > 0
                metab_res(x,y,z) = map(y_ind,x);    
            end
        end
        if mask_res(x,y,z) > 0
            y_ind = y_ind -1;
        end
    end
 end
 metab_res = metab_res .* mask_res;

metab_res_int = zeros(vol_mask.dim);
if (vol_mask.dim(1) ~= size(metab_res,1)) && (vol_mask.dim(2) ~= size(metab_res,2))
    for z = 1 : vol_mask.dim(3)
        metab_res_int(:,:,z)= imresize(double(squeeze(metab_res(:,:,z))),[vol_mask.dim(1),vol_mask.dim(2)],'method','nearest');
    end
end
name = '';
for mm = 1 : length(nominator)
    name = [name ' ' nominator{mm}];
end
name = [name ' normalized to ' ];
for mm = 1 : length(denominator)
    name = [name ' ' denominator{mm}];
end
MRSI_map_filename            = strrep(vol_mask.fname,'VoxelMask_',[name '_']);    
vol_map.fname   = MRSI_map_filename;
vol_map.dim     = vol_mask.dim;
vol_map.dt      = vol_mask.dt;
vol_map.mat     = vol_mask.mat;

vol_map.descrip = ['metabolite map slice ' num2str(slice)];
% Write the SPM volume to disk
vol_map = spm_write_vol(vol_map,metab_res_int);
end

%%% Lookup function for metabolite relaxation times %%%
function [T1_GM, T1_WM, T2_GM, T2_WM] = lookUpRelaxTimes(metName)

% Look up table below
% T1 values for NAA, Glu, Cr, Cho, Ins from Mlynarik et al, NMR Biomed
% 14:325-331 (2001)
% T1 for GABA from Puts et al, J Magn Reson Imaging 37:999-1003 (2013)
% T2 values from Wyss et al, Magn Reson Med 80:452-461 (2018)
% T2 values are averaged between OCC and pACC for GM; and PVWM for WM
relax.Asc   = [1340 1190 (125+105)/2 172];
relax.Asp   = [1340 1190 (111+90)/2 148];
relax.Cr    = [1460 1240 (148+144)/2 166]; % 3.03 ppm resonance; 3.92 ppm signal is taken care of by -CrCH2 during fitting
relax.GABA  = [1310 1310 (102+75)/2 (102+75)/2]; % No WM estimate available; take GM estimate; both in good accordance with 88 ms reported by Edden et al
relax.Glc   = [1340 1190 (117+88)/2 155]; % Glc1: [1310 1310 (128+90)/2 156];
relax.Gln   = [1340 1190 (122+99)/2 168];
relax.Glu   = [1270 1170 (135+122)/2 124];
relax.Gly   = [1340 1190 (102+81)/2 152];
relax.GPC   = [1300 1080 (274+222)/2 218]; % This is the Choline singlet (3.21 ppm, tcho2 in the paper); glycerol is tcho: [1310 1310 (257+213)/2 182]; % choline multiplet is tcho1: [1310 1310 (242+190)/2 178];
relax.GSH   = [1340 1190 (100+77)/2 145]; % This is the cysteine signal (GSH1 in the paper), glycine is GSH: [1310 1310 (99+72)/2 145]; % glutamate is GSH2: [1310 1310 (102+76)/2 165];
relax.Lac   = [1340 1190 (110+99)/2 159];
relax.Ins   = [1230 1010 (244+229)/2 161];
relax.NAA   = [1470 1350 (253+263)/2 343]; % This is the 2.008 ppm acetyl signal (naa in the paper); aspartyl is naa1: [1310 1310 (223+229)/2 310];
relax.NAAG  = [1340 1190 (128+107)/2 185]; % This is the 2.042 ppm acetyl signal (naag in the paper); aspartyl is naag1: [1310 1310 (108+87)/2 180]; % glutamate is NAAG2: [1310 1310 (110+78)/2 157];
relax.PCh   = [1300 1080 (274+221)/2 213]; % This is the singlet (3.20 ppm, tcho4 in the paper); multiple is tcho3: [1310 1310 (243+191)/2 178];
relax.PCr   = [1460 1240 (148+144)/2 166]; % 3.03 ppm resonance; 3.92 ppm signal is taken care of by -CrCH2 during fitting; same as Cr
relax.PE    = [1340 1190 (119+86)/2 158];
relax.Scy   = [1340 1190 (125+107)/2 170];
relax.Tau   = [1340 1190 (123+102)/2 (123+102)/2]; % No WM estimate available; take GM estimate
relax.tNAA  = [1495 1385 155 155]; % Mean values from NAA + NAAG
relax.tCr  = [1740 1780 107 107]; % The singlet peak ar 3 ppm. 3.9 ppm peak values are [1240 1190 94 94] %T1 from Mlynarik et al. (2012)
relax.tCho  = [1510 1320  153 153]; % Entire molecule; T1 from Mlynarik et al. (2012) 10.1002/mrm.24352
relax.Glx  = [1625 1745 107 112]; % Mean values from Glu + Glx
    
T1_GM = [];
T1_WM = [];
T2_GM = [];
T2_WM = [];
for mm = 1 : length(metName)
% Check if metabolite name is in the look-up table
    if isfield(relax, metName{mm})
        T1_GM = [T1_GM relax.(metName{mm})(1) * 1e-3];
        T1_WM = [T1_WM relax.(metName{mm})(2) * 1e-3];
        T2_GM = [T2_GM relax.(metName{mm})(3) * 1e-3];
        T2_WM = [T2_WM relax.(metName{mm})(4) * 1e-3];
    else
        % If not, use an average
        T1_GM = [T1_GM 1340 * 1e-3];
        T1_WM = [T1_WM 1190 * 1e-3];
        T2_GM = [T2_GM 140 * 1e-3];
        T2_WM = [T2_WM 169 * 1e-3];
    end
end
T1_GM = mean(T1_GM);
T1_WM = mean(T1_WM);
T2_GM = mean(T2_GM);
T2_WM = mean(T2_WM);
end   