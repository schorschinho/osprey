% CSI data processing
% BY: Liangjie Lin
% PLACE: Johns Hopkins

%% add path and prepare
clc; close all; clear all;
path1 = 'ExampleData.SDAT'; % data with water suppression
path2 = 'ExampleData_NWS.SDAT'; % reference data without water suppression

ZF = 2048; % zerofill points
EP = 0.5;  % echo position 0~1

%% read Philips MRSI data
fidraw1 = mrs_readSDAT(path1); %% read Philips MRS/MRSI data
fidraw1 = permute(fidraw1,[3,2,1,4]); %% data rotate
fidraw2 = mrs_readSDAT(path2);
fidraw2 = permute(fidraw2,[3,2,1,4]);

lip1 = generator_lipid(2000, 512, ZF, 360, 480, 10, 30, ZF, EP); % generation of lipid basic matrix
wat1 = generator_water(2000, 512, ZF, 60, 10, 30, ZF, EP); % generation of water basic matrix

%% setting processing parameters
wf = 1; % window function/line broadening
ws = 1; % water shift/B0 correction
ls = 1; % lipid suppression
md = 1; % mask add
cs = 1; % chemical shift correction
bc = 1; % baseline correction

for ii = 2 % slice selection
    %% select one slice
    fid1 = fidraw1(:,:,:,ii);
    fid2 = fidraw2(:,:,:,ii);
    
    %% line boardening
    if wf == 1
        fid1 = mrs_apod3(fid1,2000,3,128);
        fid2 = mrs_apod3(fid2,2000,3,128);
    end
    
    %% Fourier transform to the spectral domain
    mrs1 = flip(fftshift(fft(fid1,ZF,3),3),3);
    mrs2 = flip(fftshift(fft(fid2,ZF,3),3),3);
    
    %% water frequency correction with water reference scan
    if  ws == 1
        [mrs2 mrs1] = watershift(mrs2,mrs1);
    end
    
    %% L2 water suppression
    mrs1 = watersup_sim(mrs1, real(wat1), 3);
    
    %% add mask to the csi data
    if md == 1
        meta_mask = maskmade(mrs1,80);
        N = size(mrs2);
        mrs1 = mrs1.*repmat(meta_mask,[1,1,N(3)]);
    end
    
    %% lipid suppression
    if ls == 1
        mrs1 = watersup_sim(mrs1, real(lip1), 3);
    end
    
    %% baseline correction
    if bc == 1
        mrs1 = mrs_baselinecorre(mrs1, ZF/5);
    end
    
    %% chemical shift correction after water and lipid suppression
    if cs == 1
        mrs1 = shiftcorre(abs(mrs1));
    end
    
    %% spectra plot
    xx = [10 15 15 15 20];
    yy = [8  8  12 16 12];
    for jj=1:5
        spectra_plot(xx(jj), yy(jj), 3, 3, 1, 5, 1, mrs1);
    end  
end


