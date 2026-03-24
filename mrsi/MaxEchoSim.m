% Simulation Scenario 2

load('/Volumes/Samsung/working/OspreyMRSIbeta/basissets/BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA140_noMM.mat')
load('/Volumes/Samsung/working/MaxEcho/tumor-01/derivatives/ShortTE/jobShortTE.mat')
[resBasisSet] = fit_resampleBasis(MRSCont.processed.metab{1, 1}, BASIS);
dummyDataFID = MRSCont.processed.metab{1, 1};
dummyDataFID.fids = resBasisSet.fids(:,17);
dummyDataFID.specs = resBasisSet.specs(:,17);
op_plotspec(dummyDataFID)
op_plotfid(dummyDataFID)
dummyDataFID.te = 140;

load('/Volumes/Samsung/working/OspreyMRSIbeta/basissets/BASIS_Philips_UnEdited_se_MRSI_PRESS_GABA140_noMM_maxEcho.mat');
[resBasisSet] = fit_resampleBasis(MRSCont.processed.metab{1, 1}, BASIS);
dummyDataECHO = MRSCont.processed.metab{1, 1};
dummyDataECHO.fids = resBasisSet.fids(:,17);
dummyDataECHO.specs = resBasisSet.specs(:,17);
tleft = 140-75;        % How much time before echo
tleft = tleft/1000;   % Convert to ms
pointsToPick = (tleft/resBasisSet.dwelltime)-1;
dummyDataECHO.t = dummyDataECHO.t - resBasisSet.dwelltime*pointsToPick;  % Get time vector
dummyDataECHO.te = 140;
op_plotspec(dummyDataECHO)
op_plotfid(dummyDataECHO)

%%


%% Let's try adding some LB
noise_sd = 0.0000001;
T2 = 300/1000; %NAA
R2_L = 1/T2; % or R2p_L = sigma-R2_L

sigma = 6*pi;
R2p_L = sigma-R2_L;


DataFID_LB_LG = op_echo_filter(dummyDataFID, R2_L, R2p_L, sigma);
DataFID_LB_LG = op_addNoise(DataFID_LB_LG,noise_sd);
[SNR_FID]=op_getSNR(DataFID_LB_LG);

figure
tiledlayout(2,1,'TileSpacing','tight')
nexttile
plot(DataFID_LB_LG.ppm,abs(DataFID_LB_LG.specs))
nexttile
plot(DataFID_LB_LG.t,real(DataFID_LB_LG.fids))

DataECHO_LB_LG = op_echo_filter(dummyDataECHO, R2_L, R2p_L, sigma);
DataECHO_LB_LG = op_addNoise(DataECHO_LB_LG,noise_sd);
[SNR_ECHO]=op_getSNR(DataECHO_LB_LG);

figure
tiledlayout(2,1,'TileSpacing','tight')
nexttile
plot(DataECHO_LB_LG.ppm,abs(DataECHO_LB_LG.specs))
nexttile
plot(DataECHO_LB_LG.t,real(DataECHO_LB_LG.fids))

%% Now lets run a matrix of different combinations
noise_sd = 0.0000001;
sigma_values = linspace(2,8,100)*pi;
T2_values = linspace(250,300,100)/1000;
R2_L_values = 1./T2_values;
R2p_L_values = sigma_values - R2_L_values;
T2p_L_values = 1./R2p_L_values;


EchoToFIDSNR = zeros(100,100);

for rr = 1 : 100
    for kk = 1 : 100

        DataFID_LB_LG = op_echo_filter(dummyDataFID, R2_L_values(kk), R2p_L_values(kk), sigma_values(rr));
        DataFID_LB_LG = op_addNoise(DataFID_LB_LG,noise_sd);
        NAAwindow=DataFID_LB_LG.specs(DataFID_LB_LG.ppm>1.9 & DataFID_LB_LG.ppm<2.1);
        ppmwindow=DataFID_LB_LG.ppm(DataFID_LB_LG.ppm>1.9 & DataFID_LB_LG.ppm<2.1);
        
        maxNAA_index=find(abs(NAAwindow)==max(abs((NAAwindow))));
        maxNAA=abs(NAAwindow(maxNAA_index));
    
        SNR_FID(kk)=maxNAA/noise_sd;

        DataECHO_LB_LG = op_echo_filter(dummyDataECHO, R2_L_values(kk), R2p_L_values(kk), sigma_values(rr));
        DataECHO_LB_LG = op_addNoise(DataECHO_LB_LG,noise_sd);
        NAAwindow=DataECHO_LB_LG.specs(DataECHO_LB_LG.ppm>1.9 & DataECHO_LB_LG.ppm<2.1);
        ppmwindow=DataECHO_LB_LG.ppm(DataECHO_LB_LG.ppm>1.9 & DataECHO_LB_LG.ppm<2.1);
        
        maxNAA_index=find(abs(NAAwindow)==max(abs((NAAwindow))));
        maxNAA=abs(NAAwindow(maxNAA_index));
    
        SNR_ECHO(kk)=maxNAA/noise_sd;

        EchoToFIDSNR(kk,rr) = SNR_ECHO/SNR_FID;
    
    end
end

EchoToFIDSNR = flip(EchoToFIDSNR,1);
EchoToFIDSNR = flip(EchoToFIDSNR,2);
T2_values = flip(T2_values,2);
T2p_L_values = flip(T2p_L_values,2);
sigma_values = flip(sigma_values)/pi;
%%
figure, imagesc(EchoToFIDSNR)
colormap(viridis)
cbar = colorbar;
axis image
set(gca,'yticklabels',num2cell(round(T2_values(10:10:100),2)*1000))
set(gca,'xticklabels',num2cell(round(T2p_L_values(10:10:100),2)*1000))
% yticks(T2_values(10:10:100)*1000)
ylabel('T_2 (ms)')
xlabel('T_2'' = 1/(\pi * FWHM - R_2) (ms)')
cbar.Label.String = 'SNR_{Max Echo}/SNR_{FID}';
set(gca,'TickDir','out');
set(gcf,'Color',[1 1 1])
set(gcf,'renderer','painters')
saveas(gcf,fullfile('/Volumes/Samsung/abstracts/ISMRM2026/MaxEcho','SimulationScenario_2.pdf'),'pdf')


figure, imagesc(EchoToFIDSNR)
colormap(viridis)
cbar = colorbar;
axis image
set(gca,'yticklabels',num2cell(round(T2_values(10:10:100),2)*1000))
set(gca,'xticklabels',num2cell(round(sigma_values(10:10:100),2)))
% yticks(T2_values(10:10:100)*1000)
ylabel('T_2 (ms)')
xlabel('FWHM (Hz)')
cbar.Label.String = 'SNR_{Max Echo}/SNR_{FID}';
set(gca,'TickDir','out');
set(gcf,'Color',[1 1 1])
set(gcf,'renderer','painters')
saveas(gcf,fullfile('/Volumes/Samsung/abstracts/ISMRM2026/MaxEcho','SimulationScenario_2_FWHM.pdf'),'pdf')
%% Now some SNR and tNAA/MM plots
% Load the 15 ms spin echo data here for some amplitudes etc
load('/Volumes/Samsung/T7Shield/working/MRSI/derivatives/conv/sub_6/visit_1/final_metab_mets_cc_spline_optim_m14_ind_lor.mat');
load('/Volumes/Samsung/working/MaxEcho/tumor-01/derivatives/ShortTE/jobShortTE.mat')
center_voxel = final_metab_mets_cc{7,8,1};
%%
TE_values = [15:center_voxel.Data.DwellTime*1000:350];

for kk = 1 : length(TE_values)
    TE = TE_values(kk);

    noise_sd = 0.000000003;
    NAA_fid = center_voxel.BasisSets.fids(:,13);
    MM_fid = center_voxel.BasisSets.fids(:,24);
    
    NAA_ampl = (center_voxel.Model{1, 1}.parsOut.metAmpl(12) + center_voxel.Model{1, 1}.parsOut.metAmpl(12)) .* exp(-(TE-15)/300);
    MM_ampl = center_voxel.Model{1, 1}.parsOut.metAmpl(19).* exp(-(TE-15)/14);
    
    NAA_ampl_val(kk) = NAA_ampl;
    MM_ampl_val(kk) = MM_ampl;
    ratioNAAtoMM(kk) = max(NAA_ampl)/max(MM_ampl);
    
    T2 = 300/1000; %NAA
    R2_L = 1/T2; % or R2p_L = sigma-R2_L
    sigma = 5*pi;
    R2p_L = sigma-R2_L;
    
    dummyDataFID_tNAA = MRSCont.processed.metab{1, 1};
    dummyDataFID_tNAA.fids = NAA_fid * NAA_ampl;
    dummyDataFID_tNAA.specs = fftshift(fft(dummyDataFID_tNAA.fids,[],dummyDataFID_tNAA.dims.t),dummyDataFID_tNAA.dims.t);
    dummyDataFID_tNAA.te = TE;
    DataFID_tNAA_filter = op_echo_filter(dummyDataFID_tNAA, R2_L, R2p_L, sigma);
    DataFID_tNAA_filter = op_addNoise(DataFID_tNAA_filter,noise_sd);
    
    
    
    T2 = 14/1000; %MM20
    R2_L = 1/T2; % or R2p_L = sigma-R2_L
    sigma = 5*pi;
    R2p_L = sigma-R2_L;
    
    dummyDataFID_MM = MRSCont.processed.metab{1, 1};
    dummyDataFID_MM.fids = MM_fid * MM_ampl;
    dummyDataFID_MM.specs = fftshift(fft(dummyDataFID_MM.fids,[],dummyDataFID_MM.dims.t),dummyDataFID_MM.dims.t);
    dummyDataFID_MM.te = TE;
    DataFID_MM_filter = op_echo_filter(dummyDataFID_MM, R2_L, R2p_L, sigma);
    DataFID_MM_filter = op_addNoise(DataFID_MM_filter,noise_sd);
    
    
    DataFID = DataFID_tNAA_filter;
    DataFID.fids = DataFID.fids + DataFID_MM_filter.fids;
    DataFID.specs = DataFID.specs + DataFID_MM_filter.specs;
    
    NAAwindow=DataFID.specs(DataFID.ppm>1.9 & DataFID.ppm<2.1);
    ppmwindow=DataFID.ppm(DataFID.ppm>1.9 & DataFID.ppm<2.1);
    
    maxNAA_index=find(abs(NAAwindow)==max(abs((NAAwindow))));
    maxNAA=abs(NAAwindow(maxNAA_index));

    SNR_FID(kk)=maxNAA/noise_sd;
    
    % figure, plot(DataFID_tNAA_filter.ppm,abs(DataFID_tNAA_filter.specs)); hold on
    % plot(DataFID_MM_filter.ppm,abs(DataFID_MM_filter.specs));
    % plot(DataFID.ppm,abs(DataFID.specs));
    % 
    % figure, plot(DataFID.t,real(DataFID.fids)); 
    
    % Now for a MaxEcho readout
    % The maximal time before echo is TE/2 is the gradient duration of 5 ms
    
    
    tadd = TE/2-5;
    tadd = tadd/1000;
    
    pointsToAdd = round(tadd/center_voxel.Data.DwellTime);
    
    NAA_echo = [conj(flipud(NAA_fid(2:pointsToAdd)));NAA_fid(1:end-(pointsToAdd-1))];
    MM_echo = [conj(flipud(MM_fid(2:pointsToAdd)));MM_fid(1:end-(pointsToAdd-1))];

    
    T2 = 300/1000; %NAA
    R2_L = 1/T2; % or R2p_L = sigma-R2_L
    sigma = 5*pi;
    R2p_L = sigma-R2_L;
    
    dummyDataECHO_tNAA = MRSCont.processed.metab{1, 1};
    dummyDataECHO_tNAA.fids = NAA_echo * NAA_ampl;
    dummyDataECHO_tNAA.specs = fftshift(fft(dummyDataECHO_tNAA.fids,[],dummyDataECHO_tNAA.dims.t),dummyDataECHO_tNAA.dims.t);
    dummyDataECHO_tNAA.te = TE;
    dummyDataECHO_tNAA.t = dummyDataECHO_tNAA.t-tadd;
    DataECHO_tNAA_filter = op_echo_filter(dummyDataECHO_tNAA, R2_L, R2p_L, sigma);
    DataECHO_tNAA_filter = op_addNoise(DataECHO_tNAA_filter,noise_sd);
    
    
    T2 = 14/1000; %MM20
    R2_L = 1/T2; % or R2p_L = sigma-R2_L
    sigma = 5*pi;
    R2p_L = sigma-R2_L;
    
    dummyDataECHO_MM = MRSCont.processed.metab{1, 1};
    dummyDataECHO_MM.fids = MM_echo * MM_ampl;
    dummyDataECHO_MM.specs = fftshift(fft(dummyDataECHO_MM.fids,[],dummyDataECHO_MM.dims.t),dummyDataECHO_MM.dims.t);
    dummyDataECHO_MM.te = TE;
    dummyDataECHO_MM.t = dummyDataECHO_MM.t-tadd;
    DataECHO_MM_filter = op_echo_filter(dummyDataECHO_MM, R2_L, R2p_L, sigma);
    DataECHO_MM_filter = op_addNoise(DataECHO_MM_filter,noise_sd);
    
    DataECHO = DataECHO_tNAA_filter;
    DataECHO.fids = DataECHO.fids + DataECHO_MM_filter.fids;
    DataECHO.specs = DataECHO.specs + DataECHO_MM_filter.specs;
    
    
    NAAwindow=DataECHO.specs(DataECHO.ppm>1.9 & DataECHO.ppm<2.1);
    ppmwindow=DataECHO.ppm(DataECHO.ppm>1.9 & DataECHO.ppm<2.1);
    
    maxNAA_index=find(abs(NAAwindow)==max(abs((NAAwindow))));
    maxNAA=abs(NAAwindow(maxNAA_index));

    SNR_ECHO(kk)=maxNAA/noise_sd;
    
    
    % figure, plot(DataECHO_tNAA_filter.ppm,abs(DataECHO_tNAA_filter.specs)); hold on
    % plot(DataECHO_MM_filter.ppm,abs(DataECHO_MM_filter.specs));
    % plot(DataECHO.ppm,abs(DataECHO.specs));
    % 
    % figure, plot(DataECHO.t,real(DataECHO.fids)); 
    % close all
end
%%
figure
yyaxis left
plot(TE_values,SNR_FID/SNR_FID(1)), hold on
plot(TE_values,SNR_ECHO/SNR_ECHO(1))
plot(TE_values,NAA_ampl_val/NAA_ampl_val(1))
plot(TE_values,MM_ampl_val/MM_ampl_val(1))
ylabel('SNR or amplitude relative to TE = 15 ms')
yyaxis right
plot(TE_values,SNR_ECHO./SNR_FID)
ylabel('SNR_{ECHO}/SNR_{FID}')
set(gca,'XLim',[15 TE_values(end)],'TickDir','out')
xlabel('TE (ms)')
set(gcf,'Color',[1 1 1])
set(gcf,'renderer','painters')
saveas(gcf,fullfile('/Volumes/Samsung/abstracts/ISMRM2026/MaxEcho','SimulationScenario_1.pdf'),'pdf')