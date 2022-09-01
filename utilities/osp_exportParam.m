function osp_exportParam(MRSCont,path)
=========|=========|=========|=========|=========|=========|=========|=========|
% osp_exportParam is used to export the fitting parameters for each metabolite
% and the spectrum as a whole.
  
names = {'amplitude','lorentz','gauss','fshift','ph0','ph1','SNR','SNRstd','ecc'};
num_basis_fcns = shape(MRSCont.fit.basisSet{?})
% Define storage matrix as: [nDatasets,metab_amp,Lorentzians,Gaussians,fshifts,phi0,phi1,SNR,ecc]
% index_dict stores names/indices
% 18 + 8 = 26 basis functions
parameters = zeros(MRSCont.nDatasets(1), 26 * 3 + fshifts + phi0 + phi1 + SNR + ecc)

for kk = 1:MRSCont.nDatasets(1)
    parameters(kk,1:26)    = fit.results.metab.fitParams{?, kk}.amp;
    parameters(kk,27:52)   = fit.results.metab.fitParams{?, kk}.lorentzLB;
    parameters(kk,53:70)   = fit.results.metab.fitParams{?, kk}.gaussLB;
    parameters(kk,70:78)   = fit.results.metab.fitParams{?, kk}.gaussLBMM;
    parameters(kk,79:104)  = fit.results.metab.fitParams{?, kk}.fshifts; % rad2deg()
    parameters(kk,105)     = fit.results.metab.fitParams{?, kk}.ph0; % rad2deg()
    parameters(kk,106)     = fit.results.metab.fitParams{?, kk}.ph1; % rad2deg()
    parameters(kk,107)     = fit.results.metab.fitParams{?, kk}.SNR;
    parameters(kk,107)     = fit.results.metab.fitParams{?, kk}.SNRstd;
    parameters(kk,108)     = fit.results.metab.fitParams{?, kk}.ecc;
    
%     MMdef(kk,1) = FullDefSin.fit.results.metab.fitParams{1, kk}.ph1;
%     MMupd(kk,1) = FullUpdSep.fit.results.metab.fitParams{1, kk}.ph1;
%     MM1G(kk,1)  = FullDefSin.fit.results.metab.fitParams{2, kk}.ph1;
%     MM2G(kk,1)  = FullUpdSep.fit.results.metab.fitParams{2, kk}.ph1;
%     iMM1G(kk,1) = SingleGauss.fit.results.metab.fitParams{1, kk}.ph1;
%     iMM2G(kk,1) = SeparateGauss.fit.results.metab.fitParams{1, kk}.ph1;

%     MMdef(kk,2) = FullDefSin.fit.results.metab.fitParams{1, kk}.gaussLB;
%     MMupd(kk,2) = FullUpdSep.fit.results.metab.fitParams{1, kk}.gaussLB;
%     MM1G(kk,2)  = FullDefSin.fit.results.metab.fitParams{2, kk}.gaussLB;
%     MM2G(kk,2)  = FullUpdSep.fit.results.metab.fitParams{2, kk}.gaussLB;
%     iMM1G(kk,2) = SingleGauss.fit.results.metab.fitParams{1, kk}.gaussLB;
%     iMM2G(kk,2) = SeparateGauss.fit.results.metab.fitParams{1, kk}.gaussLB;
    
%     MMdef(kk,3) = 0;
%     MMupd(kk,3) = 0;
%     MM1G(kk,3)  = 0;
%     MM2G(kk,3)  = FullUpdSep.fit.results.metab.fitParams{2, kk}.gaussLBMM;
%     iMM1G(kk,3) = 0;
%     iMM2G(kk,3) = SeparateGauss.fit.results.metab.fitParams{1, kk}.gaussLBMM;

%     MMdef(kk,4) = mean(FullDefSin.fit.results.metab.fitParams{1, kk}.lorentzLB(1:18));
%     MMupd(kk,4) = mean(FullUpdSep.fit.results.metab.fitParams{1, kk}.lorentzLB(1:18));
%     MM1G(kk,4)  = mean(FullDefSin.fit.results.metab.fitParams{2, kk}.lorentzLB(1:18));
%     MM2G(kk,4)  = mean(FullUpdSep.fit.results.metab.fitParams{2, kk}.lorentzLB(1:18));
%     iMM1G(kk,4) = mean(SingleGauss.fit.results.metab.fitParams{1, kk}.lorentzLB(1:18));
%     iMM2G(kk,4) = mean(SeparateGauss.fit.results.metab.fitParams{1, kk}.lorentzLB(1:18));

%     MMdef(kk,5) = mean(FullDefSin.fit.results.metab.fitParams{1, kk}.lorentzLB(18:end));
%     MMupd(kk,5) = mean(FullUpdSep.fit.results.metab.fitParams{1, kk}.lorentzLB(18:end));
%     MM1G(kk,5)  = mean(FullDefSin.fit.results.metab.fitParams{2, kk}.lorentzLB(18:end));
%     MM2G(kk,5)  = mean(FullUpdSep.fit.results.metab.fitParams{2, kk}.lorentzLB(18:end));
%     iMM1G(kk,5) = mean(SingleGauss.fit.results.metab.fitParams{1, kk}.lorentzLB(18:end));
%     iMM2G(kk,5) = mean(SeparateGauss.fit.results.metab.fitParams{1, kk}.lorentzLB(18:end));

%     MMdef(kk,6) = FullDefSin.fit.results.metab.fitParams{1, kk}.ph0;
%     MMupd(kk,6) = FullUpdSep.fit.results.metab.fitParams{1, kk}.ph0;
%     MM1G(kk,6)  = FullDefSin.fit.results.metab.fitParams{2, kk}.ph0;
%     MM2G(kk,6)  = FullUpdSep.fit.results.metab.fitParams{2, kk}.ph0;
%     iMM1G(kk,6) = SingleGauss.fit.results.metab.fitParams{1, kk}.ph0;
%     iMM2G(kk,6) = SeparateGauss.fit.results.metab.fitParams{1, kk}.ph0;

end

% csv_paths = fullfile(path,{'def.csv','upd.csv','M_sin.csv','M_sep.csv','iM_sin.csv','iM_sep.csv'})

    
% writetable(array2table(MMdef,'VariableNames',names),csv_paths{1},'Delimiter',','); % Write table with tab delimiter
% writetable(array2table(MMupd,'VariableNames',names),csv_paths{2},'Delimiter',','); % Write table with tab delimiter
% writetable(array2table(MM1G, 'VariableNames',names),csv_paths{3},'Delimiter',','); % Write table with tab delimiter
% writetable(array2table(MM2G, 'VariableNames',names),csv_paths{4},'Delimiter',','); % Write table with tab delimiter
% writetable(array2table(iMM1G,'VariableNames',names),csv_paths{5},'Delimiter',','); % Write table with tab delimiter
% writetable(array2table(iMM2G,'VariableNames',names),csv_paths{6},'Delimiter',','); % Write table with tab delimiter

writetable(array2table(parameters, 'VariableNames', names), fullfile(path,'parameters.csv'),'Delimiter',',');
save(fullfile(path,'parameters.mat'),'names',names,'parameters',parameters','linenames',basisfcn_names)

  
end
