function [MRSCont] =  osp_import_LCModel_spectra(MRSCont)

%% Get sequence type
if MRSCont.flags.isUnEdited
    name = 'off_A';
end

PathName = uipickfiles('FilterSpec',MRSCont.outputFolder,'REFilter', '\','NumFiles',1,'Prompt','Select the folder containing all LCModel .coord-files.');
PathName = PathName{1};
d1 = dir(fullfile(PathName));
files = length(d1);

j = 1;
for i = 1 : files
    if ~strcmp (d1(i,1).name, '.')
                if ~strcmp (d1(i,1).name, '..')
                    if ~strcmp (d1(i,1).name(1), '.')
                        if strcmp(d1(i,1).name(end), 'd')
                            if ~exist('metab')
                                [ spectra, ~, x_ppm, ~ ] = mrs_readLcmodelCOORD( fullfile(PathName,d1(i,1).name) );
                                temp_sz(1,j)= length(spectra);
                                temp_ppm{j} = x_ppm;
                                j = j + 1;
                            end
                        end
                    end
                end
    end
end
[max_point,max_ind] = max(temp_sz);
j = 1;
for i = 1 : files
    if ~strcmp (d1(i,1).name, '.')
                if ~strcmp (d1(i,1).name, '..')
                    if ~strcmp (d1(i,1).name(1), '.')
                        if strcmp(d1(i,1).name(end), 'd')
                            if ~exist('metab')
                                [ data, indivMetabs, ppm_x, info ] = mrs_readLcmodelCOORD( fullfile(PathName,d1(i,1).name) );
                                [ inital ] = mrs_readLcmodelPRINT( strrep(fullfile(PathName,d1(i,1).name),'.coord','.print') );
                                dataNames = info.metabolites;
                                idx_1  = find(strcmp(dataNames,'-CrCH2'));
                                if ~isempty(idx_1)
                                    dataNames{idx_1} = 'CrCH2';
                                end
                                scale = 1;
                                if length(data(:,1)) < max_point
                                    ppmRangeData        = temp_ppm{max_ind};
                                    ppmRangeDataToInt       = ppm_x;
                                    ppmIsInDataRange    = (ppmRangeDataToInt < ppmRangeData(1)) & (ppmRangeDataToInt > ppmRangeData(end));
                                    if sum(ppmIsInDataRange) == 0
                                        ppmIsInDataRange    = (ppmRangeDataToInt > ppmRangeData(1)) & (ppmRangeDataToInt < ppmRangeData(end));
                                    end
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.data      = interp1(ppmRangeDataToInt(ppmIsInDataRange), data(ppmIsInDataRange,1), ppmRangeData, 'pchip', 'extrap');
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.fit      = interp1(ppmRangeDataToInt(ppmIsInDataRange), data(ppmIsInDataRange,2), ppmRangeData, 'pchip', 'extrap');
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.baseline      = interp1(ppmRangeDataToInt(ppmIsInDataRange), data(ppmIsInDataRange,3), ppmRangeData, 'pchip', 'extrap');
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.ppm = ppmRangeData';
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.res = MRSCont.overview.LCModel.all_models.(name){1,j}.data-MRSCont.overview.LCModel.all_models.(name){1,j}.fit;
                                    for n = 1 : length(dataNames)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{n}])  = interp1(ppmRangeDataToInt(ppmIsInDataRange), indivMetabs(ppmIsInDataRange,n), ppmRangeData, 'pchip', 'extrap')-MRSCont.overview.LCModel.all_models.(name){1,j}.baseline;
                                    end
                                else
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.ppm = ppm_x';
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.data = data(:,1);
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.fit = data(:,2);
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.baseline = data(:,3);
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.res = MRSCont.overview.LCModel.all_models.(name){1,j}.data-MRSCont.overview.LCModel.all_models.(name){1,j}.fit;
                                    for n = 1 : length(dataNames)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{n}])  = indivMetabs(:,n)-data(:,3);
                                    end
                                end
                                if MRSCont.opts.fit.fitMM == 1
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM  = 0;
                                    idx_1  = find(strcmp(dataNames,'Lip13a'));
                                    if ~isempty(idx_1)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM  = MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    end
                                    idx_1  = find(strcmp(dataNames,'Lip13b'));
                                    if ~isempty(idx_1)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM  = MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    end
                                    idx_1  = find(strcmp(dataNames,'Lip09'));
                                    if ~isempty(idx_1)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM  = MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    end
                                    idx_1  = find(strcmp(dataNames,'MM09'));
                                    if ~isempty(idx_1)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM  = MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    end
                                    idx_1  = find(strcmp(dataNames,'Lip20'));
                                    if ~isempty(idx_1)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM  = MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    end
                                    idx_1  = find(strcmp(dataNames,'MM20'));
                                    if ~isempty(idx_1)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM  = MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    end
                                    idx_1  = find(strcmp(dataNames,'MM12'));
                                    if ~isempty(idx_1)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM  = MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    end
                                    idx_1  = find(strcmp(dataNames,'MM14'));
                                    if ~isempty(idx_1)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM  = MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    end
                                    idx_1  = find(strcmp(dataNames,'MM17'));
                                    if ~isempty(idx_1)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM  = MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    end
                                    idx_1  = find(strcmp(dataNames,'MMexp'));
                                    if ~isempty(idx_1)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM  = MRSCont.overview.LCModel.all_models.(name){1,j}.fittMM + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    end
                                end
                                idx_1  = find(strcmp(dataNames,'NAA'));
                                idx_2  = find(strcmp(dataNames,'NAAG'));
                                if ~isempty(idx_1) && ~isempty(idx_2)
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.fittNAA  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]) + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                else if ~isempty(idx_1)
                                     MRSCont.overview.LCModel.all_models.(name){1,j}.fittNAA  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    else if ~isempty(idx_2)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittNAA  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                        end
                                    end
                                end
                                idx_1  = find(strcmp(dataNames,'Cr'));
                                idx_2  = find(strcmp(dataNames,'PCr'));
                                if ~isempty(idx_1) && ~isempty(idx_2)
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.fittCr  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]) + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                else if ~isempty(idx_1)
                                     MRSCont.overview.LCModel.all_models.(name){1,j}.fittCr  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    else if ~isempty(idx_2)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittCr  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                        end
                                    end
                                end
                                idx_1  = find(strcmp(dataNames,'GPC'));
                                idx_2  = find(strcmp(dataNames,'PCh'));
                                if ~isempty(idx_1) && ~isempty(idx_2)
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.fittCho  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]) + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                else if ~isempty(idx_1)
                                     MRSCont.overview.LCModel.all_models.(name){1,j}.fittCho  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    else if ~isempty(idx_2)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fittCho  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                        end
                                    end
                                end
                                idx_1  = find(strcmp(dataNames,'Glu'));
                                idx_2  = find(strcmp(dataNames,'Gln'));
                                if ~isempty(idx_1) && ~isempty(idx_2)
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.fitGlx  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]) + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                else if ~isempty(idx_1)
                                     MRSCont.overview.LCModel.all_models.(name){1,j}.fitGlx = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    else if ~isempty(idx_2)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fitGlx  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                        end
                                    end
                                end
                                idx_1  = find(strcmp(dataNames,'Lip13a'));
                                idx_2  = find(strcmp(dataNames,'Lip13b'));
                                if ~isempty(idx_1) && ~isempty(idx_2)
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.fitLip13  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]) + MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                else if ~isempty(idx_1)
                                     MRSCont.overview.LCModel.all_models.(name){1,j}.fitLip13 = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_1}]);
                                    else if ~isempty(idx_2)
                                        MRSCont.overview.LCModel.all_models.(name){1,j}.fitLip13  = MRSCont.overview.LCModel.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                        end
                                    end
                                end                                
                                MRSCont.overview.LCModel.all_models.(name){1,j}.data = MRSCont.overview.LCModel.all_models.(name){1,j}.data/scale;
                                MRSCont.overview.LCModel.all_models.(name){1,j}.fit = MRSCont.overview.LCModel.all_models.(name){1,j}.fit/scale;
                                MRSCont.overview.LCModel.all_models.(name){1,j}.baseline = MRSCont.overview.LCModel.all_models.(name){1,j}.baseline/scale;
                                MRSCont.overview.LCModel.all_models.(name){1,j}.res = MRSCont.overview.LCModel.all_models.(name){1,j}.res/scale;
                                modelNames = fields(MRSCont.overview.LCModel.all_models.(name){1,j});
                                for n = 6 : length(fields(MRSCont.overview.LCModel.all_models.(name){1,j}))
                                    MRSCont.overview.LCModel.all_models.(name){1,j}.(modelNames{n})  = MRSCont.overview.LCModel.all_models.(name){1,j}.(modelNames{n})/scale;
                                end
                                %Find the ppm of the maximum peak magnitude within the given range:
                                ppmindex=find(MRSCont.overview.LCModel.all_models.(name){1,j}.data(MRSCont.overview.LCModel.all_models.(name){1,j}.ppm>1.9 & MRSCont.overview.LCModel.all_models.(name){1,j}.ppm<2.1)==max(MRSCont.overview.LCModel.all_models.(name){1,j}.data(MRSCont.overview.LCModel.all_models.(name){1,j}.ppm>1.9 & MRSCont.overview.LCModel.all_models.(name){1,j}.ppm<2.1)));
                                ppmrange=MRSCont.overview.LCModel.all_models.(name){1,j}.ppm(MRSCont.overview.LCModel.all_models.(name){1,j}.ppm>1.9 & MRSCont.overview.LCModel.all_models.(name){1,j}.ppm<2.1);
                                ppmmax=ppmrange(ppmindex);
                                refShift=(ppmmax-2.013);
                                MRSCont.overview.LCModel.all_models.(name){1,j}.ph0 = info.ph0;
                                MRSCont.overview.LCModel.all_models.(name){1,j}.ph1 = info.ph1;
                                MRSCont.overview.LCModel.all_models.(name){1,j}.iniph0 = inital.iniph0;
                                MRSCont.overview.LCModel.all_models.(name){1,j}.iniph1 = inital.iniph1;
                                MRSCont.overview.LCModel.all_models.(name){1,j}.iniFWHM = inital.iniFWHM;
                                 MRSCont.overview.LCModel.all_models.(name){1,j}.ppm = MRSCont.overview.LCModel.all_models.(name){1,j}.ppm - refShift;
                                j = j + 1;
                            end
                        end
                    end
                end
    end
end


%%
%%% 4. SORTING DATA  %%%

for kk = 1 : MRSCont.nDatasets
    data = MRSCont.overview.LCModel.all_models.(name){1,kk};
    modelNames = fields(data);
    Cr_height = max(data.data(data.ppm>2.9 & data.ppm<3.1));
    for n = 1 : length(modelNames)
        if ~strcmp('ppm',modelNames{n}) && ~strcmp('ph0',modelNames{n}) && ~strcmp('ph1',modelNames{n}) && ~strcmp('iniFWHM',modelNames{n}) && ~strcmp('iniph0',modelNames{n}) && ~strcmp('iniph1',modelNames{n})
            MRSCont.overview.LCModel.all_models.(name){1,kk}.(modelNames{n})  = MRSCont.overview.LCModel.all_models.(name){1,kk}.(modelNames{n})/Cr_height;
        end
    end    
end

%Exclude datasets
% if isfield(MRSCont, 'exclude')
%     if~isempty(MRSCont.exclude)
%         MRSCont.overview.groups(MRSCont.exclude) = [];
%     end
% end
for g = 1 : MRSCont.overview.NoGroups
    MRSCont.overview.LCModel.sort_fit.(['g_' num2str(g)]).(name) = MRSCont.overview.LCModel.all_models.(name)(1,MRSCont.overview.groups == g);
end
MRSCont.overview.LCModel.sort_fit.GMean.(name) = MRSCont.overview.LCModel.all_models.(name)(1,MRSCont.overview.groups > 0);



fits = fields(MRSCont.overview.LCModel.sort_fit.GMean.(name){1,1});
for kk = 2 : length(MRSCont.overview.LCModel.sort_fit.GMean.(name))
    newI = ismember(fields(MRSCont.overview.LCModel.sort_fit.GMean.(name){1,kk}),fits);
    if sum(newI) < length(newI)
        tempFits = fields(MRSCont.overview.LCModel.sort_fit.GMean.(name){1,kk});
        for f = 1 : length(newI)
            if newI(f) == 0
                fits{length(fits)+1} = tempFits{f};
            end
        end
    end
end

newI = ismember(fields(MRSCont.overview.LCModel.sort_fit.GMean.(name){1,kk}),{'fitCrCH2'});
if sum(newI) < length(newI)
    fits{length(fits)+1} = 'fitCrCH2';
end

names = fields(MRSCont.overview.LCModel.sort_fit);
for g = 1 : length(names)
        tempSubSpec = zeros(length(MRSCont.overview.LCModel.sort_fit.(names{g}).(name){1}),length(MRSCont.overview.LCModel.sort_fit.(names{g}).(name){1}.ppm));
        tempSubBaseline = tempSubSpec;
        tempSubRes = tempSubSpec;
        tempSubdata = tempSubSpec;
        for kk = 1 : length(MRSCont.overview.LCModel.sort_fit.(names{g}).(name))
          tempSubSpec(kk,:) = MRSCont.overview.LCModel.sort_fit.(names{g}).(name){1,kk}.fit;
          tempSubRes(kk,:) = MRSCont.overview.LCModel.sort_fit.(names{g}).(name){1,kk}.res;
          tempSubdata(kk,:) = MRSCont.overview.LCModel.sort_fit.(names{g}).(name){1,kk}.data;              
          if ~(strcmp(name, 'ref_ref') || strcmp(name, 'w_w'))
            tempSubBaseline(kk,:) = MRSCont.overview.LCModel.sort_fit.(names{g}).(name){1,kk}.baseline;            
             for f = 6 : length(fits)
                 if isfield(MRSCont.overview.LCModel.sort_fit.(names{g}).(name){1,kk},fits{f})
                    tempInidivMetab.(fits{f})(kk,:)= MRSCont.overview.LCModel.sort_fit.(names{g}).(name){1,kk}.(fits{f});
                 else
                    tempInidivMetab.(fits{f})(kk,:)=ones(1,length(MRSCont.overview.LCModel.sort_fit.(names{g}).(name){1,kk}.(fits{1}))) .* nan;
                 end
             end
          end
        end
        MRSCont.overview.LCModel.sort_fit.(names{g}).(['mean_' name]) = nanmean(real(tempSubSpec),1);
        MRSCont.overview.LCModel.sort_fit.(names{g}).(['sd_' name]) = nanstd(real(tempSubSpec),1);
        MRSCont.overview.LCModel.sort_fit.(names{g}).(['mean_res_' name]) = nanmean(real(tempSubRes),1);
        MRSCont.overview.LCModel.sort_fit.(names{g}).(['sd_res_' name]) = nanstd(real(tempSubRes),1);
        MRSCont.overview.LCModel.sort_fit.(names{g}).(['mean_data_' name]) = nanmean(real(tempSubdata),1);
        MRSCont.overview.LCModel.sort_fit.(names{g}).(['sd_data_' name]) = nanstd(real(tempSubdata),1);            
        MRSCont.overview.LCModel.sort_fit.(names{g}).(['mean_baseline_' name]) = nanmean(real(tempSubBaseline),1);
        MRSCont.overview.LCModel.sort_fit.(names{g}).(['sd_baseline_' name]) = nanstd(real(tempSubBaseline),1);
        for f = 6 : length(fits)
                MRSCont.overview.LCModel.sort_fit.(names{g}).(['mean_' fits{f} '_' name]) = nanmean(real(tempInidivMetab.(fits{f})),1);
                MRSCont.overview.LCModel.sort_fit.(names{g}).(['sd_' fits{f} '_' name]) = nanstd(real(tempInidivMetab.(fits{f})),1);
        end
        MRSCont.overview.LCModel.sort_fit.(names{g}).(['ppm_fit_' name]) = MRSCont.overview.LCModel.sort_fit.(names{g}).(name){1,1}.ppm;
        %Find the ppm of the maximum peak magnitude within the given range:
        ppmindex=find(MRSCont.overview.LCModel.sort_fit.(names{g}).(['mean_data_' name])(MRSCont.overview.LCModel.sort_fit.(names{g}).(['ppm_fit_' name])>1.9 & MRSCont.overview.LCModel.sort_fit.(names{g}).(['ppm_fit_' name])<2.1)==max(MRSCont.overview.LCModel.sort_fit.(names{g}).(['mean_data_' name])(MRSCont.overview.LCModel.sort_fit.(names{g}).(['ppm_fit_' name])>1.9 & MRSCont.overview.LCModel.sort_fit.(names{g}).(['ppm_fit_' name])<2.1)));
        ppmrange=MRSCont.overview.LCModel.sort_fit.(names{g}).(['ppm_fit_' name])(MRSCont.overview.LCModel.sort_fit.(names{g}).(['ppm_fit_' name])>1.9 & MRSCont.overview.LCModel.sort_fit.(names{g}).(['ppm_fit_' name])<2.1);
        ppmmax=ppmrange(ppmindex);
        refShift=(ppmmax-2.013);
         MRSCont.overview.LCModel.sort_fit.(names{g}).(['ppm_fit_' name]) = MRSCont.overview.LCModel.sort_fit.(names{g}).(['ppm_fit_' name]) - refShift;
        
end


end
