function [MRSCont] =  osp_import_Tarquin_spectra(MRSCont)
%% Get sequence type
if MRSCont.flags.isUnEdited
    name = 'off_A';
end

PathName = uipickfiles('FilterSpec',MRSCont.outputFolder,'REFilter', '\','NumFiles',1,'Prompt','Select the folder containing all Tarquin fit .txt- and .csv- files.');
PathName = PathName{1};
d1 = dir(fullfile(PathName));
files = length(d1);

j = 1;
for i = 1 : files
    if ~strcmp (d1(i,1).name, '.')
                if ~strcmp (d1(i,1).name, '..')
                    if ~strcmp (d1(i,1).name(1), '.')
                        if strcmp(d1(i,1).name(end-4), 't')
                            if ~exist('metab')
                                data = osp_import_Tarquin_data_fit_baseline(fullfile(PathName,d1(i,1).name));
                                temp_sz(1,j)= length(data(:,1));
                                temp_ppm{j} = data(:,1);
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
                        if strcmp(d1(i,1).name(end-4), 't')
                            if ~exist('metab')
                                data = osp_import_Tarquin_fit_ouput(fullfile(PathName,d1(i,1).name));
                                dataNames = cellstr(osp_import_Tarquin_metabolite_names(fullfile(PathName,d1(i,1).name)));
                                data = table2array(data(:,1:length(dataNames)));
                                 waterAmp = osp_importfile_waterAmpl_Tarquin(strrep(fullfile(PathName,d1(i,1).name),'_fit.txt','.csv'));
                                 phase = osp_importfile_ph0ph1_Tarquin(strrep(fullfile(PathName,d1(i,1).name),'_fit.txt','.csv'));
                                 scale = waterAmp;
                                if length(data(:,1)) < max_point
                                    ppmRangeData        = temp_ppm{max_ind};
                                    ppmRangeDataToInt       = data(:,1);
                                    ppmIsInDataRange    = (ppmRangeDataToInt < ppmRangeData(1)) & (ppmRangeDataToInt > ppmRangeData(end));
                                    if sum(ppmIsInDataRange) == 0
                                        ppmIsInDataRange    = (ppmRangeDataToInt > ppmRangeData(1)) & (ppmRangeDataToInt < ppmRangeData(end));
                                    end
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.data      = interp1(ppmRangeDataToInt(ppmIsInDataRange), data(ppmIsInDataRange,2), ppmRangeData, 'pchip', 'extrap');
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.fit      = interp1(ppmRangeDataToInt(ppmIsInDataRange), data(ppmIsInDataRange,3), ppmRangeData, 'pchip', 'extrap') + interp1(ppmRangeDataToInt(ppmIsInDataRange), data(ppmIsInDataRange,4), ppmRangeData, 'pchip', 'extrap');
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.baseline      = interp1(ppmRangeDataToInt(ppmIsInDataRange), data(ppmIsInDataRange,4), ppmRangeData, 'pchip', 'extrap');
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.ppm = ppmRangeData';
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.res = MRSCont.overview.Tarquin.all_models.(name){1,j}.data-MRSCont.overview.Tarquin.all_models.(name){1,j}.fit;
                                    for n = 5 : length(dataNames)
                                        MRSCont.overview.Tarquin.all_models.(name){1,j}.(['fit' dataNames{n}])  = interp1(ppmRangeDataToInt(ppmIsInDataRange), data(ppmIsInDataRange,n), ppmRangeData, 'pchip', 'extrap');
                                    end
                                else
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.ppm = data(:,1)';
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.data = data(:,2);
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.fit = data(:,3) + data(:,4);
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.baseline = data(:,4);
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.res = MRSCont.overview.Tarquin.all_models.(name){1,j}.data-MRSCont.overview.Tarquin.all_models.(name){1,j}.fit;
                                    for n = 5 : length(dataNames)
                                        MRSCont.overview.Tarquin.all_models.(name){1,j}.(['fit' dataNames{n}])  = data(:,n);
                                    end
                                end
                                if MRSCont.opts.fit.fitMM == 1
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.fittMM  = 0;
                                    idx_1  = find(strcmp(dataNames,'MM09'));
                                    if isempty(idx_1)
                                        idx_1  = find(strcmp(dataNames,'MMexp'));
                                    end
                                    for f = idx_1 : length(dataNames)
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.fittMM  = MRSCont.overview.Tarquin.all_models.(name){1,j}.fittMM + MRSCont.overview.Tarquin.all_models.(name){1,j}.(['fit' dataNames{f}]);
                                    end
                                end
                                idx_1  = find(strcmp(dataNames,'NAA'));
                                idx_2  = find(strcmp(dataNames,'NAAG'));
                                MRSCont.overview.Tarquin.all_models.(name){1,j}.fittNAA  = MRSCont.overview.Tarquin.all_models.(name){1,j}.(['fit' dataNames{idx_1}]) + MRSCont.overview.Tarquin.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                idx_Cr  = find(strcmp(dataNames,'Cr'));
                                idx_PCr  = find(strcmp(dataNames,'PCr'));
                                MRSCont.overview.Tarquin.all_models.(name){1,j}.fittCr  = MRSCont.overview.Tarquin.all_models.(name){1,j}.(['fit' dataNames{idx_Cr}]) + MRSCont.overview.Tarquin.all_models.(name){1,j}.(['fit' dataNames{idx_PCr}]);
                                idx_1  = find(strcmp(dataNames,'GPC'));
                                idx_2  = find(strcmp(dataNames,'PCh'));
                                MRSCont.overview.Tarquin.all_models.(name){1,j}.fittCho  = MRSCont.overview.Tarquin.all_models.(name){1,j}.(['fit' dataNames{idx_1}]) + MRSCont.overview.Tarquin.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                idx_1  = find(strcmp(dataNames,'Glu'));
                                idx_2  = find(strcmp(dataNames,'Gln'));
                                MRSCont.overview.Tarquin.all_models.(name){1,j}.fitGlx  = MRSCont.overview.Tarquin.all_models.(name){1,j}.(['fit' dataNames{idx_1}]) + MRSCont.overview.Tarquin.all_models.(name){1,j}.(['fit' dataNames{idx_2}]);
                                MRSCont.overview.Tarquin.all_models.(name){1,j}.data = MRSCont.overview.Tarquin.all_models.(name){1,j}.data/scale;
                                MRSCont.overview.Tarquin.all_models.(name){1,j}.fit = MRSCont.overview.Tarquin.all_models.(name){1,j}.fit/scale;
                                MRSCont.overview.Tarquin.all_models.(name){1,j}.baseline = MRSCont.overview.Tarquin.all_models.(name){1,j}.baseline/scale;
                                MRSCont.overview.Tarquin.all_models.(name){1,j}.res = MRSCont.overview.Tarquin.all_models.(name){1,j}.res/scale;
                                modelNames = fields(MRSCont.overview.Tarquin.all_models.(name){1,j});
                                for n = 6 : length(fields(MRSCont.overview.Tarquin.all_models.(name){1,j}))
                                    MRSCont.overview.Tarquin.all_models.(name){1,j}.(modelNames{n})  = MRSCont.overview.Tarquin.all_models.(name){1,j}.(modelNames{n})/scale;
                                end
                                %Find the ppm of the maximum peak magnitude within the given range:
                                ppmindex=find(MRSCont.overview.Tarquin.all_models.(name){1,j}.data(MRSCont.overview.Tarquin.all_models.(name){1,j}.ppm>1.9 & MRSCont.overview.Tarquin.all_models.(name){1,j}.ppm<2.1)==max(MRSCont.overview.Tarquin.all_models.(name){1,j}.data(MRSCont.overview.Tarquin.all_models.(name){1,j}.ppm>1.9 & MRSCont.overview.Tarquin.all_models.(name){1,j}.ppm<2.1)));
                                ppmrange=MRSCont.overview.Tarquin.all_models.(name){1,j}.ppm(MRSCont.overview.Tarquin.all_models.(name){1,j}.ppm>1.9 & MRSCont.overview.Tarquin.all_models.(name){1,j}.ppm<2.1);
                                ppmmax=ppmrange(ppmindex);
                                refShift=(ppmmax-2.013);
                                MRSCont.overview.Tarquin.all_models.(name){1,j}.ppm = MRSCont.overview.Tarquin.all_models.(name){1,j}.ppm - refShift;
                                MRSCont.overview.Tarquin.all_models.(name){1,j}.ph0 = phase(1);
                                MRSCont.overview.Tarquin.all_models.(name){1,j}.ph1 = phase(2);
                                j = j + 1;
                                
                            end
                        end
                    end
                end
    end
end

for kk = 1 : MRSCont.nDatasets
    data = MRSCont.overview.Tarquin.all_models.(name){1,kk};
    modelNames = fields(data);
    Cr_height = max(data.data(data.ppm>2.9 & data.ppm<3.1));
    for n = 1 : length(modelNames)
        if ~strcmp('ppm',modelNames{n}) && ~strcmp('ph0',modelNames{n}) && ~strcmp('ph1',modelNames{n})
            MRSCont.overview.Tarquin.all_models.(name){1,kk}.(modelNames{n})  = MRSCont.overview.Tarquin.all_models.(name){1,kk}.(modelNames{n})/Cr_height;
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
    MRSCont.overview.Tarquin.sort_fit.(['g_' num2str(g)]).(name) = MRSCont.overview.Tarquin.all_models.(name)(1,MRSCont.overview.groups == g);
end
MRSCont.overview.Tarquin.sort_fit.GMean.(name) = MRSCont.overview.Tarquin.all_models.(name)(1,MRSCont.overview.groups > 0);


%Fits
names = fields(MRSCont.overview.Tarquin.sort_fit);
for g = 1 : length(names)
        tempSubSpec = zeros(length(MRSCont.overview.Tarquin.sort_fit.(names{g}).(name){1}),length(MRSCont.overview.Tarquin.sort_fit.(names{g}).(name){1}.ppm));
        tempSubBaseline = tempSubSpec;
        tempSubRes = tempSubSpec;
        tempSubdata = tempSubSpec;
        for kk = 1 : length(MRSCont.overview.Tarquin.sort_fit.(names{g}).(name))
          tempSubSpec(kk,:) = MRSCont.overview.Tarquin.sort_fit.(names{g}).(name){1,kk}.fit;
          tempSubRes(kk,:) = MRSCont.overview.Tarquin.sort_fit.(names{g}).(name){1,kk}.res;
          tempSubdata(kk,:) = MRSCont.overview.Tarquin.sort_fit.(names{g}).(name){1,kk}.data;              
          if ~(strcmp(name, 'ref_ref') || strcmp(name, 'w_w'))
            tempSubBaseline(kk,:) = MRSCont.overview.Tarquin.sort_fit.(names{g}).(name){1,kk}.baseline;
            fits = fields(MRSCont.overview.Tarquin.sort_fit.(names{g}).(name){1,kk});
             for f = 6 : length(fits)
                    tempInidivMetab.(fits{f})(kk,:)= MRSCont.overview.Tarquin.sort_fit.(names{g}).(name){1,kk}.(fits{f});
             end
          end
        end
        MRSCont.overview.Tarquin.sort_fit.(names{g}).(['mean_' name]) = nanmean(real(tempSubSpec),1);
        MRSCont.overview.Tarquin.sort_fit.(names{g}).(['sd_' name]) = nanstd(real(tempSubSpec),1);
        MRSCont.overview.Tarquin.sort_fit.(names{g}).(['mean_res_' name]) = nanmean(real(tempSubRes),1);
        MRSCont.overview.Tarquin.sort_fit.(names{g}).(['sd_res_' name]) = nanstd(real(tempSubRes),1);
        MRSCont.overview.Tarquin.sort_fit.(names{g}).(['mean_data_' name]) = nanmean(real(tempSubdata),1);
        MRSCont.overview.Tarquin.sort_fit.(names{g}).(['sd_data_' name]) = nanstd(real(tempSubdata),1);            
        MRSCont.overview.Tarquin.sort_fit.(names{g}).(['mean_baseline_' name]) = nanmean(real(tempSubBaseline),1);
        MRSCont.overview.Tarquin.sort_fit.(names{g}).(['sd_baseline_' name]) = nanstd(real(tempSubBaseline),1);
        for f = 6 : length(fits)
                MRSCont.overview.Tarquin.sort_fit.(names{g}).(['mean_' fits{f} '_' name]) = nanmean(real(tempInidivMetab.(fits{f})),1);
                MRSCont.overview.Tarquin.sort_fit.(names{g}).(['sd_' fits{f} '_' name]) = nanstd(real(tempInidivMetab.(fits{f})),1);
        end
        MRSCont.overview.Tarquin.sort_fit.(names{g}).(['ppm_fit_' name]) = MRSCont.overview.Tarquin.sort_fit.(names{g}).(name){1,1}.ppm;
        %Find the ppm of the maximum peak magnitude within the given range:
        ppmindex=find(MRSCont.overview.Tarquin.sort_fit.(names{g}).(['mean_data_' name])(MRSCont.overview.Tarquin.sort_fit.(names{g}).(['ppm_fit_' name])>1.9 & MRSCont.overview.Tarquin.sort_fit.(names{g}).(['ppm_fit_' name])<2.1)==max(MRSCont.overview.Tarquin.sort_fit.(names{g}).(['mean_data_' name])(MRSCont.overview.Tarquin.sort_fit.(names{g}).(['ppm_fit_' name])>1.9 & MRSCont.overview.Tarquin.sort_fit.(names{g}).(['ppm_fit_' name])<2.1)));
        ppmrange=MRSCont.overview.Tarquin.sort_fit.(names{g}).(['ppm_fit_' name])(MRSCont.overview.Tarquin.sort_fit.(names{g}).(['ppm_fit_' name])>1.9 & MRSCont.overview.Tarquin.sort_fit.(names{g}).(['ppm_fit_' name])<2.1);
        ppmmax=ppmrange(ppmindex);
        refShift=(ppmmax-2.013);
         MRSCont.overview.Tarquin.sort_fit.(names{g}).(['ppm_fit_' name]) = MRSCont.overview.Tarquin.sort_fit.(names{g}).(['ppm_fit_' name]) - refShift;
        
end

end
