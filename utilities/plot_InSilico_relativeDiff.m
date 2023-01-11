function out = plot_InSilico_relativeDiff(MRSCont,step,par,metab)
    % Check that OspreySeg has been run before
    if ~isfield(MRSCont,'in_silico')
        error('No in sislico data included.')
    end
    % Check that OspreySeg has been run before
    if ~isfield(MRSCont,'fit')
        error('Trying to plot model results, but data has not been modeled yet. Run Ospreyfit first.')
    end
    
    if nargin < 4
        metab = [];
    end

    [cb] = cbrewer('qual', 'Dark2', 12, 'pchip');
    temp = cb(3,:);
    cb(3,:) = cb(4,:);
    cb(4,:) = temp;

    out =[];

    if strcmp(par,'metAmpl') || strcmp(par,'lorentzLB') || strcmp(par,'freqShift')
        for mm = 1: length(metab)
            metabs_in_model  = MRSCont.fit{1}.BasisSets.names;
            metabs_in_model(~MRSCont.fit{1}.BasisSets.includeInFit(step,:)) = [];
            idx = find(strcmp(metabs_in_model,metab{mm}));
            for kk = 1 : MRSCont.nDatasets(1)
                model_res(kk) = MRSCont.fit{kk}.Model{step}.parsOut.(par)(idx);
            end
            metabs_in_sim  = MRSCont.in_silico.par_full.metabolite_names;
            idx_sim = find(strcmp(metabs_in_sim,metab{mm}));
            relative_diff(:,mm) = (model_res' - MRSCont.in_silico.par_full.(par)(:,idx_sim))./MRSCont.in_silico.par_full.(par)(:,idx_sim) * 100;
            scatter((0.1 * ones(1,MRSCont.nDatasets(1))*mm) -0.05,relative_diff(:,mm),5, cb(1,:),'filled'), hold on
            
        end
        boxplot(relative_diff,'PlotStyle','compact', 'Colors', cb(1,:), 'Positions',[1:mm]/10,'Labels',metab);
    else if ~strcmp(par,'baseAmpl')

            for kk = 1 : MRSCont.nDatasets(1)
                model_res(kk) = MRSCont.fit{kk}.Model{step}.parsOut.(par);
            end
            relative_diff(:) = (model_res' - MRSCont.in_silico.par_full.(par))./MRSCont.in_silico.par_full.(par) * 100;
            relative_diff(isinf(relative_diff(:)))= model_res(isinf(relative_diff(:)));
            scatter((0.1 * ones(1,MRSCont.nDatasets(1))) -0.05,relative_diff(:),5, cb(1,:),'filled'), hold on
                
            boxplot(relative_diff,'PlotStyle','compact', 'Colors', cb(1,:), 'Positions',1/10,'Labels',{par});
    else
         for mm = 1: length(metab)
            for kk = 1 : MRSCont.nDatasets(1)
                model_res(kk) = MRSCont.fit{kk}.Model{step}.parsOut.(par);
            end
            relative_diff(:) = (model_res' - MRSCont.in_silico.par_full.(par))./MRSCont.in_silico.par_full.(par) * 100;
            relative_diff(isinf(relative_diff(:)))= model_res(isinf(relative_diff(:)));
            scatter((0.1 * ones(1,MRSCont.nDatasets(1))) -0.05,relative_diff(:),5, cb(1,:),'filled'), hold on
                
            boxplot(relative_diff,'PlotStyle','compact', 'Colors', cb(1,:), 'Positions',1/10,'Labels',{par});
         end
        end
    end
end