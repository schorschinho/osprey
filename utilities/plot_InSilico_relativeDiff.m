function out = plot_InSilico_relativeDiff(MRSCont,step,par,metab,label)
    % Check that OspreySeg has been run before
    if ~isfield(MRSCont,'in_silico')
        error('No in sislico data included.')
    end
    % Check that OspreySeg has been run before
    if ~isfield(MRSCont,'fit')
        error('Trying to plot model results, but data has not been modeled yet. Run Ospreyfit first.')
    end

    [cb] = cbrewer('qual', 'Dark2', 12, 'pchip');
    temp = cb(3,:);
    cb(3,:) = cb(4,:);
    cb(4,:) = temp;

    out =[];
%     out = figure;
%     set(gcf, 'Color', 'w');
    
    for mm = 1: length(metab)
        for pp = 1 length(par)
            for kk = 1 : MRSCont.nDatasets(1)
                model_res(kk) = MRSCont.fit{kk}.Model{step}.parsOut.(par)(metab(mm));
            end
            relative_diff(:,mm) = (model_res - MRSCont.in_silico.par_full.(par)(metab(mm)))./MRSCont.in_silico.par_full.(par)(metab(mm)) * 100;
        end
        scatter((0.1 * ones(1,MRSCont.nDatasets(1))*mm) -0.01,relative_diff(:,mm),5, cb(1,:),'filled'), hold on
        
    end
    boxplot(relative_diff,'PlotStyle','compact', 'Colors', cb(1,:), 'Positions',[1:mm]/10,'Labels',label);
end