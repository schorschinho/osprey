
ys = zeros(MRSCont.nDatasets,length(MRSCont.fit.basisSet.name));
for kk = 1:MRSCont.nDatasets
ys(kk,:) = MRSCont.quantify.tCr{kk}';
end

h=boxplot(ys,MRSCont.fit.basisSet.name,'factorgap',[10],'labelverbosity','minor','notch','off');
% h =  boxplot(metab,group);
set(h,'LineWidth',2);


for i = 1 : length(MRSCont.fit.basisSet.name)
set(h(:,i),'Color','k');
patch(get(h(5,i),'XData'),get(h(5,i),'YData'),'k','FaceAlpha',.2);
end;

h = findobj(gcf,'tag','Outliers');
for iH = 1:length(h)
    h(iH).MarkerEdgeColor = 'k';
end

ylabel('[metabolite] / tCr');

%%
metab_name = repmat({'.' '..' '...'},1,6);
groups = [repmat({'GABA'},1,3), repmat({'Gln'},1,3),repmat({'mI'},1,3),repmat({'Asp'},1,3),repmat({'GSH'},1,3),repmat({'GPC'},1,3)];

h=boxplot(ys,{groups,metab_name},'factorgap',[10 0.5],'labelverbosity','minor','notch','off');
% h =  boxplot(metab,group);
set(h,'LineWidth',2);


for i = 1 : 3 : 16
set(h(:,i),'Color','k');
patch(get(h(5,i),'XData'),get(h(5,i),'YData'),'k','FaceAlpha',.2);
set(h(:,i+1),'Color','g');
patch(get(h(5,i+1),'XData'),get(h(5,i+1),'YData'),'g','EdgeColor','g','FaceAlpha',.2);
set(h(:,i+2),'Color','r');
patch(get(h(5,i+2),'XData'),get(h(5,i+2),'YData'),'r','EdgeColor','r','FaceAlpha',.2);
end;

h = findobj(gcf,'tag','Outliers');
for iH = 1:3 : length(h)-2
    h(iH).MarkerEdgeColor = 'r';
    h(iH+1).MarkerEdgeColor = 'g';
    h(iH+2).MarkerEdgeColor = 'k';
end

ylabel('[metabolite] / tCr');

