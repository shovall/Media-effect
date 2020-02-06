function Plot_pie_chart_medias(ccle_metadata,plot_folder)
MIN_PERCENT = 7;

[grouping,labels] = findgroups(ccle_metadata.CellCultureMedia);

counts = zeros(length(labels),1);
for i=1:length(labels)
    counts(i) = length(find(grouping==i));
end

T = table(labels,counts,'VariableNames',{'Media','count'});
file_name = sprintf('%s\\ Media_counts.xlsx',plot_folder);
writetable(T, file_name);

selected = find(counts>MIN_PERCENT/100*sum(counts));
counts_selected = counts(selected);
labels_selected = labels(selected);

[counts_selected,order] = sort(counts_selected,'descend');
labels_selected = labels_selected(order);

counts_selected(end+1) = sum(counts)-sum(counts_selected);
labels_selected{end+1} = 'Others';

pie(counts_selected);
colormap([      0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330])
      
lgd = legend(labels_selected);
lgd.FontSize = 14;
lgd.Location = 'best';
set(gca,'fontsize',14);
filePath = fullfile(plot_folder,sprintf('%s.pdf','pie'));
print(gcf,filePath,'-dpdf','-bestfit');
end