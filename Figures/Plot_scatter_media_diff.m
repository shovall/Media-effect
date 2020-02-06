function T = Plot_scatter_media_diff(res,media_info,plot_folder)

locs1 = find(~isnan(media_info.RPMI));
locs2 = find(~isnan(media_info.DMEM));

locs_media_info = intersect(locs1,locs2);
mets  = media_info.mets(locs_media_info);

[exists,locs_ccle_res] = ismember(mets,res.metabolites);
locs_media_info = locs_media_info(exists);
mets = mets(exists);

locs_ccle_res = locs_ccle_res(exists);
X = log2(media_info.RPMI(locs_media_info)./media_info.DMEM(locs_media_info));
Y = log2(res.fold_change(locs_ccle_res));

T = table(2.^X,2.^Y,mets,'VariableNames',{'media','cells','met'});

x_label = 'Fold change in media composition, log_2';
y_label = 'Fold change measured in cell lines, log_2';

Plot_scatter(X,Y, x_label, y_label,plot_folder);
end
