function T = Plot_scatter_one_vs_across_cells(res_one_line,res_across_lines,plot_folder)

[exists,locs] = ismember(lower(res_one_line.metabolites),lower(res_across_lines.metabolites));

mets = res_one_line.metabolites(exists);
X = log2(res_one_line.fold_change(exists));
Y = log2(res_across_lines.fold_change(locs(exists)));
T = table(2.^X,2.^Y,mets,'VariableNames',{'PANC1','across_lines','met'});

x_label = 'Fold change in PANC1 cells, RPMI vs DMEM, log_2';
y_label = 'Fold change across cell lines, log_2';

Plot_scatter(X,Y, x_label, y_label,plot_folder);
end
