function Plot_volcano(res, plotFolder,labels,file_name)
Q_THRESH_LABLELS = 30;
Y = -log10(res.pFDR);
X = log2(res.fold_change);

scatter(X,Y)
set(gca,'fontsize',14);
grid on;
box on;

xlabel('log_2 fold change, RPMI vs DMEM');
ylabel('-log_1_0 (q-value)');

line(xlim, [-log10(0.05),-log10(0.05)],'LineStyle','--','Color','black','HandleVisibility','off','LineWidth',1.5);

locs = find(Y>Q_THRESH_LABLELS);
texts_arr = labels(locs);
for i=1:length(locs)
    label = sprintf(' %s \\rightarrow ',texts_arr{i});
    text(X(locs(i)),Y(locs(i)),label,'HorizontalAlignment','right');
end

filePath = fullfile(plotFolder,sprintf('volcano %s.pdf',file_name));
print(gcf,filePath,'-dpdf','-bestfit');
hold off;
end

