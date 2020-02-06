function Plot_two_hists(data1,data2,plotFolder,fileName,legendTitles,xLabel)
hold off;
h1 = histogram(data1,'Normalization','probability');
hold on;
h2 = histogram(data2,'Normalization','probability');
h2.BinWidth = h1.BinWidth;
xlabel(xLabel);
ylabel('Frequency');
set(gca,'fontsize',14);
grid on;
box on;
lgd = legend(legendTitles);
lgd.FontSize = 14;
filePath = fullfile(plotFolder,sprintf('%s.pdf',fileName));
print(gcf,filePath,'-dpdf','-bestfit');
end