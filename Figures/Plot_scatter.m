function Plot_scatter(X,Y,x_label,y_label,plot_folder)
scatter(X,Y)
box on;
grid on;
xlabel(x_label);
ylabel(y_label);
set(gca,'fontsize',14);

lm = fitlm(X,Y);
hold on;
xPlotRegress = linspace(min(X),max(X),10);
yPlotRegress = feval(lm,xPlotRegress);
plot(xPlotRegress,yPlotRegress,'LineWidth',1.5,'Color',[0 0.4470 0.7410] ,'LineStyle','--','HandleVisibility','off');

[r,p] = corr(X,Y);
fprintf('r = %.2f, p = %e \n',r,p);

fileName = sprintf('scatter r =%.2f, p = %e',r,p);
filePath = fullfile(plot_folder,sprintf('%s.pdf',fileName));
print(gcf,filePath,'-dpdf','-bestfit');
hold off;
end