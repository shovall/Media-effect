function  [res] = Check_CCLE_media_effect(ccle,ccle_metadata,media_info, field_str,plotFolder,toPlot)

[cellsRPMI,cellsDMEM] = Identify_cells_media(ccle_metadata);
[~, locs] = ismember(upper(cellsRPMI),upper(ccle.celllines));
locs1 = locs(locs~=0);

[~, locs] = ismember(upper(cellsDMEM),upper(ccle.celllines));
locs2 = locs(locs~=0);

num_mets = length(ccle.(field_str));
p = zeros(length(ccle.(field_str)),1);
p_control_cancer = zeros(length(ccle.(field_str)),1);
fold_change = zeros(length(ccle.(field_str)),1);
var_explained = zeros(length(ccle.(field_str)),1);

res = struct;
res.(field_str) = ccle.(field_str);
res.plot_links = cell(size(res.(field_str)));
legend_titles = {'RPMI','DMEM'};

cell_used = [upper(ccle.celllines(locs1)); upper(ccle.celllines(locs2))];
[~, locs] = ismember(cell_used,upper(ccle_metadata.NameWithTissueOrigins));
[grouping,labels] = findgroups(ccle_metadata.Classifications(locs));

mat = zeros(length(cell_used),length(labels));
for i=1:length(labels)
    cells_in_cancer = find(grouping==i);
    mat(cells_in_cancer,i) = 1;
end


for i=1:num_mets    
    data1 = ccle.data(locs1,i);
    data2 = ccle.data(locs2,i);
    
    current = res.(field_str){i};
    fileName = matlab.lang.makeValidName(current);
    if toPlot
        Plot_two_hists(log10(data1),log10(data2),plotFolder,fileName,legend_titles,'log 10 abundance');
    end
    res.plot_links{i} = sprintf('=HYPERLINK("%s","figure")',strcat(fileName,'.pdf')); 
    p(i) = ranksum(data1,data2);
    fold_change(i) = median(data1)/median(data2);
    
    y = [data1;data2];
    x = zeros(length(y),1);
    x(1:length(data1)) = 1;

    [~,p_control_cancer(i)] = partialcorr(x,y,mat,'type','Spearman');
end

pFDR = Calc_FDR(p);
[exists,locs] = ismember(res.metabolites,media_info.mets);
res.RPMI = zeros(size(res.metabolites))+NaN;
res.DMEM = zeros(size(res.metabolites))+NaN;
res.RPMI(exists) = media_info.RPMI(locs(exists));
res.DMEM(exists) = media_info.DMEM(locs(exists));
res.pFDR = pFDR;
res.p = p;
res.fold_change = fold_change;
res.p_control_cancer = p_control_cancer;
res.p_control_cancer_FDR = Calc_FDR(p_control_cancer);
end
