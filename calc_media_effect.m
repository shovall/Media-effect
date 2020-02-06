%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Config
% Folder to plot analysis results
plot_folder = '.\Media';

% CCLE metabolomics metadata, available at:
% Li, H., Ning, S., Ghandi, M. et al.
% The landscape of cancer cell line metabolism. 
% Nat Med 25, 850–860 (2019) doi:10.1038/s41591-019-0404-8
% Supplementary Table 1: cell line annotations
ccle_metabolomics_metadata_path = '.\Data\ccle_metabolomics_annotations.xlsx';
% CCLE metabolomics data file, available at:
% "https://portals.broadinstitute.org/ccle"
ccle_metabolomics_path = '.\Data\CCLE_metabolomics_20190502.csv';

% Media information file path, includes RPMI and DMEM
media_info_path = '.\Data\media-info.xlsx';

% Data file locations, PANC1 cells cultured in RPMI and DMEM
files_path = {'.\Data\PANC1_NEG.xlsx','.\Data\PANC1_POS.xlsx'};
pcv_counts_path = '.\Data\PANC1_pcv_counts.xlsx';

% Bool value to plot histograms showing media effect in each metabolite
plot_histograms = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% Analyzes metabolomics profiling from the CCLE
%%%%%%%%%
% Loads CCLE metabolomics data
ccle_metadata = Load_ccle_metabolomics_metadata(ccle_metabolomics_metadata_path);
ccle_metabolomics = Load_CCLE_metabolomics(ccle_metabolomics_path);
% Loads media information
media_info = readtable(media_info_path);

% Plots pie chart showing common media types
Plot_pie_chart_medias(ccle_metadata,plot_folder);

% Checks for media effect for each metabolite 
[res] = Check_CCLE_media_effect(ccle_metabolomics,ccle_metadata,media_info,'metabolites',plot_folder,plot_histograms);
T_media_effect_CCLE = table(res.metabolites,res.p,res.pFDR,res.fold_change,...
res.p_control_cancer,res.p_control_cancer_FDR,res.RPMI, res.DMEM, res.plot_links,...
'VariableNames',{'met','p','p_FDR','median_fold_change','p_controlled','p_controlled_FDR',...
'RPMI','DMEM','figure'});

writetable(T_media_effect_CCLE,sprintf('%s//CCLE_Metabolomics_media_effect.xlsx',plot_folder));
Plot_volcano(res, plot_folder,res.metabolites,'CCLE metabolomics media effect');

% Calculates fold change difference correlation between media to
% intracellular levels
T_scatter = Plot_scatter_media_diff(res,media_info,plot_folder);
writetable(T_scatter,sprintf('%s//Scatter_media_vs_cells.xlsx',plot_folder));

%%%%%%%%%
% Analyzes PANC1 LC-MS metabolic profiling results
%%%%%%%%%
% Loads data
T_LCMS = Load_LCMS(files_path,pcv_counts_path);
writetable(T_LCMS,sprintf('%s//PANC1_LCMS.xlsx',plot_folder));

% Checks for media effect for each metabolite 
res_LCMS = Check_media_effect_LCMS(T_LCMS);
T_res_LCMS = table(res_LCMS.metabolites,res_LCMS.p,res_LCMS.pFDR,res_LCMS.fold_change,'VariableNames',{'met','p','p_FDR','fold_change'});
writetable(T_res_LCMS,sprintf('%s//PANC1_media_effect.xlsx',plot_folder));
Plot_volcano(res_LCMS, plot_folder,res_LCMS.metabolites,'PANC1');

% Calculates fold change difference correlation between PANC1 cells
% cultured in RPMI and DMEM to median fold change across cell lines from
% the CCLE
T_scatter_LCMS = Plot_scatter_one_vs_across_cells(res_LCMS,res,plot_folder);
writetable(T_scatter_LCMS,sprintf('%s//Scatter_PANC1_vs_across_lines.xlsx',plot_folder));
