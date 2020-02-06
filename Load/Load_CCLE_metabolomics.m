function ccleMetabolomics = Load_CCLE_metabolomics(ccleMetabolomicsPath)

ccleMetabolomics = struct;
[r,d] = xlsread(ccleMetabolomicsPath);
ccleMetabolomics.celllines = d(2:end,1);
ccleMetabolomics.metabolites = d(1,3:end)';
ccleMetabolomics.data = 10.^r;
end

