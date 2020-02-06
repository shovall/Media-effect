function ccle_metadata = Load_ccle_metabolomics_metadata(ccle_metadata_path)

ccle_metadata = readtable(ccle_metadata_path);
ccle_metadata.TissueOrigin = ccle_metadata.NameWithTissueOrigins;

for i=1:length(ccle_metadata.TissueOrigin)
    current = ccle_metadata.TissueOrigin{i};
    name = ccle_metadata.Name{i};
    ccle_metadata.TissueOrigin{i} = current(length(name)+2:end);
end
end
