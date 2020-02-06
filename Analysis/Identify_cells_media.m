function [cellsRPMI,cellsDMEM] = Identify_cells_media(ccle_metadata)

locsRPMI = find(strcmpi(ccle_metadata.CellCultureMedia,'RPMI 1640 + 10% FBS'));
cellsRPMI = ccle_metadata.NameWithTissueOrigins(locsRPMI);
locsDMEM = find(strcmpi(ccle_metadata.CellCultureMedia,'DMEM + 10% FBS'));
cellsDMEM = ccle_metadata.NameWithTissueOrigins(locsDMEM);
end
