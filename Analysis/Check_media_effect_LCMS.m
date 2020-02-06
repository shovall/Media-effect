function res = Check_media_effect_LCMS(T_LCMS)

table_vars = fieldnames(T_LCMS);
RPMI_ids = find(contains(table_vars,'RPMI'));
DMEM_ids = find(contains(table_vars,'DMEM'));

fold_change = zeros(length(T_LCMS.compound),1)+NaN;
p_value = zeros(length(T_LCMS.compound),1)+NaN;
for i=1:length(T_LCMS.compound)
    vals1 = T_LCMS{i,RPMI_ids};
    vals2 = T_LCMS{i,DMEM_ids};
    if sum(vals1) ==0
        continue;
    end
    if sum(vals2) ==0
        continue;
    end

    [~,p] = ttest2(vals1,vals2);
    p_value(i) = p;
    fold_change(i) = mean(vals1)/mean(vals2);
end

res = struct;
res.metabolites = T_LCMS.compound;
res.p = p_value;
res.fold_change = fold_change;
res.pFDR = Calc_FDR(res.p);
end
