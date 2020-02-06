function T = Load_LCMS(files_path,pcv_counts_path)
INT_FIELD = 'peakAreaTop';
pcv_counts = readtable(pcv_counts_path);
pcv_fields = fieldnames(pcv_counts);
sample_vars = pcv_fields(1:end-3);
res = cell(2,1);
for i=1:2
    cur_res = struct;
    cur_T = readtable(files_path{i});
    cur_res.compound= unique(cur_T.compoundName);

    for sample = sample_vars'
        cur_res.(sample{1}) = zeros(length(cur_res.compound),1)+NaN;
    end
    
    for j=1:length(cur_res.compound)
        comp = cur_res.compound{j};
        locs = find(strcmp(cur_T.compoundName,comp));

        for loc=locs'
            cur_res.(cur_T.sample{loc})(j) = cur_T.(INT_FIELD)(loc);
        end
    end
    res{i} = cur_res;
end

T1 = struct2table(res{1});
T2 = struct2table(res{2});

[all_compounds,~,~] = union(T1.compound,T2.compound);

[~,ia] = ismember(all_compounds,T1.compound);
[~,ib] = ismember(all_compounds,T2.compound);

selected1 = [];
selected2 = [];
for i=1:length(ia)
    if ia(i) == 0
        selected2(end+1) = ib(i);
        continue;
    elseif ib(i) == 0
        selected1(end+1) = ia(i);
        continue;
    end

    sample_signal1 = T1{ia(i),sample_vars};
    sample_signal2 = T2{ib(i),sample_vars};
    
    val1 = median(sample_signal1);
    val2 = median(sample_signal2);
    if val1>val2
        selected1(end+1) = ia(i);
    else
        selected2(end+1) = ib(i);
    end
end
T = [T1(selected1,:);T2(selected2,:)];
for field = pcv_fields(1:end-3)'
    pcv = pcv_counts.(field{:});
    T.(field{:}) = T.(field{:})./pcv*median(pcv_counts{1,:});
end
end