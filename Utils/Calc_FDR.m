function [ adjustedP ] = Calc_FDR(p)
origSize = size(p);
p_as_1D_vector = reshape(p,[(size(p,2)*size(p,1)),1]);

adjustedFDR = mafdr(p_as_1D_vector,'BHFDR',true);
adjustedP = reshape(adjustedFDR,origSize);
end

