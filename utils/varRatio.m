function [r, Rmeans] = varRatio(NVec, vertPotCell, SolCell)

var = 0;
idx = cumsum(NVec);
idx = [1 idx];
Rmeans = cell(length(NVec),1);
for i = 1:length(NVec)
    [varadd, Rmeans{i}] = var_SO3(SolCell(idx(i):idx(i+1)), vertPotCell(idx(i):idx(i+1)) );
    var = var + varadd;
end

r = var / var_SO3(SolCell, vertPotCell);