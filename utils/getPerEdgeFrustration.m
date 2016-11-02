function [perEdgeFrustVec, perEdgeFrustMat] = getPerEdgeFrustration(G, GCL_W, d, vertPot)
%GETPEREDGEFRUSTRATION collection frustration on each edge of graph G
%
% INPUTS:
%   G -------------------------- graph
%   GCL_W ---------------------- prescribed edge potential
%   d -------------------------- dimension of each block of GCL_W
%   vertPot -------------------- vertex potential of which the per edge
%                                frustration is to be collected, stored in
%                                cell array format
%
% OUTPUTS:
%   perEdgeFrustVec ------------ vector of length nnz(perEdgeFrustMat)/2
%                                storing per edge frustration
%   perEdgeFrustMat ------------ matrix of the same size as G.adjMat
%                                encoding edge frustration at the
%                                corresponding edge in G.adjMat
%   
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Oct 14, 2016
%

[rIdx, cIdx] = find(triu(G.adjMat));
perEdgeFrustVec = zeros(length(rIdx),1);
for j=1:length(rIdx)
    rIdx_j = rIdx(j);
    cIdx_j = cIdx(j);
    rowIdx = ((rIdx_j-1)*d+1):(rIdx_j*d);
    colIdx = ((cIdx_j-1)*d+1):(cIdx_j*d);
    perEdgeFrustVec(j) = ...
      norm(vertPot{rIdx(j)}-GCL_W(rowIdx,colIdx)*vertPot{cIdx(j)},'fro')^2;
end

perEdgeFrustMat = sparse(rIdx,cIdx,perEdgeFrustVec,size(G.adjMat,1),size(G.adjMat,2));
perEdgeFrustMat = perEdgeFrustMat + perEdgeFrustMat';

end

