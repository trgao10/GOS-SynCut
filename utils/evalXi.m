function xi = evalXi(G, clusterLabel, d, perEdgeFrustration)
%EVALXI evaluate total frustration for given partition and vertex potential
%   
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Oct 14, 2016
%

numClusters = length(clusterLabel);
adjMatMask = G.adjMat;
perClusterVolume = zeros(1,numClusters);
degVec = sum(G.adjMat);
%%%%%%%%%%% this is completely wrong, should be redone!!!
for j=1:numClusters
    adjMatMask(clusterLabel{j},clusterLabel{j}) = 0;
    perClusterVolume(j) = sum(degVec(clusterLabel{j}));
%     perClusterVolume(j) = sum(sum(G.adjMat(clusterLabel{j},clusterLabel{j})))/2;
end
adjMatMask = G.adjMat - adjMatMask;
adjMatMask = double(adjMatMask > 0);

xi = full((sum(sum(perEdgeFrustration.*adjMatMask)))/2 * ...
      sum(1./perClusterVolume)/(2*d*sum(sum(degVec))));

end

