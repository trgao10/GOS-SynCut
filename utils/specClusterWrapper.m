function [clusterLabel,embedPts] = specClusterWrapper(weightedAdj, numClusters, numKmeans, bandwidth, type)
%SPECTRALCLUSTERINGWRAPPER run spectral clustering multiple times and pick
%                          the best result
%   
% INPUTS:
%   weightedAdj ---------------------------- weighted adjacency matrix
%                                            (must be symmetric)
%   numClusters ---------------------------- number of clusters
%   numKmeans ------------------------------ number of runs
%   bandwidth ------------------------------ diffusion bandwidth paramerter
%                                            for the anisotropic kernel
%                                            (not useful if type == 'sim')
%   type ----------------------------------- indicating whether weightedAdj
%                                            contains similarity ('sim') or
%                                            dissimilarity ('dis') scores
%   debugFlag ------------------------------ print to screen
%
% OUTPUTS:
%   clusterLabel --------------------------- cell array of length
%                                            numClusters
%   
%
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Oct 25, 2016
% 

[rIdx,cIdx,vals] = find(triu(weightedAdj));
if strcmpi(type, 'dis')
    if isnumeric(bandwidth)
        simScore = exp(-vals/bandwidth);
    elseif ischar(bandwidth)
%         rowMins = zeros(size(weightedAdj,1),1);
%         for j=1:length(rowMins)
%             currRow = weightedAdj(j,:);
%             rowMins(j) = min(currRow(currRow > 0));
%         end        
%         simScore = exp(-vals/max(rowMins));
        simScore = exp(-vals/mean(vals));
    end
elseif strcmpi(type, 'sim')
    if (max(vals)-min(vals)) > 0
        simScore = (vals-min(vals))/(max(vals)-min(vals));
    else
        simScore = vals / max(vals);
    end
else
    error('unknown type');
end

W = sparse(rIdx,cIdx,simScore,size(weightedAdj,1),size(weightedAdj,2));
W = W+W';
% W = W + baseAdjWeights;
% W = min(W,baseAdjWeights);

Dvec = sum(W);
L = diag(1./sqrt(Dvec))*W*diag(1./sqrt(Dvec));
L = (L+L')/2;

[evecs, evals] = eig(L);
devals = 1-diag(evals);
[devals,evalidx] = sort(devals,'ascend');
% [~,evalidx] = sort(abs(diag(evals)),'descend');
evecs = evecs(:,evalidx);
% embedPts = sign(diag(1./sqrt(Dvec))*evecs(:,2)*diag(sqrt(devals(2))));
% embedPts = diag(1./sqrt(Dvec))*evecs(:,2:(numClusters+1));
embedPts = diag(1./sqrt(Dvec))*evecs(:,2:(numClusters+1))*sqrt(diag(1-devals(2:(numClusters+1))));
% embedPtsFull = diag(1./sqrt(Dvec))*evecs;
% embedPtsFull = diag(1./sqrt(Dvec))*evecs*sqrt(diag(1-devals));
% embedPts = diag(1./sqrt(Dvec))*evecs(:,2:(numClusters+1))*diag(sqrt(devals(2:(numClusters+1))));
% embedPts = diag(1./sqrt(Dvec))*evecs(:,1:(numClusters+1))*diag(sqrt(devals(1:(numClusters+1))));

% sumd = Inf;
% sumdList = zeros(numKmeans,1);
% if debugFlag
%     cback = 0;
% end
% for j=1:numKmeans
%     [t_idx, ~, t_sumd] = kmeans(embedPts,numClusters,'MaxIter',1000,'Replicates',5);
%     sumdList(j) = sum(t_sumd);
%     if sum(t_sumd) < sum(sumd)
%         cluster_idx = t_idx;
%         sumd = t_sumd;
%     end
%     
%     if debugFlag
%         for cc=1:cback
%             fprintf('\b');
%         end
%         cback = fprintf('kmeans clustering %4d/%d done.\n',j,numKmeans);
%     end
% end
cluster_idx = kmeans(embedPts,numClusters,'MaxIter',1000,'Replicates',numKmeans,'Display','off');
% if debugFlag
%     cluster_idx = kmeans(embedPts,numClusters,'MaxIter',1000,'Replicates',numKmeans,'Display','final');
% else
%     cluster_idx = kmeans(embedPts,numClusters,'MaxIter',1000,'Replicates',numKmeans,'Display','off');
% end
clusterLabel = cell(1,numClusters);
for j=1:numClusters
    clusterLabel{j} = find(cluster_idx == j);
end

end
