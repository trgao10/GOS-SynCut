function G = genWeakExpander(n, numClusters, numLinks, degBounds, ccType)
%GENWEAKEXPANDER generate a weak expander, consisting of a few dense
%                clusters and relatively sparse inter-cluster links
%   INPUTS：
%   n -------------- number of vertices within each cluster
%   numClusters ---- number of clusters
%   numLinks ------- number of links between each pair of clusters
%   degBounds ------ vertex degree [lower,upper] bounds within each cluster
%   ccType --------- in-cluster connectivity type ['nn'|'unif']
%
%   OUTPUTS：
%   G -------------- struct containing all information about the generated
%                    weak expander
%   |--> V --------- Vertices of the graph (2-by-(n x numClusters) matrix)
%                    by default, cluster j lies in the unit square
%                    shifted j units to the right
%   |--> adjMat ---- sparse adjacency matrix of the graph
%   |--> specGap --- spectral gap of the generated graph
%   |--> ccRowIdx -- row (start) indices of cross-cluster links
%   |--> ccColIdx -- column (end) indices of cross-cluster links
%   |--> numClusters
%   |--> numLinks
%   |--> degBounds
%
%
%   Tingran Gao (trgao10@math.duke.edu)
%   last modified: Oct 17, 2016
%

V = zeros(n*numClusters, 2);
adjMat = zeros(n*numClusters);
for j=1:numClusters
    blockIdx = ((j-1)*n+1):(j*n);
    V(blockIdx,:) = rand([n,2]);
    V(blockIdx,1) = V(blockIdx,1)+(j-1);    
    degVec = randi(degBounds, [n,1]);
    if mod(sum(degVec),2) == 1
        degVec(find(degVec>1, 1)) = degVec(find(degVec>1, 1))-1;
    end
%     while ~checkErdoesGallai(degVec)
%         disp('re-generate degVec due to violation of Erdoes-Gallai');
%         degVec = randi(degBounds, [n,1]);
%     end
    blockAdj = zeros(n);
    while true
        if strcmpi(ccType, 'nn')
            distMat = squareform(pdist(V(blockIdx,:))) + diag(Inf(n,1));
            [~,idx] = sort(distMat,2);
            for jj=1:n
                blockAdj(jj, idx(jj,1:degVec(jj))) = 1;
            end
            blockAdj = min(blockAdj,blockAdj');
            isolatedIdx = find(sum(blockAdj) < degVec');
            if ~isempty(isolatedIdx)
                for jj=1:length(isolatedIdx)
                    blockAdj(isolatedIdx(jj),...
                        idx(isolatedIdx(jj),1:degVec(isolatedIdx(jj)))) = 1;
                end
            end
            blockAdj = max(blockAdj,blockAdj');
        elseif strcmpi(ccType, 'unif')
            origDegVec = degVec;
            while any(degVec)
                remIdx = find(degVec ~= 0);
                [~,sortRemIdx] = sort(degVec(remIdx), 'descend');
                remIdx = remIdx(sortRemIdx);
                if degVec(remIdx(1)) > (length(remIdx)-1)
                    disp('error: single degree too large; re-generate');
                    degVec = origDegVec;
                    continue
                end
                newLinkEnds = remIdx(datasample(2:length(remIdx),degVec(remIdx(1)),'Replace',false));
                for k=1:length(newLinkEnds)
                    blockAdj(remIdx(1),newLinkEnds(k)) = 1;
                    blockAdj(newLinkEnds(k),remIdx(1)) = 1;
                    degVec([remIdx(1),newLinkEnds(k)]) = degVec([remIdx(1),newLinkEnds(k)])-1;
                end
            end
        else
            error('unknown ccType');
        end
        
        %%%% if a component is not connected, re-generate
        connBins = conncomp(graph(blockAdj));
        if length(unique(connBins)) == 1
            break
        end
    end
    adjMat(blockIdx,blockIdx) = blockAdj;
end

%%%% shrink cluster to its centroid for better visualization quality
for j=1:numClusters
    blockIdx = ((j-1)*n+1):(j*n);
    Vblock = V(blockIdx,:);
    VCenter = mean(Vblock,1);
    Vblock = (Vblock-repmat(VCenter,size(Vblock,1),1))*0.8 + ...
        repmat(VCenter,size(Vblock,1),1);
    V(blockIdx,:) = Vblock;
end

ccRowIdx = zeros(1,numLinks*nchoosek(numClusters,2));
ccColIdx = zeros(1,numLinks*nchoosek(numClusters,2));
count = 1;
for j=1:numClusters
    blockIdxRow = ((j-1)*n+1):(j*n);
    for k=(j+1):numClusters
        blockIdxCol = ((k-1)*n+1):(k*n);
        nzRowIdx = randi([1,n], [numLinks,1]);
        nzColIdx = randi([1,n], [numLinks,1]);
        ccRowIdx(((count-1)*numLinks+1):(count*numLinks)) = nzRowIdx+(j-1)*n;
        ccColIdx(((count-1)*numLinks+1):(count*numLinks)) = nzColIdx+(k-1)*n;
        blockAdj = sparse(nzRowIdx, nzColIdx, ones(1,numLinks), n, n, numLinks);
        adjMat(blockIdxRow,blockIdxCol) = full(blockAdj);
        count = count+1;
    end
end

adjMat = double(max(adjMat, adjMat') > 0);
if max(adjMat(:)) > 1
    error('adjMat should only contain 0 or 1 as entries');
end

Dvec = sum(adjMat);
GL = diag(1./sqrt(Dvec))*adjMat*diag(1./sqrt(Dvec));
evals = eig(GL);
evals = 1-evals;
evals = sort(evals);
specGap = evals(2);

G = struct('V', V, 'adjMat', adjMat, 'specGap', specGap,...
           'ccRowIdx', ccRowIdx, 'ccColIdx', ccColIdx,...
           'numClusters', numClusters, 'numLinks', numLinks,...
           'degBounds', degBounds);

end

function flag = checkErdoesGallai(degSeq)

sortDegSeq = sort(degSeq, 'descend');
cumsumDegSeq = cumsum(sortDegSeq);
targetSeq = (1:length(sortDegSeq)).*(0:(length(sortDegSeq)-1));
for k=1:length(targetSeq)
    targetSeq(k) = targetSeq(k) + sum(min(k,sortDegSeq((k+1):end)));
end
flag = all(cumsumDegSeq <= targetSeq');

end

