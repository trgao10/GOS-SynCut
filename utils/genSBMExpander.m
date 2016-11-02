function G = genSBMExpander(N,p,q,debugFlag)
%GENSBMEXPANDER Generate a graph representing a random graph drawn from a
%               Stochastic Block Model (SBM)
%   INPUT:
%       N --------------------- (vector) the numbers of vertices within
%                               each cluster
%       p --------------------- (scalar) probability of connectivity within
%                               each community
%       q --------------------- (scalar) probabiiltiy of connectivity
%                               between communities
%   OUTPUT:
%       G --------------------- struct containing all information about the
%                               expander generated
%           |--> V --------- Vertices of the graph (2-by-(n x numClusters) matrix)
%                            by default, cluster j lies in the unit square
%                            shifted j units to the right
%           |--> adjMat ---- sparse adjacency matrix of the graph
%           |--> specGap --- spectral gap of the generated graph
%           |--> ccRowIdx -- row (start) indices of cross-cluster links
%           |--> ccColIdx -- column (end) indices of cross-cluster links
%           |--> numClusters
%
%
%   
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Nov 2, 2016
%

if size(N,1) > size(N,2)
    N = N';
end

numClusters = length(N);
adjMat = tril(rand(sum(N)) < q);
posMarker = [0,cumsum(N)];
V = zeros(sum(N), 2);
for k = 1:numClusters;
    blockIdx = (posMarker(k)+1):posMarker(k+1);
    adjMat(blockIdx,blockIdx) = tril(rand(N(k)) < p);
    V(blockIdx,:) = rand([N(k),2]);
    V(blockIdx,1) = V(blockIdx,1)+(k-1);
end
adjMat = adjMat + adjMat';
adjMat = adjMat - diag(diag(adjMat));

if debugFlag
    fprintf('counting inter-cluster links...\n');
end
ccAdjMat = adjMat;
for k=1:numClusters
    blockIdx = (posMarker(k)+1):posMarker(k+1);
    ccAdjMat(blockIdx,blockIdx) = 0;
end
if debugFlag
    fprintf('Done.\n');
end

[ccRowIdx, ccColIdx] = find(triu(ccAdjMat,1));

if debugFlag
    fprintf('computing spectral gap.....\n');
end
Dvec = sum(adjMat);
GL = diag(1./sqrt(Dvec))*adjMat*diag(1./sqrt(Dvec));
evals = eig(GL);
evals = 1-evals;
evals = sort(evals);
specGap = evals(2);
if debugFlag
    fprintf('Done.\n');
end

G = struct('V', V, 'adjMat', adjMat, 'specGap', specGap,...
           'ccRowIdx', ccRowIdx, 'ccColIdx', ccColIdx,...
           'numClusters', numClusters);

if debugFlag
   imagesc(adjMat);
   axis equal
   axis square
   axis([1,sum(N),1,sum(N)]);
   title(sprintf('Spectral Gap = %.2f, CCE/TTE = %d/%d',...
        G.specGap, length(G.ccRowIdx), sum(G.adjMat(:))/2),'Interpreter','latex');
end

end

