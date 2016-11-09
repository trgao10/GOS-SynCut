function [rslt] = SynCut(G, edgePotCell, params)
%SYNCUT_TWOWAY_HDM
%   INPUTS:
%     G ---------------------------- graph
%     edgePotCell ------------------ edge potential in cell array format
%     params ----------------------- struct holding all parameters
%
%   OUTPUTS:
%     rslt ------------------------- struct holding all outputs
%
%
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Oct 26, 2016
%

debugFlag = getoptions(params, 'debugFlag', false);
d = getoptions(params, 'd', Inf);
numClusters = getoptions(params, 'numClusters', 2);
tol = getoptions(params, 'tol', 1e-8);
maxIter = getoptions(params, 'maxIter', 10);
numKmeans = getoptions(params, 'numKmeans', 200);
bandwidth = getoptions(params, 'bandwidth', 1);
adjType = getoptions(params, 'adjType', 'dis');
syncRoutine = getoptions(params, 'syncRoutine', @syncSpecRelax);

xi = Inf;
iterCounter = 0;
wAdjMat = G.adjMat;
while true
    iterCounter = iterCounter+1;
    if debugFlag
        fprintf('++++++++ Iteration %d +++++++++++++++++++++++++++++\n',...
                iterCounter);
    end
    
    xi_old = xi;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 1. Synchronization (global)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [wGCL, wGCL_Dvec, wGCL_W] = assembleGCL(wAdjMat, edgePotCell, d);
    [RelaxSolCell, RelaxSolMat] = syncRoutine(wGCL, d, wGCL_Dvec);
%     [RelaxSolCell, RelaxSolMat] = syncSpecRelax(wGCL, d, wGCL_Dvec);
    [RelaxSolPerEdgeFrustVec, RelaxSolPerEdgeFrustMat] =...
        getPerEdgeFrustFromEdgePot(G.adjMat, edgePotCell, RelaxSolCell);
    if iterCounter == 1
        origRelaxSolCell = RelaxSolCell;
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 2. Spectral Clustering, Pass 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while true
        clusterLabel =...
            specClusterWrapper(RelaxSolPerEdgeFrustMat,numClusters,...
            numKmeans,bandwidth,adjType);
        %%%%% we expect spectral clustering to partition the graph into two
        %%%%% connected components
        connTest = zeros(size(clusterLabel));
        numConnComp = zeros(size(clusterLabel));
        for j=1:length(connTest)
            connTest(j) = ~isempty(find(sum(G.adjMat(clusterLabel{j},clusterLabel{j})) == 0, 1));
            numConnComp(j) = length(unique(conncomp(graph(G.adjMat(clusterLabel{j},clusterLabel{j})))));
        end
        if (sum(connTest) == 0) && (all(numConnComp == 1))
            break
        end
        error('error: spectral clustering does not produce connected components');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 3. Synchronization (local)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rescaled_relaxSol_cluster = cell(1,numClusters);
    for j=1:numClusters
        adjMat_cluster = G.adjMat(clusterLabel{j},clusterLabel{j});
        edgePotCell_cluster = edgePotCell(clusterLabel{j},clusterLabel{j});
        [GCL_cluster, GCL_Dvec_cluster, GCL_W_cluster] = assembleGCL(adjMat_cluster, edgePotCell_cluster, d);
        [~, rescaled_relaxSol_cluster{j}] = syncRoutine(GCL_cluster, d, GCL_Dvec_cluster);        
%         [~, rescaled_relaxSol_cluster{j}] = syncSpecRelax(GCL_cluster, d, GCL_Dvec_cluster);        
    end
    
    combinedSolCell = cell(1,size(G.adjMat,1));
    combinedSolMat = zeros(size(G.adjMat,1)*d,d);
    for j=1:numClusters
        for k=1:length(clusterLabel{j})
            tmpIdx = ((clusterLabel{j}(k)-1)*d+1):(clusterLabel{j}(k)*d);
            tmpIdxInCluster = ((k-1)*d+1):(k*d);
            tmpBlock = rescaled_relaxSol_cluster{j}(tmpIdxInCluster,:);
            [U,S,V] = svd(tmpBlock);
            combinedSolCell{clusterLabel{j}(k)} = U*V';
            combinedSolMat(tmpIdx,:) = combinedSolCell{clusterLabel{j}(k)};
        end
    end
    
    [combinedSolPerEdgeFrustVec, combinedSolPerEdgeFrustMat] =...
            getPerEdgeFrustFromEdgePot(G.adjMat, edgePotCell, combinedSolCell);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 4. Collage
    % [TODO] Try implementing minimum instead of average?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CollageSolCell = combinedSolCell;
    CollageSolMat = combinedSolMat;
    
    %%%%%%%%% identify cross-cluster edges
    crossClusterAdjMat = G.adjMat;
    for j=1:numClusters
        crossClusterAdjMat(clusterLabel{j},clusterLabel{j}) = 0;
    end
    % [ccRowIdx, ccColIdx] = find(triu(crossClusterAdjMat));
    
    %%%%%%%%% form cross-partition edge potential
    ccEdgePotCell = cell(numClusters, numClusters);
    ccAdjMat = zeros(numClusters,numClusters);
    for j=1:numClusters
        for k=(j+1):numClusters
            localCCAdj = crossClusterAdjMat(clusterLabel{j},clusterLabel{k});
            [ccRowIdx, ccColIdx] = find(localCCAdj);
            if ~isempty(ccRowIdx)
                ccAdjMat(j,k) = 1;
                ccAdjMat(k,j) = 1;
            end
            localCCRowIdx = clusterLabel{j}(ccRowIdx);
            localCCColIdx = clusterLabel{k}(ccColIdx);
            for jj=1:length(localCCRowIdx)
                if isempty(ccEdgePotCell{j,k})
                    ccEdgePotCell{j,k} = combinedSolPerEdgeFrustMat(localCCRowIdx(jj),localCCColIdx(jj))*combinedSolCell{localCCRowIdx(jj)}'*edgePotCell{localCCRowIdx(jj),localCCColIdx(jj)}*combinedSolCell{localCCColIdx(jj)};
                else
                    ccEdgePotCell{j,k} = ccEdgePotCell{j,k}+combinedSolPerEdgeFrustMat(localCCRowIdx(jj),localCCColIdx(jj))*combinedSolCell{localCCRowIdx(jj)}'*edgePotCell{localCCRowIdx(jj),localCCColIdx(jj)}*combinedSolCell{localCCColIdx(jj)};
                end
            end
            ccEdgePotCell{k,j} = ccEdgePotCell{j,k}';
        end
    end
    
    %%%%%%%%% synchronize the ccGraph
    [ccGCL, ccGCL_Dvec, ccGCL_W] = assembleGCL(ccAdjMat, ccEdgePotCell, d);
    [ccSolCell, ccSolMat] = syncRoutine(ccGCL, d, ccGCL_Dvec);
%     [ccSolCell, ccSolMat] = syncSpecRelax(ccGCL, d, ccGCL_Dvec);
    
    %%%%%%%%% perform the collage
    for j=1:numClusters
        for k=1:length(clusterLabel{j})
            CollageSolCell{clusterLabel{j}(k)} = CollageSolCell{clusterLabel{j}(k)}*ccSolCell{j};
            tmpIdx = ((clusterLabel{j}(k)-1)*d+1):(clusterLabel{j}(k)*d);
            CollageSolMat(tmpIdx,:) = CollageSolMat(tmpIdx,:)*ccSolCell{j};
        end
    end
        
    [CollageSolPerEdgeFrustVec, CollageSolPerEdgeFrustMat] =...
        getPerEdgeFrustFromEdgePot(G.adjMat, edgePotCell, CollageSolCell);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 5. Spectral Clustering, Pass 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [clusterLabel,embedPts] =...
        specClusterWrapper(CollageSolPerEdgeFrustMat,numClusters,...
                           numKmeans,bandwidth,adjType);
    xi = evalXi(G, clusterLabel, d, CollageSolPerEdgeFrustMat);
        
    [wRowIdx,wColIdx,wVals] = find(CollageSolPerEdgeFrustMat);
    wAdjMat = sparse(wRowIdx,wColIdx,exp(-wVals/mean(wVals)),size(wAdjMat,1),size(wAdjMat,2));
    
    if (abs(xi) < tol) || (iterCounter >= maxIter) || ((xi_old < Inf) && (abs(xi-xi_old) < tol*xi_old))
        break
    end
end

%%%%%% very ad-hoc error counts and error rates computation for K=2
rslt = struct('G', G, 'params', params, 'iterCounter', iterCounter);
[rslt.edgePotCell] = deal(edgePotCell);
[rslt.clusterLabel] = deal(clusterLabel);
[rslt.CollageSolCell] = deal(CollageSolCell);
[rslt.RelaxSolCell] = deal(origRelaxSolCell);
rslt.embedPts = embedPts;

if exist('combinedSolCell', 'var')
    [rslt.combinedSolCell] = deal(combinedSolCell);
end

end
