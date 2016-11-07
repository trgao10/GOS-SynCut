function [rslt] = SynCut_TwoWay(G, edgePotCell, params)
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

if numClusters > 2
    warning('For numClusters > 2, use the routine "SynCut" instead.');
    error('This routine only demos the details of SynCut using K=2.');
end

if debugFlag    
    try
        vertPotCell = params.vertPotCell;
    catch
        error('with debugFlag=true, must provide ground truth vertPotCell');
    end
    try
        colorList = getoptions(params, 'colorList', {'r','b','k','m'});
    catch
        error('with debugFlag=true, must provide ground truth colorList');
    end
    hsv = getoptions(params, 'hsv', rgb2hsv(winter));
    
    [GCL, GCL_Dvec, GCL_W] = assembleGCL(G.adjMat, edgePotCell, d);
    
    figure('Position',[30,550,560,420]);
    plot(graph(G.adjMat), 'XData', G.V(:,1), 'YData', G.V(:,2));
    axis equal
    axis([0,numClusters,0,1]);
    hold on
    for j=1:length(G.ccRowIdx)
        staPtCoords = G.V(G.ccRowIdx(j),:);
        endPtCoords = G.V(G.ccColIdx(j),:);
        line([staPtCoords(:,1);endPtCoords(:,1)],...
            [staPtCoords(:,2);endPtCoords(:,2)],'Color','r');
    end
    title(sprintf('Spectral Gap = %.2f, CCE/TTE = %d/%d',...
        G.specGap, length(G.ccRowIdx), sum(G.adjMat(:))/2),'Interpreter','latex');

    vertPotMat = cat(1, vertPotCell{:});
    GCL_unnormalized = diag(GCL_Dvec) - GCL_W;
    [GroundTruthPerEdgeFrustVec, GroundTruthPerEdgeFrustMat] =...
        getPerEdgeFrustration(G, GCL_W, d, vertPotCell);
    rescaled_vertPotMat = diag(sqrt(GCL_Dvec))*vertPotMat;
    fprintf('[GroundTruth] Rayleigh quotient (normalized GCL) = %f\n',...
        trace(rescaled_vertPotMat'*GCL*rescaled_vertPotMat));
    fprintf('[GroundTruth] Rayleigh quotient (unnormalized GCL) = %f\n',...
        trace(vertPotMat'*GCL_unnormalized*vertPotMat));
    fprintf('[GroundTruth] total frustration = %f\n',...
        sum(GroundTruthPerEdgeFrustVec));
end

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
    [RelaxSolCell, RelaxSolMat] = syncSpecRelax(wGCL, d, wGCL_Dvec);
    [RelaxSolPerEdgeFrustVec, RelaxSolPerEdgeFrustMat] =...
        getPerEdgeFrustFromEdgePot(G.adjMat, edgePotCell, RelaxSolCell);
    
    if debugFlag
        rescaled_RelaxSolMat = diag(sqrt(GCL_Dvec))*RelaxSolMat;
        fprintf('[RelaxSol] Rayleigh quotient (normalized GCL) = %f\n',...
            trace(rescaled_RelaxSolMat'*GCL*rescaled_RelaxSolMat));
        fprintf('[RelaxSol] Rayleigh quotient (unnormalized GCL) = %f\n',...
            trace(RelaxSolMat'*GCL_unnormalized*RelaxSolMat));
        fprintf('[RelaxSol] total frustration = %f\n',...
            sum(RelaxSolPerEdgeFrustVec));        
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
    
    if debugFlag
        %%%% for more consistent visualization effect, always set left
        %%%% cluster to be red an right cluster to be blue
        %%%% Note this will only work for K=2 (binary clustering)!
        NVec = G.NVec;
%         n = size(G.adjMat,1) / numClusters;
        if sum(clusterLabel{1} <= NVec(1)) < sum(clusterLabel{2} <= NVec(1))
            tmp = clusterLabel{1};
            clusterLabel{1} = clusterLabel{2};
            clusterLabel{2} = tmp;
            clear tmp
        end
        
        %%%% plot spectral clustering results
        if ~exist('specClusteringFigure', 'var')
            specClusteringFigure = figure('Position',[1200,550,1000,400]);
        else
            figure(specClusteringFigure);
        end
        subplot(1,2,1);
        for j=1:length(clusterLabel)
            scatter(G.V(clusterLabel{j},1),G.V(clusterLabel{j},2),...
                    20,colorList{j},'filled');
            if j==1
                hold on
            end
        end
        axis equal
        axis([0,numClusters,0,1]);
        line([1,1],[0,1],'Color','g');
        title(sprintf('Spectral Clustering in Iteration %d, Pass 1',iterCounter));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 3. Synchronization (local)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rescaled_relaxSol_cluster = cell(1,numClusters);
    for j=1:numClusters
        adjMat_cluster = G.adjMat(clusterLabel{j},clusterLabel{j});
        edgePotCell_cluster = edgePotCell(clusterLabel{j},clusterLabel{j});
        [GCL_cluster, GCL_Dvec_cluster, GCL_W_cluster] = assembleGCL(adjMat_cluster, edgePotCell_cluster, d);
        [~, rescaled_relaxSol_cluster{j}] = syncSpecRelax(GCL_cluster, d, GCL_Dvec_cluster);        
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
    
    if debugFlag
        rescaled_combinedSolMat = diag(sqrt(GCL_Dvec))*combinedSolMat;
        fprintf('[combinedSol] Rayleigh quotient (normalized Laplacian) = %f\n',...
            trace(rescaled_combinedSolMat'*GCL*rescaled_combinedSolMat));
        fprintf('[combinedSol] Rayleigh quotient (unnormalized Laplacian) = %f\n',...
            trace(combinedSolMat'*GCL_unnormalized*combinedSolMat));
        fprintf('[combinedSol] total frustration = %f\n',...
            sum(combinedSolPerEdgeFrustVec));
    end

        
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
    [ccSolCell, ccSolMat] = syncSpecRelax(ccGCL, d, ccGCL_Dvec);
    
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
    
    if debugFlag
        [ccRowIdx, ccColIdx] = find(triu(crossClusterAdjMat));
        inClusterTotalFrust = zeros(1,2);
        bwClusterTotalFrust = zeros(1,2);
        bwClusterTotalFrustVecs = zeros(2,length(ccRowIdx));
        for k=1:length(ccRowIdx)
            bwClusterTotalFrust(1) = bwClusterTotalFrust(1)+combinedSolPerEdgeFrustMat(ccRowIdx(k),ccColIdx(k));
            bwClusterTotalFrustVecs(1,k) = combinedSolPerEdgeFrustMat(ccRowIdx(k),ccColIdx(k));
            bwClusterTotalFrust(2) = bwClusterTotalFrust(2)+CollageSolPerEdgeFrustMat(ccRowIdx(k),ccColIdx(k));
            bwClusterTotalFrustVecs(2,k) = CollageSolPerEdgeFrustMat(ccRowIdx(k),ccColIdx(k));
        end
        inClusterTotalFrustCells = cell(2,length(clusterLabel));
        for j=1:length(clusterLabel)
            inClusterTotalFrust(1) = inClusterTotalFrust(1)+sum(sum(combinedSolPerEdgeFrustMat(clusterLabel{j},clusterLabel{j})))/2;
            [~,~,inClusterTotalFrustCells{1,j}] = find(triu(combinedSolPerEdgeFrustMat(clusterLabel{j},clusterLabel{j})));
            inClusterTotalFrust(2) = inClusterTotalFrust(2)+sum(sum(CollageSolPerEdgeFrustMat(clusterLabel{j},clusterLabel{j})))/2;
            [~,~,inClusterTotalFrustCells{2,j}] = find(triu(CollageSolPerEdgeFrustMat(clusterLabel{j},clusterLabel{j})));
        end
        
        rescaled_combinedSolMat = diag(sqrt(GCL_Dvec))*combinedSolMat;
        fprintf('[combinedSol] Rayleigh quotient (normalized Laplacian) = %f\n',...
            trace(rescaled_combinedSolMat'*GCL*rescaled_combinedSolMat));
        fprintf('[combinedSol] Rayleigh quotient (unnormalized Laplacian) = %f\n',...
            trace(combinedSolMat'*GCL_unnormalized*combinedSolMat));
        fprintf('[combinedSol] total frustration = %f\n',...
            sum(combinedSolPerEdgeFrustVec));
        fprintf('[combinedSol] (in+bw clusters) frustration = %f\n',...
            sum(inClusterTotalFrust(1)+bwClusterTotalFrust(1)));
        fprintf('[combinedSol] in-cluster frustration = %f\n',...
            sum(inClusterTotalFrust(1)));
        fprintf('[combinedSol] bw-cluster frustration = %f\n',...
            sum(bwClusterTotalFrust(1)));
        
        rescaled_CollageSolMat = diag(sqrt(GCL_Dvec))*CollageSolMat;
        fprintf('[CollageSol] Rayleigh quotient (normalized Laplacian) = %f\n',...
            trace(rescaled_CollageSolMat'*GCL*rescaled_CollageSolMat));
        fprintf('[CollageSol] Rayleigh quotient (unnormalized Laplacian) = %f\n',...
            trace(CollageSolMat'*GCL_unnormalized*CollageSolMat));
        fprintf('[CollageSol] total frustration = %f\n',...
            sum(CollageSolPerEdgeFrustVec));
        fprintf('[CollageSol] (in+bw clusters) frustration = %f\n',...
            sum(inClusterTotalFrust(2)+bwClusterTotalFrust(2)));
        fprintf('[CollageSol] in-cluster frustration = %f\n',...
            sum(inClusterTotalFrust(2)));
        fprintf('[CollageSol] bw-cluster frustration = %f\n',...
            sum(bwClusterTotalFrust(2)));
        
%         if ~exist('inCluserFigure', 'var')
%             inCluserFigure = figure('Name', 'in-cluster edge-wise frustrations should coincide');
%         else
%             figure(inCluserFigure);
%         end
%         subplot(1,2,1);
%         hist(cat(1,inClusterTotalFrustCells{1,:}));
%         title(sprintf('combinedSol, %.2f',sum(cat(1,inClusterTotalFrustCells{1,:}))));
%         subplot(1,2,2);
%         hist(cat(1,inClusterTotalFrustCells{2,:}));
%         title(sprintf('CollageSol, %.2f',sum(cat(1,inClusterTotalFrustCells{2,:}))));
        
%         if ~exist('bwClusterFigure', 'var')
%             bwClusterFigure = figure('Name', 'bw-cluster edge-wise frustrations');
%         else
%             figure(bwClusterFigure);
%         end
%         subplot(1,2,1);
%         hist(bwClusterTotalFrustVecs(1,:));
%         title(sprintf('combinedSol, %.2f',sum(bwClusterTotalFrustVecs(1,:))));
%         subplot(1,2,2);
%         hist(bwClusterTotalFrustVecs(2,:));
%         title(sprintf('CollageSol, %.2f',sum(bwClusterTotalFrustVecs(2,:))));
        
%         if ~exist('perEdgeFrustFigure', 'var')
%             perEdgeFrustFigure = figure('Position',[50,100,800,600]);
%         else
%             figure(perEdgeFrustFigure);
%         end
%         subplot(2,2,1);
%         plotPerEdgeFrustration_enhance_cc(G,GroundTruthPerEdgeFrustMat,hsv);
%         title(['GroundTruth total frustration = '...
%             num2str(sum(GroundTruthPerEdgeFrustVec))]);
%         subplot(2,2,2);
%         plotPerEdgeFrustration_enhance_cc(G,RelaxSolPerEdgeFrustMat,hsv);
%         title(['RelaxSol total frustration = '...
%             num2str(sum(RelaxSolPerEdgeFrustVec))]);
%         subplot(2,2,3);
%         plotPerEdgeFrustration_enhance_cc(G,combinedSolPerEdgeFrustMat,hsv);
%         title(['combinedSol total frustration = '...
%             num2str(sum(combinedSolPerEdgeFrustVec))]);
%         subplot(2,2,4);
%         plotPerEdgeFrustration_enhance_cc(G,CollageSolPerEdgeFrustMat,hsv);
%         title(['CollageSol total frustration = '...
%             num2str(sum(CollageSolPerEdgeFrustVec))]);
            
        if ~exist('perNodeFrobDistFigure', 'var')
            perNodeFrobDistFigure = figure('Position', [50,100,1600,400], 'Name','Frobenius dist b/t solution and groundtruth','NumberTitle','off');
        else
            figure(perNodeFrobDistFigure);
        end
        frobdist = @(x)cellfun(@(x,y)norm(x-y,'fro'), params.vertPotCell, x);
%         varRatio = @(x)(var(x(1:G.NVec(1)))*G.NVec(1) + var(x(G.NVec(1)+1:G.NVec(1)+G.NVec(2)))*G.NVec(2))/var(x)/length(x);
        varRatio = @(x)((var_SO3(params.vertPotCell(1:G.NVec(1)), x(1:G.NVec(1))) + var_SO3(params.vertPotCell(G.NVec(1)+1:end), x(G.NVec(1)+1:end)))/var_SO3(params.vertPotCell, x));
        RelaxSolFD = frobdist(RelaxSolCell);
        combinedSolFD = frobdist(combinedSolCell);
        CollageSolFD = frobdist(CollageSolCell);
        
        subplot(1, 3, 1);
        scatter(G.V( : , 1 ), G.V( : , 2 ), [], RelaxSolFD, 'filled');
        colormap(winter); colorbar;
        title( ['RelaxSol, in-class/total variation: ', num2str( varRatio(RelaxSolCell) )] );
        subplot(1, 3, 2);
        scatter(G.V( : , 1 ), G.V( : , 2 ), [], combinedSolFD, 'filled');
        colormap(winter); colorbar;
        title( ['combinedSol, in-class/total variation: ', num2str( varRatio(combinedSolCell) )] );
        subplot(1, 3, 3);
        scatter(G.V( : , 1 ), G.V( : , 2 ), [], CollageSolFD, 'filled');
        colormap(winter); colorbar;
        title( ['CollageSol, in-class/total variation: ', num2str( varRatio(CollageSolCell) )] );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 5. Spectral Clustering, Pass 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [clusterLabel,embedPts] =...
        specClusterWrapper(CollageSolPerEdgeFrustMat,numClusters,...
                           numKmeans,bandwidth,adjType);
    xi = evalXi(G, clusterLabel, d, CollageSolPerEdgeFrustMat);
    
    if debugFlag
        %%%% for more consistent visualization effect, always set left
        %%%% cluster to be red an right cluster to be blue
        %%%% Note this will only work for K=2 (binary clustering)!
        NVec = G.NVec;
%         n = size(G.adjMat,1) / numClusters;
        if sum(clusterLabel{1} <= NVec(1)) < sum(clusterLabel{2} <= NVec(1))
            tmp = clusterLabel{1};
            clusterLabel{1} = clusterLabel{2};
            clusterLabel{2} = tmp;
            clear tmp
        end
        
        %%%% plot spectral clustering results
        figure(specClusteringFigure);
        subplot(1,2,2);
        for j=1:length(clusterLabel)
            scatter(G.V(clusterLabel{j},1),G.V(clusterLabel{j},2),...
                    20,colorList{j},'filled');
            if j==1
                hold on
            end
        end
        axis equal
        axis([0,numClusters,0,1]);
        line([1,1],[0,1],'Color','g');
        title(sprintf('Spectral Clustering in Iteration %d, Pass 2',iterCounter));
        
        fprintf('iteration %d, xi = %4d\n', iterCounter, xi);
        fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++\n');
        pause();
    end
        
    [wRowIdx,wColIdx,wVals] = find(CollageSolPerEdgeFrustMat);
    wAdjMat = sparse(wRowIdx,wColIdx,exp(-wVals/mean(wVals)),size(wAdjMat,1),size(wAdjMat,2));
    
%     if (abs(xi) < tol) || (iterCounter >= maxIter) || ((xi_old < Inf) && (abs(xi-xi_old) < tol*xi_old))
%         break
%     end
    break;

end

rslt = struct('G', G, 'params', params, 'iterCounter', iterCounter);
[rslt.edgePotCell] = deal(edgePotCell);
[rslt.clusterLabel] = deal(clusterLabel);
[rslt.CollageSolCell] = deal(CollageSolCell);
rslt.embedPts = embedPts;

end
