function investigateHighErrorRates(rslt)
%INVESTIGATEHIGHERRORRATES
%
%
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Oct 16, 2016
%

n = size(rslt.G.adjMat,1)/rslt.params.numClusters;
tAdjMat = zeros(size(rslt.G.adjMat));
for j=1:rslt.params.numClusters
    blockIdx = ((j-1)*n+1):(j*n);
    tAdjMat(blockIdx,blockIdx) = rslt.G.adjMat(blockIdx,blockIdx);
end

figure;
plot(graph(tAdjMat), 'XData', rslt.G.V(:,1), 'YData', rslt.G.V(:,2), 'NodeLabel', {});
axis equal
axis([0,rslt.params.numClusters,0,1]);
hold on
for j=1:length(rslt.G.ccRowIdx)
    staPtCoords = rslt.G.V(rslt.G.ccRowIdx(j),:);
    endPtCoords = rslt.G.V(rslt.G.ccColIdx(j),:);
    line([staPtCoords(:,1);endPtCoords(:,1)],...
         [staPtCoords(:,2);endPtCoords(:,2)],...
         'Color', 'r', 'LineStyle', ':');
end
title(sprintf('Spectral Gap = %.2f, CCE/TTE = %d/%d',...
    rslt.G.specGap, length(rslt.G.ccRowIdx), sum(rslt.G.adjMat(:))/2),'Interpreter','latex');

[~, ~, GCL_W] = assembleGCL(rslt.G.adjMat, rslt.edgePotCell, rslt.params.d);

[GroundTruthPerEdgeFrustVec, GroundTruthPerEdgeFrustMat] =...
    getPerEdgeFrustration(rslt.G, GCL_W, rslt.params.d, rslt.vertPotCell);
[RelaxSolPerEdgeFrustVec, RelaxSolPerEdgeFrustMat] =...
    getPerEdgeFrustration(rslt.G, GCL_W, rslt.params.d, rslt.RelaxSolCell);
[combinedSolPerEdgeFrustVec, combinedSolPerEdgeFrustMat] =...
    getPerEdgeFrustration(rslt.G, GCL_W, rslt.params.d, rslt.combinedSolCell);
[CollageSolPerEdgeFrustVec, CollageSolPerEdgeFrustMat] =...
    getPerEdgeFrustration(rslt.G, GCL_W, rslt.params.d, rslt.CollageSolCell);

figure('Position',[50,100,800,600]);
subplot(2,2,1);
plotPerEdgeFrustration(rslt.G,GroundTruthPerEdgeFrustMat,rslt.params.hsv);
title(sprintf('GroundTruth total frustration = %.2f', sum(GroundTruthPerEdgeFrustVec)),'Interpreter','latex');
subplot(2,2,2);
plotPerEdgeFrustration(rslt.G,RelaxSolPerEdgeFrustMat,rslt.params.hsv);
title(sprintf('RelaxSol total frustration = %.2f', sum(RelaxSolPerEdgeFrustVec)),'Interpreter','latex');
subplot(2,2,3);
plotPerEdgeFrustration(rslt.G,combinedSolPerEdgeFrustMat,rslt.params.hsv);
title(sprintf('combinedSol total frustration = %.2f', sum(combinedSolPerEdgeFrustVec)),'Interpreter','latex');
subplot(2,2,4);
plotPerEdgeFrustration(rslt.G,CollageSolPerEdgeFrustMat,rslt.params.hsv);
title(sprintf('CollageSol total frustration = %.2f', sum(CollageSolPerEdgeFrustVec)),'Interpreter','latex');

figure('Position',[900,150,1000,400]);
for j=1:length(rslt.clusterLabel)
    scatter(rslt.G.V(rslt.clusterLabel{j},1),rslt.G.V(rslt.clusterLabel{j},2),...
        20,rslt.params.colorList{j},'filled');
    if j==1
        hold on
    end
end
axis equal
axis([0,rslt.params.numClusters,0,1]);
line([1,1],[0,1],'Color','g');
title(sprintf('Spectral Clustering errRate = %.2f', rslt.errRate),'Interpreter','latex');

end

