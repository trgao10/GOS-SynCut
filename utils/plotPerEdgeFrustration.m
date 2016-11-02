function [] = plotPerEdgeFrustration(G,perEdgeFrustMat,hsv)
%PLOTPEREDGEFRUSTRATION plot the edge colored proportionally to the
%                       per-edge frustration
%   
% Tingran Gao (trgao10@math.duke.edu)
% last modified: Oct 14, 2016
%

[rFrustMat,cFrustMat,vFrustMat] = find(triu(perEdgeFrustMat));
vFrustMat_rescaled = (vFrustMat-min(vFrustMat))/(max(vFrustMat)-min(vFrustMat));
cm_data = interp1(linspace(0,1,size(hsv,1)),hsv,vFrustMat_rescaled);
cm_data = hsv2rgb(cm_data);
for j=1:length(rFrustMat)
    staPtCoords = G.V(rFrustMat(j),:);
    endPtCoords = G.V(cFrustMat(j),:);
    line([staPtCoords(:,1);endPtCoords(:,1)],...
        [staPtCoords(:,2);endPtCoords(:,2)],'Color',cm_data(j,:));
    if j==1
        hold on
    end
end
axis equal
axis([0,G.numClusters,0,1]);

end

