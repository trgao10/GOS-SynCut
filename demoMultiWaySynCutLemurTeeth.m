%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup paths and parameters
data_path = '/home/trgao10/Work/MATLAB/DATA/HDM/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
sample_path = [data_path 'samples/'];
result_path = '/home/trgao10/Work/MATLAB/ArchivedResults/HDM/cPdist/';
dist_path = [result_path 'cPDistMatrix.mat'];
maps_path = [result_path 'cPMapsMatrix.mat'];
GroupLevel = 'Genus';
GroupNames = {'Alouatta','Ateles','Brachyteles','Callicebus','Saimiri'};
DisplayLayout = [5,10];
visFlag = true;

BNN = 7; %%% perfect spot for separating three clusters
colorsList = [228,26,28;0,126,204;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191]/255;
markerList = {'+','o','x','^','d'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% parameters for SynCut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numClusters = 3;
d = 3;
maxIter = 10;
tol = 1e-8;
numKmeans = 200;
bandwidth = 'auto'; %%%% ['auto' | some number]
adjType = 'dis';
colorList = {'r','b','k','m'};
hsv = rgb2hsv(winter);
close(gcf);

%% visualization options
options = struct('sample_path', sample_path,...
                 'DisplayLayout', DisplayLayout,...
                 'DisplayOrient', 'Horizontal',...
                 'boundary', 'on', 'names', 'off',...
                 'linkCamera', 'on', 'Shading', 'Smooth');

%% load mapsMatrix
load(maps_path);
mapsMatrix = cPMapsMatrix;

%% load taxa codes
taxa_file = [data_path 'hdm_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;

%% parse GroupNames
[~,ClTable,~] = xlsread(spreadsheet_path);
Names = {};
NamesByGroup = cell(1,length(GroupNames));
TaxaByGroup = cell(1,length(GroupNames));
for j=1:length(GroupNames)
    NamesJ = ClTable(strcmpi(ClTable(1:end,strcmpi(ClTable(1,:),GroupLevel)),GroupNames{j}),1);
    Names = [Names,NamesJ{:}];
    NamesByGroup{j} = NamesJ;
    TaxaByGroup{j} = cellfun(@(x) find(strcmpi(taxa_code, x)), NamesJ);
end
GroupSize = length(Names);

%% load meshes
TAXAinds = zeros(GroupSize,1);
meshList = cell(1,GroupSize);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,Names{j}));
    load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']);
    meshList{j} = G;
end
Names = taxa_code(TAXAinds); % match upper/lower cases

LocalByGroup = TaxaByGroup;
for j=1:length(TaxaByGroup)
    for k=1:length(TaxaByGroup{j})
        TaxaByGroup{j}(k) = find(TAXAinds == TaxaByGroup{j}(k));
    end
end

foliIdx = [TaxaByGroup{1};TaxaByGroup{3}];
frugIdx = [TaxaByGroup{2};TaxaByGroup{4}];
inscIdx = TaxaByGroup{5};
DietLabels = zeros(GroupSize,1);
DietLabels(foliIdx) = 1;
DietLabels(frugIdx) = 2;
DietLabels(inscIdx) = 3;

%% visualize un-aligned meshes
if visFlag
    drawMeshList(meshList, options);
    set(gcf, 'Name', 'un-aligned');
end

%% build weighted graph from cP distance matrix
clear G
%%%%% build KNN graph
load(dist_path);
BaseDistMatrix = cPDistMatrix(TAXAinds,TAXAinds);
BaseDistMatrix = BaseDistMatrix-diag(diag(BaseDistMatrix));

if BNN > (GroupSize-1)
    warning('BNN larger than allowed by sample size');
    BNN = min(BNN,GroupSize-1);
end

%%% only connect BNN-nearest-neighbors
[sDists,rowNNs] = sort(BaseDistMatrix,2);
sDists = sDists(:,2:min(1+BNN,GroupSize));
rowNNs = rowNNs(:,2:min(1+BNN,GroupSize));
% BaseWeights = sparse(repmat((1:GroupSize)',1,BNN),rowNNs,exp(-sDists.^2/BaseEps));
BaseWeights = sparse(repmat((1:GroupSize)',1,BNN),rowNNs,sDists);
BaseWeights = max(BaseWeights, BaseWeights');

tD = BaseDistMatrix;
tD = tD+diag(Inf(GroupSize,1));
epsilon = mean(min(tD, [] ,2));
W = exp(-tD.^2/epsilon^2);
W = W.*double(BaseWeights>0);
D = sum(W,2);
L = diag(1./sqrt(D))*W*diag(1./sqrt(D));
L = (L+L')/2;

[evecs, evals] = eig(L);
devals = 1-diag(evals);
[devals,evalidx] = sort(devals,'ascend');
evecs = evecs(:,evalidx);
Ydm = diag(1./sqrt(D))*evecs(:,2:(numClusters+5))*sqrt(diag(1-devals(2:(numClusters+5))));

embedPts = Ydm(:,1:3);
three_cluster_idx_dm = kmeans(embedPts,numClusters,'MaxIter',1000,'Replicates',numKmeans,'Display','off');
three_clusterLabel_dm = cell(1,numClusters);
for j=1:numClusters
    three_clusterLabel_dm{j} = find(three_cluster_idx_dm == j);
end
five_cluster_idx_dm = kmeans(embedPts,5,'MaxIter',1000,'Replicates',numKmeans,'Display','off');
five_clusterLabel_dm = cell(1,5);
for j=1:5
    five_clusterLabel_dm{j} = find(five_cluster_idx_dm == j);
end

planeFigure = figure();
subplot(1,2,1);
scatter3(Ydm(foliIdx,1),Ydm(foliIdx,2),Ydm(foliIdx,3),20,'r','filled');
hold on
scatter3(Ydm(frugIdx,1),Ydm(frugIdx,2),Ydm(frugIdx,3),20,'b','filled');
scatter3(Ydm(inscIdx,1),Ydm(inscIdx,2),Ydm(inscIdx,3),20,'g','filled');
title('Diffusion Maps');

%%%%% build kNN graph
G.adjMat = full(BaseWeights);
G.adjMat = double(G.adjMat > 0);
% G.adjMat = G.adjMat.*W;

adjMask = G.adjMat;
G.V = rand(GroupSize,2);
for j=1:length(TaxaByGroup)
    blockIdx = TaxaByGroup{j};
    adjMask(blockIdx,blockIdx) = 0;
    G.V(blockIdx,1) = G.V(blockIdx,1)+(j-1);
    %%%% shrink cluster to its centroid for better visual quality
    Vblock = G.V(blockIdx,:);
    VCenter = mean(Vblock,1);
    Vblock = (Vblock-repmat(VCenter,size(Vblock,1),1))*0.8 + ...
        repmat(VCenter,size(Vblock,1),1);
    G.V(blockIdx,:) = Vblock;
end

[G.ccRowIdx,G.ccColIdx] = find(triu(adjMask));

Dvec = sum(G.adjMat);
GL = diag(1./sqrt(Dvec))*G.adjMat*diag(1./sqrt(Dvec));
evals = eig(GL);
evals = 1-evals;
evals = sort(evals);
G.specGap = evals(2);
G.numClusters = numClusters;

%% load rigid motions output from cPDist
R = cell(GroupSize,GroupSize);
cback = 0;
for j=1:GroupSize
    for k = 1:GroupSize
        [~,R{k,j},~] = MapToDist(meshList{j}.V,meshList{k}.V,...
            mapsMatrix{TAXAinds(j),TAXAinds(k)},meshList{j}.Aux.VertArea);
    end
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf(['%4d/' num2str(GroupSize) ' done.\n'],j);
end

[wGCL, wGCL_Dvec, wGCL_W] = assembleGCL(G.adjMat, R, d);
if ~issymmetric(wGCL)
    warning('wGCL not symmetric!');
    wGCL = triu(wGCL,1);
    wGCL = wGCL+wGCL'+eye(size(wGCL));
end
if ~issymmetric(wGCL_W)
    warning('wGCL_W not symmetric!');
    wGCL_W = triu(wGCL_W,1);
    wGCL_W = wGCL_W+wGCL_W';
end
[RelaxSolCell, RelaxSolMat] = syncSpecRelax(wGCL, d, wGCL_Dvec);

%%%%% transpose each element to obtain the correct solution
SCSolCell = cellfun(@(x) x', RelaxSolCell, 'UniformOutput', false);

%%%%% make edge potential
RCell = cell(size(R));
for j=1:size(RCell,1)
    for k=(j+1):size(RCell,2)
        RCell{j,k} = R{j,k};
        RCell{k,j} = RCell{j,k}';
    end
end

vertPotCell = R(1,:);
GTTotalFrust = 0;
RSTotalFrust = 0;
SCTotalFrust = 0;
for j=1:size(RCell,1)
    for k=(j+1):size(RCell,2)
        GTTotalFrust = GTTotalFrust+G.adjMat(j,k)*norm(vertPotCell{j}'*vertPotCell{k}-RCell{j,k}, 'fro')^2;
        RSTotalFrust = RSTotalFrust+G.adjMat(j,k)*norm(RelaxSolCell{j}*RelaxSolCell{k}'-RCell{j,k}, 'fro')^2;
        SCTotalFrust = SCTotalFrust+G.adjMat(j,k)*norm(SCSolCell{j}'*SCSolCell{k}-RCell{j,k}, 'fro')^2;
    end
end

fprintf('[GroundTruth] GTTotal = %f\n', GTTotalFrust);
fprintf('[RelaxSol] RSTotal = %f\n', RSTotalFrust);
fprintf('[SCSol] SCTotal = %f\n', SCTotalFrust);

[GTSolPerEdgeFrustVec, GTSolPerEdgeFrustMat] =...
    getPerEdgeFrustFromEdgePot(G.adjMat, RCell, cellfun(@(x) x', vertPotCell, 'UniformOutput', false));
[RSSolPerEdgeFrustVec, RSSolPerEdgeFrustMat] =...
    getPerEdgeFrustFromEdgePot(G.adjMat, RCell, RelaxSolCell);
[SCSolPerEdgeFrustVec, SCSolPerEdgeFrustMat] =...
    getPerEdgeFrustFromEdgePot(G.adjMat, RCell, cellfun(@(x) x', SCSolCell, 'UniformOutput', false));
fprintf('[GroundTruth] verify GTTotal = %f\n', sum(GTSolPerEdgeFrustVec));
fprintf('[RelaxSol] verify RSTotal = %f\n', sum(RSSolPerEdgeFrustVec));
fprintf('[SCSol] verify SCTotal = %f\n', sum(SCSolPerEdgeFrustVec));

%% visualize aligned meshes
visCell = SCSolCell;
if visFlag
    for j=1:GroupSize
        meshList{j}.V = visCell{j} * meshList{j}.V;
        if det(visCell{j}) < 0
            meshList{j}.F = meshList{j}.F([1 3 2], :);
        end
    end
    drawMeshList(meshList, options);
    set(gcf, 'Name', 'aligned');
end

%% run multi-way SynCut
debugFlag = false;

params = struct('debugFlag', debugFlag,...
                'd', d,...
                'numClusters', numClusters,...
                'tol', tol,...
                'maxIter', maxIter,...
                'numKmeans', numKmeans,...
                'bandwidth', bandwidth,...
                'adjType', adjType,...
                'hsv', hsv);
[params.vertPotCell] = deal(vertPotCell);

rslt = SynCut(G, RCell, params);

five_cluster_idx = kmeans(rslt.embedPts,5,'MaxIter',1000,'Replicates',numKmeans,'Display','off');
rslt.five_clusterLabel = cell(1,5);
for j=1:5
    rslt.five_clusterLabel{j} = find(five_cluster_idx == j);
end


%% generate statistics for clustering results
figure(planeFigure);
subplot(1,2,2);
scatter3(rslt.embedPts(foliIdx,1),rslt.embedPts(foliIdx,2),rslt.embedPts(foliIdx,3),20,'r','filled');
hold on
scatter3(rslt.embedPts(frugIdx,1),rslt.embedPts(frugIdx,2),rslt.embedPts(frugIdx,3),20,'b','filled');
scatter3(rslt.embedPts(inscIdx,1),rslt.embedPts(inscIdx,2),rslt.embedPts(inscIdx,3),20,'g','filled');
title('SynCut');
embedPts = tsne(rslt.embedPts(:,1:3),[],rslt.embedPts(:,1:3),10);
% embedPts = tsne(rslt.embedPts(:,1:3),[],rslt.embedPts(:,1:3),10);
% embedPts = embedPts(:,[1,3,2]);

compareFig = figure('Position', [662,638,1261,726]);
subplot(1,2,2);
plot(embedPts(foliIdx,1), embedPts(foliIdx,2), '+', 'Color', colorsList(1,:));
hold on
plot(embedPts(frugIdx,1), embedPts(frugIdx,2), 'o', 'Color', colorsList(2,:));
plot(embedPts(inscIdx,1), embedPts(inscIdx,2), 'x', 'Color', colorsList(3,:));
axis equal
axis square
legend('Folivore', 'Frugivore', 'Insectivore', 'Location', 'NorthWest');
title('\textbf{(b) SynCut}', 'Interpreter', 'Latex', 'FontSize', 18);

speciesFig = figure('Position', [662,638,1261,726]);
subplot(1,2,2);
for j=1:length(GroupNames)
    plot(embedPts(TaxaByGroup{j},1), embedPts(TaxaByGroup{j},2), markerList{j}, 'Color', colorsList(j,:));
    if j==1
        hold on
    end
end
axis equal
axis square
legend(GroupNames, 'Location', 'NorthWest');
title('\textbf{(b) SynCut}', 'Interpreter', 'Latex', 'FontSize', 18);

%% compare with diffusion maps
embedPts = tsne(Ydm(:,1:3),[],rslt.embedPts,10);
% embedPts = tsne(Ydm(:,1:2),[],Ydm(:,1:3),10);

figure(compareFig);
subplot(1,2,1);
plot(embedPts(foliIdx,1), embedPts(foliIdx,2), '+', 'Color', colorsList(1,:));
hold on
plot(embedPts(frugIdx,1), embedPts(frugIdx,2), 'o', 'Color', colorsList(2,:));
plot(embedPts(inscIdx,1), embedPts(inscIdx,2), 'x', 'Color', colorsList(3,:));
axis equal
axis square
legend('Folivore', 'Frugivore', 'Insectivore', 'Location', 'NorthWest');
title('\textbf{(a) Diffusion Maps}', 'Interpreter', 'Latex', 'FontSize', 18);

figure(speciesFig);
subplot(1,2,1);
for j=1:length(GroupNames)
    plot(embedPts(TaxaByGroup{j},1), embedPts(TaxaByGroup{j},2), markerList{j}, 'Color', colorsList(j,:));
    if j==1
        hold on
    end
end
axis equal
axis square
legend(GroupNames, 'Location', 'NorthWest');
title('\textbf{(a) Diffusion Maps}', 'Interpreter', 'Latex', 'FontSize', 18);
