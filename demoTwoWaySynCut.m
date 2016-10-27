%% setup parameters
clearvars
close all
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

syncableFlag = false;
debugFlag    = true;
newCaseFlag  = true;

if ~newCaseFlag
    casePath = 'data/successCase05.mat';
end

N = 100;  %%% number of vertices in each cluster of the random graph
d = 5;   %%% dimension of the orthogonal group
numClusters = 2;
numLinks = randi([100,250]);  %%% number of cross-cluster links
degUB = 8;      %%% vertex degree upper bound withtin each cluster
degLB = 4;      %%% vertex degree lower bound withtin each cluster
ccType = 'unif';  %%%% ['nn' | 'unif']

maxIter = 10;
tol = 1e-8;
numKmeans = 200;
bandwidth = 1;
adjType = 'dis';
colorList = {'r','b','k','m'};
hsv = rgb2hsv(winter);
close(gcf);

params = struct('debugFlag', debugFlag,...
    'd', d,...
    'numClusters', numClusters,...
    'tol', tol,...
    'maxIter', maxIter,...
    'numKmeans', numKmeans,...
    'bandwidth', bandwidth,...
    'adjType', adjType,...
    'hsv', hsv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% generate random graph and edge potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if newCaseFlag
    G = genWeakExpander(N, numClusters, numLinks, [degLB,degUB], ccType);
        
    % generate "ground truth" vertex and edge potentials, and contaminate the
    % edge potential on cross-cluster links with random orthonormal matrices
    vertPotCell = cell(1,N*numClusters);
    for j=1:N*numClusters
        vertPotCell{j} = orth(rand(d));
    end
    
    edgePotCell = cell(size(G.adjMat));
    [rIdx, cIdx] = find(G.adjMat);
    for j=1:length(rIdx)
        edgePotCell{rIdx(j),cIdx(j)} = vertPotCell{rIdx(j)}*vertPotCell{cIdx(j)}';
    end
    if ~syncableFlag
        for j=1:length(G.ccRowIdx)
            edgePotCell{G.ccRowIdx(j),G.ccColIdx(j)} = orth(rand(d));
            edgePotCell{G.ccColIdx(j),G.ccRowIdx(j)} = edgePotCell{G.ccRowIdx(j),G.ccColIdx(j)}';
        end
    end
else
    %%%%% load previous case
    if exist('casePath', 'var')
        load(casePath);
        params.debugFlag = debugFlag;
    else
        error('casePath unspecified');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% run SynCut and collect result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if debugFlag
    [params.colorList] = deal(colorList);
    [params.vertPotCell] = deal(vertPotCell);
end

rslt = SynCut_TwoWay(G, edgePotCell, params);

rslt.errCounts = min(sum(rslt.clusterLabel{1} > N)+sum(rslt.clusterLabel{2} <= N),...
    sum(rslt.clusterLabel{1} <= N)+sum(rslt.clusterLabel{2} > N));
rslt.errRate = rslt.errCounts / size(G.adjMat, 1);


