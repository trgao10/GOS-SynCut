%% setup parameters
clearvars
close all
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

numRuns = 1000;

syncableFlag = false;
debugFlag    = false;
newCaseFlag  = true;

N = 100;  %%% number of vertices in each cluster of the random graph
d = 5;   %%% dimension of the orthogonal group
numClusters = 2;
% numLinks = 10;
numLinks = randi([100,250], [numRuns,1]);  %%% number of cross-cluster links
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

%% repeated trials
rslt = cell(1,numRuns);
cback = 0;
for k=1:numRuns
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% generate random graph
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G = genWeakExpander(N, numClusters, numLinks(k), [degLB,degUB], ccType);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% generate "ground truth" vertex and edge potential, and contaminate
    %%%% the edge potential on cross-cluster links with random orthonormal
    %%%% matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% run SynCut and collect result
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if debugFlag
        [params.colorList] = deal(colorList);
        [params.vertPotCell] = deal(vertPotCell);
    end

    try
        rslt{k} = SynCut(G, edgePotCell, params);
%         rslt{k} = SynCut_TwoWay(G, edgePotCell, params);
    catch E
        disp(E.message);
        save('data/bugCase.mat', 'G', 'vertPotCell', 'edgePotCell', 'params');
        continue
    end
    
    rslt{k}.errCounts = min(sum(rslt{k}.clusterLabel{1} > N)+sum(rslt{k}.clusterLabel{2} <= N),...
        sum(rslt{k}.clusterLabel{1} <= N)+sum(rslt{k}.clusterLabel{2} > N));
    rslt{k}.errRate = rslt{k}.errCounts / size(G.adjMat, 1);

    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('%4d/%d done.\n',k,numRuns);
end

rslt(cellfun('isempty',rslt)) = [];

%%%%% compare with normalized graph cut
cback = 0;
for j=1:length(rslt)
    NCutClusterLabel = specClusterWrapper(rslt{j}.G.adjMat,numClusters,...
        numKmeans,bandwidth,'sim');
    rslt{j}.NCutErrCounts = min(sum(NCutClusterLabel{1} > N)+sum(NCutClusterLabel{2} <= N),...
        sum(NCutClusterLabel{1} <= N)+sum(NCutClusterLabel{2} > N));
    rslt{j}.NCutErrRate = rslt{j}.NCutErrCounts / size(rslt{j}.G.adjMat,1);
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('NCut %4d/%d done.\n',j,length(rslt));
end

specGaps = cellfun(@(result) result.G.specGap, rslt);
errRates = cellfun(@(result) result.errRate, rslt);
NCutErrRates = cellfun(@(result) result.NCutErrRate, rslt);
iterCounts = cellfun(@(result) result.iterCounter, rslt);
numCCLinks = cellfun(@(result) length(result.G.ccRowIdx), rslt);


%%
NBins = 20;

figure;
subplot(1,4,1);
hist(specGaps,NBins);
title('(a) Spectral Gaps');
subplot(1,4,2);
hist(errRates,NBins);
xlim([0,0.5]);
ylim([0,numRuns]);
title('(b) SynCut Error Rates');
subplot(1,4,3);
hist(NCutErrRates,NBins);
ylim([0,numRuns]);
xlim([0,0.5]);
title('(c) NCut Error Rates');
subplot(1,4,4);
hist(iterCounts);
xlim([1,10]);
title('(d) Number of Iterations');

%%
figure;
subplot(1,3,1);
hist(numCCLinks);
title('(a) Distrib. of Num. of Inter-Comp. Links');
subplot(1,3,2);
hist(specGaps);
title('(b) Distrib. of Spectral Gaps')
subplot(1,3,3);
plot(numCCLinks, specGaps, 'bx');
xlabel('num of inter-component links');
ylabel('spectral gap');
title('(c) Num. of Inter-Comp. Links vs. Spectral Gap');

%%
validIdx = 1:length(NCutErrRates);

figure;
plot(specGaps(validIdx), NCutErrRates(validIdx), 'ro');
hold on
plot(specGaps(validIdx), errRates(validIdx), 'b+');
legend('NCut','SynCut');
xlabel('Spectral Gap');
ylabel('Error Ratio');

%%
% figure;
% hist3([specGaps(validIdx);iterCounts(validIdx)]');
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
% xlabel('Spectral Gap');
% ylabel('Number of SynCut Iterations');
% 
% for k=1:maxIter
%     currIdx = find(iterCounts <= k);
%     fprintf('%d/%d instances terminated within %d iterations\n',length(currIdx),numRuns,k);
%     fprintf('with error ratio statistics: mean = %.4f, std = %.4f\n',mean(errRates(currIdx)),std(errRates(currIdx)));
%     fprintf('same statistics for NCut: mean = %.4f, std = %.4f\n',mean(NCutErrRates(currIdx)),std(NCutErrRates(currIdx)));
% end

