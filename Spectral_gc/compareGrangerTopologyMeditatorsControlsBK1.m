function stats = compareGrangerTopologyMeditatorsControlsBK1(protocolName,badEyeCondition,badTrialVersion,pairedDataFlag,saveOutputFlag,numPermutations,bandInfo,conditionType,hubZThreshold)
% compareGrangerTopologyMeditatorsControlsBK1 Compare Granger topology between BK1 groups.
%
% This function quantifies whether the structure of directed Granger
% connections changes differently in meditators and controls. The analysis
% is performed at the 64-electrode level using weighted directed graphs
% derived from the saved Granger matrices.
%
% Node-level metrics:
%   1. Out-degree   : weighted driver strength (sum of outgoing GC)
%   2. In-degree    : weighted receiver strength (sum of incoming GC)
%   3. Asymmetry    : (out - in) ./ (out + in)
%   4. Hubness      : z-scored total degree within each subject
%
% Global topology summaries:
%   1. Out-degree Gini coefficient
%   2. In-degree Gini coefficient
%   3. Mean absolute node asymmetry
%   4. Hub mass in top 10 percent of nodes
%   5. Number of hubs above hubZThreshold
%
% Statistics:
%   - Paired sign-flip permutation test for matched meditator-control pairs
%   - Unpaired label-shuffle permutation test for unmatched groups
%   - Benjamini-Hochberg FDR correction across electrodes for each metric
%
% Figures:
%   For each band, saves a composite figure with:
%   - meditator / control / group-difference channel heatmaps
%   - topoplots of group differences in node metrics
%   - jittered group distributions of global topology summaries
%
% Examples:
%   stats = compareGrangerTopologyMeditatorsControlsBK1('M1');
%   stats = compareGrangerTopologyMeditatorsControlsBK1('M1','ep','v8',1);
%   bandInfo(1).name = 'gamma'; bandInfo(1).range = [30 80];
%   stats = compareGrangerTopologyMeditatorsControlsBK1('M1','ep','v8',0,1,5000,bandInfo,'diff',1.0);

if ~exist('protocolName','var') || isempty(protocolName);       protocolName = 'M1';    end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition); badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion); badTrialVersion = 'v8'; end
if ~exist('pairedDataFlag','var') || isempty(pairedDataFlag);   pairedDataFlag = 0;     end
if ~exist('saveOutputFlag','var') || isempty(saveOutputFlag);   saveOutputFlag = 1;     end
if ~exist('numPermutations','var') || isempty(numPermutations); numPermutations = 5000; end
if ~exist('bandInfo','var') || isempty(bandInfo);               bandInfo = getDefaultBands(); end
if ~exist('conditionType','var') || isempty(conditionType);     conditionType = 'diff'; end
if ~exist('hubZThreshold','var') || isempty(hubZThreshold);     hubZThreshold = 1;      end

conditionType = validateConditionType(conditionType);
bandInfo = validateBandInfo(bandInfo);
bandTag = getBandInfoTag(bandInfo);

projectRoot = getProjectRoot();
grangerFolder = fullfile(fileparts(mfilename('fullpath')),'savedDataGranger');
summaryFolder = fullfile(grangerFolder,'summary');
statsFolder = fullfile(grangerFolder,'topologyStats');
figureFolder = fullfile(statsFolder,'figures');

ensureFolder(summaryFolder);
ensureFolder(statsFolder);
ensureFolder(figureFolder);

summaryFileName = fullfile(summaryFolder,[protocolName '_' badEyeCondition '_' badTrialVersion '_' bandTag '_grangerBandData.mat']);
if exist(summaryFileName,'file')
    tmp = load(summaryFileName);
    results = tmp.results;
else
    results = extractGrangerBandDataAllSubjectsBK1(protocolName,badEyeCondition,badTrialVersion,1,bandInfo);
end

if pairedDataFlag
    [medIndices,ctrlIndices,pairNames] = getValidPairedIndices(results);
    comparisonLabel = 'paired';
else
    [medIndices,ctrlIndices] = getValidUnpairedIndices(results);
    pairNames = {};
    comparisonLabel = 'unpaired';
end

if isempty(medIndices) || isempty(ctrlIndices)
    error('No valid subjects available for %s comparison.',comparisonLabel);
end

if saveOutputFlag
    addEEGLABTopoplot(projectRoot);
    [montageChanlocs,electrodeLabels] = getMontageChanlocs(projectRoot);
else
    [~,electrodeLabels] = getMontageChanlocs(projectRoot);
    montageChanlocs = [];
end

stats = struct;
stats.protocolName = protocolName;
stats.badEyeCondition = badEyeCondition;
stats.badTrialVersion = badTrialVersion;
stats.bandTag = bandTag;
stats.conditionType = conditionType;
stats.conditionLabel = getConditionLabel(conditionType);
stats.hubZThreshold = hubZThreshold;
stats.pairedDataFlag = pairedDataFlag;
stats.comparisonLabel = comparisonLabel;
stats.numPermutations = numPermutations;
stats.summaryFileName = summaryFileName;
stats.outputFileName = fullfile(statsFolder,[protocolName '_' badEyeCondition '_' badTrialVersion '_' bandTag '_' comparisonLabel '_' conditionType '_grangerTopologyStats.mat']);
stats.meditatorIndices = medIndices(:);
stats.controlIndices = ctrlIndices(:);
stats.meditatorSubjects = results.subjectNameList(medIndices);
stats.controlSubjects = results.subjectNameList(ctrlIndices);
stats.pairNames = pairNames;
stats.electrodeLabels = electrodeLabels;
stats.bandInfo = results.bandInfo;
stats.numMeditators = numel(medIndices);
stats.numControls = numel(ctrlIndices);

metricSpecs = getNodeMetricSpecs();
globalSpecs = getGlobalMetricSpecs();
numBands = numel(results.bandData);

for iBand = 1:numBands
    Gpre = results.bandData(iBand).channelPre;
    Gpost = results.bandData(iBand).channelPost;
    Gdiff = results.bandData(iBand).channelDiff;
    bandTopology = computeBandTopology(Gpre,Gpost,hubZThreshold);

    bandStats = struct;
    bandStats.name = results.bandData(iBand).name;
    bandStats.range = results.bandData(iBand).range;
    bandStats.conditionType = conditionType;
    bandStats.conditionLabel = stats.conditionLabel;
    bandStats.metricSpecs = metricSpecs;
    bandStats.globalSpecs = globalSpecs;
    bandStats.figureFileName = fullfile(figureFolder,[protocolName '_' results.bandData(iBand).name '_' bandTag '_' comparisonLabel '_' conditionType '_topology.png']);

    bandStats.meditatorMeanMatrix = mean(selectMatrixCondition(Gpre,Gpost,Gdiff,medIndices,conditionType),3,'omitnan');
    bandStats.controlMeanMatrix = mean(selectMatrixCondition(Gpre,Gpost,Gdiff,ctrlIndices,conditionType),3,'omitnan');
    bandStats.groupDifferenceMatrix = bandStats.meditatorMeanMatrix - bandStats.controlMeanMatrix;

    for iMetric = 1:numel(metricSpecs)
        metricName = metricSpecs(iMetric).name;
        metricValues = bandTopology.node.(metricName).(conditionType);
        medValues = metricValues(:,medIndices);
        ctrlValues = metricValues(:,ctrlIndices);

        metricStats = struct;
        metricStats.name = metricName;
        metricStats.label = metricSpecs(iMetric).label;
        metricStats.description = metricSpecs(iMetric).description;
        metricStats.meditatorValues = medValues;
        metricStats.controlValues = ctrlValues;
        metricStats.meditatorMean = mean(medValues,2,'omitnan');
        metricStats.controlMean = mean(ctrlValues,2,'omitnan');
        metricStats.groupDifference = metricStats.meditatorMean - metricStats.controlMean;
        metricStats.sampleCount = nan(numel(electrodeLabels),1);
        metricStats.pValue = nan(numel(electrodeLabels),1);
        metricStats.effectSize = nan(numel(electrodeLabels),1);

        for iElectrode = 1:numel(electrodeLabels)
            x = medValues(iElectrode,:)';
            y = ctrlValues(iElectrode,:)';
            if pairedDataFlag
                [pVal,effectSize,nUsed] = pairedPermutationTest(x,y,numPermutations);
            else
                [pVal,effectSize,nUsed] = unpairedPermutationTest(x,y,numPermutations);
            end
            metricStats.pValue(iElectrode) = pVal;
            metricStats.effectSize(iElectrode) = effectSize;
            metricStats.sampleCount(iElectrode) = nUsed;
        end

        [metricStats.qValue,metricStats.significantMask] = bhFDR(metricStats.pValue,0.05);
        metricStats.significantElectrodes = getSignificantElectrodes(metricStats,electrodeLabels);
        bandStats.nodeStats.(metricName) = metricStats;
    end

    for iMetric = 1:numel(globalSpecs)
        metricName = globalSpecs(iMetric).name;
        globalValues = bandTopology.global.(metricName).(conditionType);
        x = globalValues(medIndices);
        y = globalValues(ctrlIndices);

        if pairedDataFlag
            [pVal,effectSize,nUsed] = pairedPermutationTest(x(:),y(:),numPermutations);
        else
            [pVal,effectSize,nUsed] = unpairedPermutationTest(x(:),y(:),numPermutations);
        end

        bandStats.globalStats(iMetric).name = metricName;
        bandStats.globalStats(iMetric).label = globalSpecs(iMetric).label;
        bandStats.globalStats(iMetric).description = globalSpecs(iMetric).description;
        bandStats.globalStats(iMetric).meditatorValues = x(:);
        bandStats.globalStats(iMetric).controlValues = y(:);
        bandStats.globalStats(iMetric).meditatorMean = mean(x,'omitnan');
        bandStats.globalStats(iMetric).controlMean = mean(y,'omitnan');
        bandStats.globalStats(iMetric).groupDifference = bandStats.globalStats(iMetric).meditatorMean - bandStats.globalStats(iMetric).controlMean;
        bandStats.globalStats(iMetric).pValue = pVal;
        bandStats.globalStats(iMetric).effectSize = effectSize;
        bandStats.globalStats(iMetric).sampleCount = nUsed;
    end

    [globalQ,globalSig] = bhFDR([bandStats.globalStats.pValue],0.05);
    for iMetric = 1:numel(globalSpecs)
        bandStats.globalStats(iMetric).qValue = globalQ(iMetric);
        bandStats.globalStats(iMetric).significant = globalSig(iMetric);
    end

    bandStats.topology = bandTopology;
    stats.bandStats(iBand) = bandStats;

    if saveOutputFlag
        hFig = createBandFigure(bandStats,montageChanlocs,electrodeLabels,comparisonLabel,stats.numMeditators,stats.numControls,pairedDataFlag);
        saveas(hFig,bandStats.figureFileName);
        close(hFig);
    end
end

if saveOutputFlag
    save(stats.outputFileName,'stats','-v7.3');
end

disp(sprintf('Saved %s topology stats for %s using %d meditators and %d controls (%s).',comparisonLabel,protocolName,stats.numMeditators,stats.numControls,stats.conditionLabel)); %#ok<DSPS>
end

function bandTopology = computeBandTopology(GpreAll,GpostAll,hubZThreshold)
numNodes = size(GpreAll,1);
numSubjects = size(GpreAll,3);
metricSpecs = getNodeMetricSpecs();
globalSpecs = getGlobalMetricSpecs();

for iMetric = 1:numel(metricSpecs)
    metricName = metricSpecs(iMetric).name;
    bandTopology.node.(metricName).pre = nan(numNodes,numSubjects);
    bandTopology.node.(metricName).post = nan(numNodes,numSubjects);
    bandTopology.node.(metricName).diff = nan(numNodes,numSubjects);
end

for iMetric = 1:numel(globalSpecs)
    metricName = globalSpecs(iMetric).name;
    bandTopology.global.(metricName).pre = nan(numSubjects,1);
    bandTopology.global.(metricName).post = nan(numSubjects,1);
    bandTopology.global.(metricName).diff = nan(numSubjects,1);
end

for iSubject = 1:numSubjects
    preMetrics = computeSingleSubjectTopology(GpreAll(:,:,iSubject),hubZThreshold);
    postMetrics = computeSingleSubjectTopology(GpostAll(:,:,iSubject),hubZThreshold);

    for iMetric = 1:numel(metricSpecs)
        metricName = metricSpecs(iMetric).name;
        bandTopology.node.(metricName).pre(:,iSubject) = preMetrics.node.(metricName);
        bandTopology.node.(metricName).post(:,iSubject) = postMetrics.node.(metricName);
        bandTopology.node.(metricName).diff(:,iSubject) = postMetrics.node.(metricName) - preMetrics.node.(metricName);
    end

    for iMetric = 1:numel(globalSpecs)
        metricName = globalSpecs(iMetric).name;
        bandTopology.global.(metricName).pre(iSubject) = preMetrics.global.(metricName);
        bandTopology.global.(metricName).post(iSubject) = postMetrics.global.(metricName);
        bandTopology.global.(metricName).diff(iSubject) = postMetrics.global.(metricName) - preMetrics.global.(metricName);
    end
end
end

function subjectTopology = computeSingleSubjectTopology(G,hubZThreshold)
G = removeSelfConnections(G);
outDegree = sumKeepNaN(G,2);
inDegree = sumKeepNaN(G,1)';
totalDegree = outDegree + inDegree;

asymmetry = nan(size(outDegree));
denom = outDegree + inDegree;
validMask = isfinite(denom) & abs(denom) > eps;
asymmetry(validMask) = (outDegree(validMask) - inDegree(validMask)) ./ denom(validMask);

subjectTopology.node.outDegree = outDegree;
subjectTopology.node.inDegree = inDegree;
subjectTopology.node.asymmetry = asymmetry;
subjectTopology.node.hubness = safeZScore(totalDegree);

subjectTopology.global.outDegreeGini = giniCoefficient(outDegree);
subjectTopology.global.inDegreeGini = giniCoefficient(inDegree);
subjectTopology.global.meanAbsAsymmetry = mean(abs(asymmetry), 'omitnan');
subjectTopology.global.hubMassTop10 = topKMass(totalDegree,0.10);
subjectTopology.global.numHubs = sum(subjectTopology.node.hubness >= hubZThreshold & isfinite(subjectTopology.node.hubness));
end

function matrixOut = selectMatrixCondition(Gpre,Gpost,Gdiff,indices,conditionType)
switch conditionType
    case 'pre'
        matrixOut = Gpre(:,:,indices);
    case 'post'
        matrixOut = Gpost(:,:,indices);
    case 'diff'
        matrixOut = Gdiff(:,:,indices);
    otherwise
        error('Unknown conditionType: %s',conditionType);
end
end

function hFig = createBandFigure(bandStats,montageChanlocs,electrodeLabels,comparisonLabel,numMeditators,numControls,pairedDataFlag)
hFig = figure('Color','w','Name',[bandStats.name ' topology']);
hFig.Position(3:4) = [1700 1050];

groupCLim = getGroupMatrixCLim(bandStats);
diffCLim = symmetricCLim(bandStats.groupDifferenceMatrix);

subplot(3,4,1);
plotHeatmap(gca,bandStats.meditatorMeanMatrix,electrodeLabels,groupCLim);
title(sprintf('Meditator Mean (%d)',numMeditators),'Interpreter','none');

subplot(3,4,2);
plotHeatmap(gca,bandStats.controlMeanMatrix,electrodeLabels,groupCLim);
title(sprintf('Control Mean (%d)',numControls),'Interpreter','none');

subplot(3,4,3);
plotHeatmap(gca,bandStats.groupDifferenceMatrix,electrodeLabels,diffCLim);
title('Meditator - Control','Interpreter','none');

subplot(3,4,4);
plotSummaryText(gca,bandStats,comparisonLabel,numMeditators,numControls);

topoCLims = getMetricCLims(bandStats);
nodeMetricNames = {'outDegree','inDegree','asymmetry','hubness'};
for iMetric = 1:numel(nodeMetricNames)
    subplot(3,4,4+iMetric);
    metricStats = bandStats.nodeStats.(nodeMetricNames{iMetric});
    plotTopologyTopoplot(gca,metricStats,montageChanlocs,topoCLims.(nodeMetricNames{iMetric}));
end

plotGlobalMetricDistribution(gcaForSubplot(3,4,9),bandStats.globalStats(1),pairedDataFlag);
plotGlobalMetricDistribution(gcaForSubplot(3,4,10),bandStats.globalStats(2),pairedDataFlag);
plotGlobalMetricDistribution(gcaForSubplot(3,4,11),bandStats.globalStats(3),pairedDataFlag);
plotGlobalMetricDistribution(gcaForSubplot(3,4,12),bandStats.globalStats(4),pairedDataFlag);

if exist('sgtitle','file') == 2
    sgtitle(sprintf('%s band (%g-%g Hz) | %s | %s topology',bandStats.name,bandStats.range(1),bandStats.range(2),comparisonLabel,bandStats.conditionLabel),'Interpreter','none');
end
end

function ax = gcaForSubplot(m,n,p)
ax = subplot(m,n,p);
end

function plotSummaryText(ax,bandStats,comparisonLabel,numMeditators,numControls)
axis(ax,'off');
numSig = zeros(1,numel(bandStats.metricSpecs));
for iMetric = 1:numel(bandStats.metricSpecs)
    metricName = bandStats.metricSpecs(iMetric).name;
    numSig(iMetric) = sum(bandStats.nodeStats.(metricName).significantMask);
end

summaryLines = {
    sprintf('Comparison: %s',comparisonLabel)
    sprintf('Condition: %s',bandStats.conditionLabel)
    sprintf('Meditators: %d',numMeditators)
    sprintf('Controls: %d',numControls)
    ' '
    'Significant electrodes (q < 0.05):'
    sprintf('Out-degree: %d',numSig(1))
    sprintf('In-degree: %d',numSig(2))
    sprintf('Asymmetry: %d',numSig(3))
    sprintf('Hubness: %d',numSig(4))
    ' '
    'Global metrics plotted below:'
    'Out-Gini, In-Gini, |Asym|, Hub-mass'
    };

summaryText = sprintf('%s\n',summaryLines{:});
text(0,1,summaryText,'Parent',ax,'VerticalAlignment','top','FontName','Helvetica','FontSize',11);
end

function plotTopologyTopoplot(ax,metricStats,montageChanlocs,cLims)
axes(ax); %#ok<LAXES>
topoplot(metricStats.groupDifference,montageChanlocs,'electrodes','on','plotrad',0.6,'headrad',0.6);
colormap(ax,parula);
colorbar(ax);
if ~isempty(cLims) && all(isfinite(cLims))
    caxis(ax,cLims);
end
sigElectrodes = find(metricStats.significantMask);
if ~isempty(sigElectrodes)
    hold(ax,'on');
    try
        topoplot(metricStats.groupDifference,montageChanlocs,'style','blank','electrodes','off','plotrad',0.6,'headrad',0.6,...
            'emarker2',{sigElectrodes,'o','k',8});
    catch
        scatterTopoplotMarkers(ax,montageChanlocs,sigElectrodes);
    end
    hold(ax,'off');
end
title(sprintf('%s\nq<0.05 electrodes: %d',metricStats.label,sum(metricStats.significantMask)),'Interpreter','none');
end

function scatterTopoplotMarkers(ax,montageChanlocs,markerIndices)
numMarkers = numel(markerIndices);
xVals = nan(numMarkers,1);
yVals = nan(numMarkers,1);

for iMarker = 1:numMarkers
    idx = markerIndices(iMarker);
    if isfield(montageChanlocs,'X') && isfield(montageChanlocs,'Y')
        xVals(iMarker) = montageChanlocs(idx).X;
        yVals(iMarker) = montageChanlocs(idx).Y;
    else
        if isfield(montageChanlocs,'theta')
            theta = montageChanlocs(idx).theta;
        else
            theta = 0;
        end
        if isfield(montageChanlocs,'radius')
            radius = montageChanlocs(idx).radius;
        else
            radius = 0;
        end
        xVals(iMarker) = radius * sind(theta);
        yVals(iMarker) = radius * cosd(theta);
    end
end

validMask = isfinite(xVals) & isfinite(yVals);
scatter(ax,xVals(validMask),yVals(validMask),64,'ko','LineWidth',1.5);
end

function plotGlobalMetricDistribution(ax,metricStats,pairedDataFlag)
axes(ax); %#ok<LAXES>
cla(ax);
hold(ax,'on');

x = metricStats.meditatorValues(:);
y = metricStats.controlValues(:);
xPos1 = 1 + 0.10*(rand(size(x))-0.5);
xPos2 = 2 + 0.10*(rand(size(y))-0.5);

if pairedDataFlag
    nPairs = min(numel(x),numel(y));
    for iPair = 1:nPairs
        if isfinite(x(iPair)) && isfinite(y(iPair))
            plot(ax,[1 2],[x(iPair) y(iPair)],'-','Color',[0.8 0.8 0.8],'LineWidth',0.75);
        end
    end
end

scatter(ax,xPos1,x,30,[0.2 0.45 0.85],'filled');
scatter(ax,xPos2,y,30,[0.85 0.4 0.25],'filled');

plotMeanAndSEM(ax,1,x,[0.1 0.2 0.55]);
plotMeanAndSEM(ax,2,y,[0.55 0.2 0.1]);

set(ax,'XLim',[0.5 2.5],'XTick',[1 2],'XTickLabel',{'Meditator','Control'},'TickDir','out');
grid(ax,'on');
title(ax,sprintf('%s\np=%.4f | q=%.4f',metricStats.label,metricStats.pValue,metricStats.qValue),'Interpreter','none');
ylabel(ax,metricStats.label,'Interpreter','none');
hold(ax,'off');
end

function plotMeanAndSEM(ax,xPos,data,colorVal)
data = data(isfinite(data));
if isempty(data)
    return;
end
mu = mean(data,'omitnan');
semVal = std(data,'omitnan') ./ sqrt(numel(data));
plot(ax,[xPos xPos],[mu-semVal mu+semVal],'-','Color',colorVal,'LineWidth',2);
scatter(ax,xPos,mu,60,colorVal,'filled');
end

function plotHeatmap(ax,matrix,tickLabels,cLims)
imagesc(ax,matrix,'AlphaData',isfinite(matrix));
set(ax,'YDir','normal');
axis(ax,'square');
colorbar(ax);
colormap(ax,parula);
set(ax,'XTick',1:numel(tickLabels),'XTickLabel',tickLabels,'YTick',1:numel(tickLabels),'YTickLabel',tickLabels,'TickDir','out');
xtickangle(ax,90);
if ~isempty(cLims) && all(isfinite(cLims))
    caxis(ax,cLims);
end
end

function metricCLims = getMetricCLims(bandStats)
nodeMetricNames = {'outDegree','inDegree','asymmetry','hubness'};
for iMetric = 1:numel(nodeMetricNames)
    metricName = nodeMetricNames{iMetric};
    metricCLims.(metricName) = symmetricCLim(bandStats.nodeStats.(metricName).groupDifference);
end
end

function cLim = getGroupMatrixCLim(bandStats)
if strcmp(bandStats.conditionType,'diff')
    combinedVals = cat(1,bandStats.meditatorMeanMatrix(:),bandStats.controlMeanMatrix(:));
    combinedVals = combinedVals(isfinite(combinedVals));
    if isempty(combinedVals)
        cLim = [];
    else
        maxAbs = max(abs(combinedVals));
        cLim = [-maxAbs maxAbs];
    end
else
    cLim = [0 max([maxFiniteValue(bandStats.meditatorMeanMatrix) maxFiniteValue(bandStats.controlMeanMatrix)])];
end
end

function significantElectrodes = getSignificantElectrodes(metricStats,electrodeLabels)
sigIdx = find(metricStats.significantMask);
significantElectrodes = cell(numel(sigIdx),5);
for iSig = 1:numel(sigIdx)
    idx = sigIdx(iSig);
    significantElectrodes{iSig,1} = electrodeLabels{idx};
    significantElectrodes{iSig,2} = metricStats.groupDifference(idx);
    significantElectrodes{iSig,3} = metricStats.pValue(idx);
    significantElectrodes{iSig,4} = metricStats.qValue(idx);
    significantElectrodes{iSig,5} = metricStats.effectSize(idx);
end
end

function [pVal,effectSize,nUsed] = pairedPermutationTest(x,y,numPermutations)
validMask = isfinite(x) & isfinite(y);
x = x(validMask);
y = y(validMask);
nUsed = numel(x);

if nUsed < 2
    pVal = NaN;
    effectSize = NaN;
    return;
end

d = x - y;
observed = mean(d);
permVals = nan(numPermutations,1);
for iPerm = 1:numPermutations
    signs = 2*(rand(nUsed,1) > 0.5) - 1;
    permVals(iPerm) = mean(d .* signs);
end

pVal = (1 + sum(abs(permVals) >= abs(observed))) / (numPermutations + 1);
effectSize = standardizedMean(d);
end

function [pVal,effectSize,nUsed] = unpairedPermutationTest(x,y,numPermutations)
x = x(isfinite(x));
y = y(isfinite(y));
nX = numel(x);
nY = numel(y);
nUsed = nX + nY;

if nX < 2 || nY < 2
    pVal = NaN;
    effectSize = NaN;
    return;
end

observed = mean(x) - mean(y);
pooled = [x; y];
permVals = nan(numPermutations,1);
for iPerm = 1:numPermutations
    order = randperm(numel(pooled));
    xPerm = pooled(order(1:nX));
    yPerm = pooled(order(nX+1:end));
    permVals(iPerm) = mean(xPerm) - mean(yPerm);
end

pVal = (1 + sum(abs(permVals) >= abs(observed))) / (numPermutations + 1);
effectSize = cohenD(x,y);
end

function effectSize = standardizedMean(x)
sdVal = std(x);
if isempty(sdVal) || ~isfinite(sdVal) || sdVal == 0
    effectSize = NaN;
else
    effectSize = mean(x) ./ sdVal;
end
end

function effectSize = cohenD(x,y)
nX = numel(x);
nY = numel(y);
if nX < 2 || nY < 2
    effectSize = NaN;
    return;
end

varX = var(x,0);
varY = var(y,0);
pooledStd = sqrt(((nX - 1)*varX + (nY - 1)*varY) / (nX + nY - 2));
if ~isfinite(pooledStd) || pooledStd == 0
    effectSize = NaN;
else
    effectSize = (mean(x) - mean(y)) / pooledStd;
end
end

function [qMatrix,sigMask] = bhFDR(pMatrix,alpha)
pVec = pMatrix(:);
qVec = nan(size(pVec));
sigMask = false(size(pVec));

validMask = isfinite(pVec);
validP = pVec(validMask);

if isempty(validP)
    qMatrix = reshape(qVec,size(pMatrix));
    sigMask = reshape(sigMask,size(pMatrix));
    return;
end

[sortedP,sortOrder] = sort(validP(:));
m = numel(sortedP);
qSorted = sortedP .* m ./ (1:m)';
qSorted = flipud(cummin(flipud(qSorted)));
qSorted(qSorted > 1) = 1;

qValid = nan(size(validP));
qValid(sortOrder) = qSorted;
qVec(validMask) = qValid;
sigMask(validMask) = qValid <= alpha;

qMatrix = reshape(qVec,size(pMatrix));
sigMask = reshape(sigMask,size(pMatrix));
end

function z = safeZScore(x)
z = nan(size(x));
validMask = isfinite(x);
if nnz(validMask) < 2
    return;
end
mu = mean(x(validMask),'omitnan');
sdVal = std(x(validMask),0,'omitnan');
if ~isfinite(sdVal) || sdVal == 0
    return;
end
z(validMask) = (x(validMask) - mu) ./ sdVal;
end

function g = giniCoefficient(x)
x = x(isfinite(x));
if isempty(x)
    g = NaN;
    return;
end

x = sort(x(:));
if any(x < 0)
    x = x - min(x);
end
if sum(x) == 0
    g = 0;
    return;
end

n = numel(x);
g = (2 * sum((1:n)' .* x)) / (n * sum(x)) - (n + 1) / n;
end

function hubMass = topKMass(x,fractionToKeep)
x = x(isfinite(x));
if isempty(x)
    hubMass = NaN;
    return;
end

if any(x < 0)
    x = x - min(x);
end
denom = sum(x);
if denom == 0
    hubMass = 0;
    return;
end

nKeep = max(1,ceil(fractionToKeep * numel(x)));
x = sort(x,'descend');
hubMass = sum(x(1:nKeep)) ./ denom;
end

function G = removeSelfConnections(G)
n = size(G,1);
G(1:n+1:end) = NaN;
end

function s = sumKeepNaN(x,dim)
s = sum(x,dim,'omitnan');
allNaNMask = all(~isfinite(x),dim);
s(allNaNMask) = NaN;
end

function metricSpecs = getNodeMetricSpecs()
metricSpecs(1).name = 'outDegree';
metricSpecs(1).label = 'Out-degree';
metricSpecs(1).description = 'Weighted outgoing GC strength';

metricSpecs(2).name = 'inDegree';
metricSpecs(2).label = 'In-degree';
metricSpecs(2).description = 'Weighted incoming GC strength';

metricSpecs(3).name = 'asymmetry';
metricSpecs(3).label = 'Asymmetry';
metricSpecs(3).description = '(Out - In) / (Out + In)';

metricSpecs(4).name = 'hubness';
metricSpecs(4).label = 'Hubness';
metricSpecs(4).description = 'Z-scored total degree';
end

function metricSpecs = getGlobalMetricSpecs()
metricSpecs(1).name = 'outDegreeGini';
metricSpecs(1).label = 'Out-degree Gini';
metricSpecs(1).description = 'Driver concentration across electrodes';

metricSpecs(2).name = 'inDegreeGini';
metricSpecs(2).label = 'In-degree Gini';
metricSpecs(2).description = 'Receiver concentration across electrodes';

metricSpecs(3).name = 'meanAbsAsymmetry';
metricSpecs(3).label = 'Mean |Asymmetry|';
metricSpecs(3).description = 'Average driver-receiver imbalance';

metricSpecs(4).name = 'hubMassTop10';
metricSpecs(4).label = 'Hub mass top 10%';
metricSpecs(4).description = 'Fraction of total degree in top 10 percent nodes';

metricSpecs(5).name = 'numHubs';
metricSpecs(5).label = 'Number of hubs';
metricSpecs(5).description = 'Electrodes with hubness above threshold';
end

function [medIndices,ctrlIndices] = getValidUnpairedIndices(results)
validMask = results.validSubjectMask(:);
groupList = results.groupNameList(:);

medIndices = find(validMask & strcmp(groupList,'Meditator'));
ctrlIndices = find(validMask & strcmp(groupList,'Control'));
end

function [medIndices,ctrlIndices,pairNames] = getValidPairedIndices(results)
projectRoot = getProjectRoot();
pairFunctionFolder = fullfile(projectRoot,'ProjectDhyaanBK1Programs','commonAnalysisCodes','informationFiles');
addpath(pairFunctionFolder);
cleanupObj = onCleanup(@() rmpath(pairFunctionFolder)); %#ok<NASGU>
pairedSubjectNameList = getPairedSubjectsBK1();

medIndices = [];
ctrlIndices = [];
pairNames = {};

for iPair = 1:size(pairedSubjectNameList,1)
    medName = pairedSubjectNameList{iPair,1};
    ctrlName = pairedSubjectNameList{iPair,2};

    medIndex = find(strcmp(results.subjectNameList,medName),1);
    ctrlIndex = find(strcmp(results.subjectNameList,ctrlName),1);

    if isempty(medIndex) || isempty(ctrlIndex)
        continue;
    end

    if ~results.validSubjectMask(medIndex) || ~results.validSubjectMask(ctrlIndex)
        continue;
    end

    medIndices(end+1,1) = medIndex; %#ok<AGROW>
    ctrlIndices(end+1,1) = ctrlIndex; %#ok<AGROW>
    pairNames{end+1,1} = sprintf('%s_vs_%s',medName,ctrlName); %#ok<AGROW>
end
end

function bandInfo = getDefaultBands()
bandInfo(1).name = 'alpha';     bandInfo(1).range = [7 10];
bandInfo(2).name = 'beta';      bandInfo(2).range = [20 32];
bandInfo(3).name = 'highgamma'; bandInfo(3).range = [30 80];
end

function bandInfo = validateBandInfo(bandInfo)
if ~isstruct(bandInfo) || isempty(bandInfo)
    error('bandInfo must be a non-empty struct array with fields name and range.');
end

for iBand = 1:numel(bandInfo)
    if ~isfield(bandInfo(iBand),'name') || ~isfield(bandInfo(iBand),'range')
        error('Each bandInfo entry must contain name and range.');
    end
    if numel(bandInfo(iBand).range) ~= 2
        error('Each bandInfo range must have exactly two values.');
    end
    bandInfo(iBand).range = reshape(bandInfo(iBand).range,1,2);
    if bandInfo(iBand).range(1) >= bandInfo(iBand).range(2)
        error('Band range must satisfy low < high.');
    end
    bandInfo(iBand).name = char(string(bandInfo(iBand).name));
end
end

function conditionType = validateConditionType(conditionType)
conditionType = lower(char(conditionType));
switch conditionType
    case {'pre','post'}
        return;
    case {'diff','post-pre'}
        conditionType = 'diff';
    otherwise
        error('conditionType must be ''pre'', ''post'', or ''diff''.');
end
end

function conditionLabel = getConditionLabel(conditionType)
switch conditionType
    case 'pre'
        conditionLabel = 'Pre';
    case 'post'
        conditionLabel = 'Post';
    case 'diff'
        conditionLabel = 'Post - Pre';
    otherwise
        conditionLabel = conditionType;
end
end

function bandTag = getBandInfoTag(bandInfo)
tagParts = cell(1,numel(bandInfo));
for iBand = 1:numel(bandInfo)
    safeName = regexprep(lower(bandInfo(iBand).name),'[^a-z0-9]+','');
    rangeVals = round(100*bandInfo(iBand).range);
    tagParts{iBand} = sprintf('%s_%g_%g',safeName,rangeVals(1),rangeVals(2));
end
bandTag = strjoin(tagParts,'__');
end

function cLim = symmetricCLim(matrix)
finiteVals = matrix(isfinite(matrix));
if isempty(finiteVals)
    cLim = [];
else
    maxAbs = max(abs(finiteVals));
    cLim = [-maxAbs maxAbs];
end
end

function v = maxFiniteValue(matrix)
finiteVals = matrix(isfinite(matrix));
if isempty(finiteVals)
    v = NaN;
else
    v = max(finiteVals);
end
end

function projectRoot = getProjectRoot()
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

function [montageChanlocs,montageLabels] = getMontageChanlocs(projectRoot)
layoutFolder = fullfile(projectRoot,'Montages','Layouts','actiCap64_UOL');
x = load(fullfile(layoutFolder,'actiCap64_UOL.mat'));
y = load(fullfile(layoutFolder,'actiCap64_UOLLabels.mat'));
montageChanlocs = x.chanlocs;
montageLabels = y.montageLabels(:,2);
end

function addEEGLABTopoplot(projectRoot)
eeglabFolder = fullfile(projectRoot,'fieldtrip-20260211','external','eeglab');
if exist('topoplot','file') ~= 2
    addpath(eeglabFolder);
end
end

function ensureFolder(folderName)
if ~exist(folderName,'dir')
    mkdir(folderName);
end
end
