function stats = compareGrangerCommunitiesMeditatorsControlsBK1(protocolName,badEyeCondition,badTrialVersion,pairedDataFlag,saveOutputFlag,numPermutations,bandInfo,conditionType,densityThreshold,gammaValue,numCommunityRuns)
% compareGrangerCommunitiesMeditatorsControlsBK1 Compare directed GC communities between BK1 groups.
%
% This function loads subject-wise Granger matrices that were already saved
% in connectivityProjectCodes/savedDataGranger, thresholds each subject's
% directed weighted graph at a fixed edge density, detects communities with
% repeated directed modularity optimization, and compares community
% organization between meditators and controls.
%
% Community metrics per subject:
%   1. Number of communities
%   2. Modularity Q
%   3. Mean out-participation coefficient
%   4. Mean in-participation coefficient
%   5. Largest-community fraction
%
% Stability metric:
%   Subject-wise mean normalized mutual information (NMI) to other subjects
%   in the same group, computed from the detected partitions.
%
% Node-level metrics:
%   1. Out-participation coefficient
%   2. In-participation coefficient
%
% Statistical testing:
%   - Paired sign-flip permutation test for matched meditator-control pairs
%   - Unpaired label-shuffle permutation test for unmatched groups
%   - Benjamini-Hochberg FDR correction across nodes for node metrics
%
% Notes:
%   - Community labels themselves are not directly compared by label value.
%   - For conditionType = 'diff', the scalar metrics are computed as
%     post - pre, and stability is computed as (post mean NMI) - (pre mean NMI).
%
% Examples:
%   stats = compareGrangerCommunitiesMeditatorsControlsBK1('M1');
%   stats = compareGrangerCommunitiesMeditatorsControlsBK1('M1','ep','v8',0,1,5000,[], 'post',0.20,1,50);

if ~exist('protocolName','var') || isempty(protocolName);          protocolName = 'M1';    end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition);    badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion);    badTrialVersion = 'v8'; end
if ~exist('pairedDataFlag','var') || isempty(pairedDataFlag);      pairedDataFlag = 0;     end
if ~exist('saveOutputFlag','var') || isempty(saveOutputFlag);      saveOutputFlag = 1;     end
if ~exist('numPermutations','var') || isempty(numPermutations);    numPermutations = 5000; end
if ~exist('bandInfo','var') || isempty(bandInfo);                  bandInfo = getDefaultBands(); end
if ~exist('conditionType','var') || isempty(conditionType);        conditionType = 'post'; end
if ~exist('densityThreshold','var') || isempty(densityThreshold);  densityThreshold = 0.20; end
if ~exist('gammaValue','var') || isempty(gammaValue);              gammaValue = 1;         end
if ~exist('numCommunityRuns','var') || isempty(numCommunityRuns);  numCommunityRuns = 50;  end

conditionType = validateConditionType(conditionType);
bandInfo = validateBandInfo(bandInfo);
bandTag = getBandInfoTag(bandInfo);
validateDensityThreshold(densityThreshold);

projectRoot = getProjectRoot();
grangerFolder = fullfile(fileparts(mfilename('fullpath')),'savedDataGranger');
summaryFolder = fullfile(grangerFolder,'summary');
statsFolder = fullfile(grangerFolder,'communityStats');
figureFolder = fullfile(statsFolder,'figures');

ensureFolder(summaryFolder);
ensureFolder(statsFolder);
ensureFolder(figureFolder);

ensureBCTOnPath(projectRoot);

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
    error('No valid subjects available for %s community comparison.',comparisonLabel);
end

selectedResultIndices = unique([medIndices(:); ctrlIndices(:)]','stable');

[allMontageChanlocs,allElectrodeLabels] = getMontageChanlocs(projectRoot);
electrodeLabels = allElectrodeLabels;
montageChanlocs = allMontageChanlocs;

localIndexMap = zeros(numel(results.subjectNameList),1);
localIndexMap(selectedResultIndices) = 1:numel(selectedResultIndices);
medLocal = localIndexMap(medIndices);
ctrlLocal = localIndexMap(ctrlIndices);

stats = struct;
stats.protocolName = protocolName;
stats.badEyeCondition = badEyeCondition;
stats.badTrialVersion = badTrialVersion;
stats.bandTag = bandTag;
stats.conditionType = conditionType;
stats.conditionLabel = getConditionLabel(conditionType);
stats.pairedDataFlag = pairedDataFlag;
stats.comparisonLabel = comparisonLabel;
stats.numPermutations = numPermutations;
stats.densityThreshold = densityThreshold;
stats.gammaValue = gammaValue;
stats.numCommunityRuns = numCommunityRuns;
stats.summaryFileName = summaryFileName;
stats.outputFileName = fullfile(statsFolder,[protocolName '_' badEyeCondition '_' badTrialVersion '_' bandTag '_' comparisonLabel '_' conditionType '_grangerCommunityStats.mat']);
stats.commonElectrodes = 1:numel(electrodeLabels);
stats.commonElectrodeLabels = electrodeLabels;
stats.meditatorIndices = medIndices(:);
stats.controlIndices = ctrlIndices(:);
stats.meditatorSubjects = results.subjectNameList(medIndices);
stats.controlSubjects = results.subjectNameList(ctrlIndices);
stats.numMeditators = numel(medIndices);
stats.numControls = numel(ctrlIndices);
stats.pairNames = pairNames;
stats.bandInfo = results.bandInfo;
stats.selectedSubjects = results.subjectNameList(selectedResultIndices);
stats.selectedGroups = results.groupNameList(selectedResultIndices);

globalSpecs = getGlobalMetricSpecs();
nodeSpecs = getNodeMetricSpecs();
numBands = numel(results.bandData);

for iBand = 1:numBands
    bandStats = struct;
    bandStats.name = results.bandData(iBand).name;
    bandStats.range = results.bandData(iBand).range;
    bandStats.conditionType = conditionType;
    bandStats.conditionLabel = stats.conditionLabel;
    bandStats.densityThreshold = densityThreshold;
    bandStats.gammaValue = gammaValue;
    bandStats.numCommunityRuns = numCommunityRuns;
    bandStats.figureFileName = fullfile(figureFolder,[protocolName '_' results.bandData(iBand).name '_' bandTag '_' comparisonLabel '_' conditionType '_community.png']);
    bandStats.globalSpecs = globalSpecs;
    bandStats.nodeSpecs = nodeSpecs;

    bandResults = analyzeBandCommunities(results.bandData(iBand),results,selectedResultIndices,densityThreshold,gammaValue,numCommunityRuns);
    bandStats.analysis = bandResults;

    medPositions = medLocal(:);
    ctrlPositions = ctrlLocal(:);

    bandStats.meditatorMeanMatrix = mean(selectConditionMatrix(bandResults.meanMatrices,medPositions,conditionType),3,'omitnan');
    bandStats.controlMeanMatrix = mean(selectConditionMatrix(bandResults.meanMatrices,ctrlPositions,conditionType),3,'omitnan');
    bandStats.groupDifferenceMatrix = bandStats.meditatorMeanMatrix - bandStats.controlMeanMatrix;

    for iMetric = 1:numel(nodeSpecs)
        metricName = nodeSpecs(iMetric).name;
        metricValues = bandResults.node.(metricName).(conditionType);
        medValues = metricValues(:,medPositions);
        ctrlValues = metricValues(:,ctrlPositions);

        metricStats = struct;
        metricStats.name = metricName;
        metricStats.label = nodeSpecs(iMetric).label;
        metricStats.description = nodeSpecs(iMetric).description;
        metricStats.meditatorValues = medValues;
        metricStats.controlValues = ctrlValues;
        metricStats.meditatorMean = mean(medValues,2,'omitnan');
        metricStats.controlMean = mean(ctrlValues,2,'omitnan');
        metricStats.groupDifference = metricStats.meditatorMean - metricStats.controlMean;
        metricStats.pValue = nan(numel(electrodeLabels),1);
        metricStats.effectSize = nan(numel(electrodeLabels),1);
        metricStats.sampleCount = nan(numel(electrodeLabels),1);

        for iNode = 1:numel(electrodeLabels)
            x = medValues(iNode,:)';
            y = ctrlValues(iNode,:)';
            if pairedDataFlag
                [pVal,effectSize,nUsed] = pairedPermutationTest(x,y,numPermutations);
            else
                [pVal,effectSize,nUsed] = unpairedPermutationTest(x,y,numPermutations);
            end
            metricStats.pValue(iNode) = pVal;
            metricStats.effectSize(iNode) = effectSize;
            metricStats.sampleCount(iNode) = nUsed;
        end

        [metricStats.qValue,metricStats.significantMask] = bhFDR(metricStats.pValue,0.05);
        metricStats.significantElectrodes = getSignificantElectrodes(metricStats,electrodeLabels);
        bandStats.nodeStats.(metricName) = metricStats;
    end

    for iMetric = 1:numel(globalSpecs)
        metricName = globalSpecs(iMetric).name;
        globalValues = bandResults.global.(metricName).(conditionType);
        x = globalValues(medPositions);
        y = globalValues(ctrlPositions);

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

    bandStats.partitionSummary = summarizePartitionMatrices(bandResults.partitionSimilarity,medPositions,ctrlPositions,conditionType);
    stats.bandStats(iBand) = bandStats;

    if saveOutputFlag
        addEEGLABTopoplot(projectRoot);
        hFig = createCommunityFigure(bandStats,montageChanlocs,electrodeLabels,comparisonLabel,stats.numMeditators,stats.numControls,pairedDataFlag);
        saveas(hFig,bandStats.figureFileName);
        close(hFig);
    end
end

if saveOutputFlag
    save(stats.outputFileName,'stats','-v7.3');
end

disp(sprintf('Saved %s community stats for %s using %d meditators and %d controls (%s).',comparisonLabel,protocolName,stats.numMeditators,stats.numControls,stats.conditionLabel)); %#ok<DSPS>
end

function bandResults = analyzeBandCommunities(bandData,results,selectedResultIndices,densityThreshold,gammaValue,numCommunityRuns)
numSubjects = numel(selectedResultIndices);
numNodes = size(bandData.channelPre,1);
globalSpecs = getGlobalMetricSpecs();
nodeSpecs = getNodeMetricSpecs();

bandResults.subjectResultIndices = selectedResultIndices(:)';
bandResults.subjectNames = results.subjectNameList(selectedResultIndices);
bandResults.groupNames = results.groupNameList(selectedResultIndices);
bandResults.commonElectrodes = 1:numNodes;
bandResults.statusPre = cell(numSubjects,1);
bandResults.statusPost = cell(numSubjects,1);
bandResults.availableElectrodes = cell(numSubjects,1);
bandResults.numNodes = numNodes;

bandResults.meanMatrices.pre = nan(numNodes,numNodes,numSubjects);
bandResults.meanMatrices.post = nan(numNodes,numNodes,numSubjects);
bandResults.meanMatrices.diff = nan(numNodes,numNodes,numSubjects);

for iMetric = 1:numel(globalSpecs)
    metricName = globalSpecs(iMetric).name;
    bandResults.global.(metricName).pre = nan(numSubjects,1);
    bandResults.global.(metricName).post = nan(numSubjects,1);
    bandResults.global.(metricName).diff = nan(numSubjects,1);
end

for iMetric = 1:numel(nodeSpecs)
    metricName = nodeSpecs(iMetric).name;
    bandResults.node.(metricName).pre = nan(numNodes,numSubjects);
    bandResults.node.(metricName).post = nan(numNodes,numSubjects);
    bandResults.node.(metricName).diff = nan(numNodes,numSubjects);
end

bandResults.partitionLabels.pre = cell(numSubjects,1);
bandResults.partitionLabels.post = cell(numSubjects,1);
bandResults.partitionSimilarity.pre = nan(numSubjects,numSubjects);
bandResults.partitionSimilarity.post = nan(numSubjects,numSubjects);
bandResults.subjectMeanNMI.pre = nan(numSubjects,1);
bandResults.subjectMeanNMI.post = nan(numSubjects,1);
bandResults.subjectMeanNMI.diff = nan(numSubjects,1);

for iSubject = 1:numSubjects
    subjectIndex = selectedResultIndices(iSubject);
    availableElectrodes = results.selectedElectrodeList{subjectIndex}(:)';
    availableElectrodes = availableElectrodes(isfinite(availableElectrodes));
    availableElectrodes = unique(availableElectrodes,'stable');
    availableElectrodes = availableElectrodes(availableElectrodes >= 1 & availableElectrodes <= numNodes);
    bandResults.availableElectrodes{iSubject} = availableElectrodes;

    if numel(availableElectrodes) < 8
        bandResults.statusPre{iSubject} = 'too_few_subject_electrodes';
        bandResults.statusPost{iSubject} = 'too_few_subject_electrodes';
        continue;
    end

    Gpre = bandData.channelPre(availableElectrodes,availableElectrodes,subjectIndex);
    Gpost = bandData.channelPost(availableElectrodes,availableElectrodes,subjectIndex);

    [subjectPre,statusPre] = analyzeSingleSubjectCommunityGraph(Gpre,densityThreshold,gammaValue,numCommunityRuns);
    [subjectPost,statusPost] = analyzeSingleSubjectCommunityGraph(Gpost,densityThreshold,gammaValue,numCommunityRuns);

    bandResults.statusPre{iSubject} = statusPre;
    bandResults.statusPost{iSubject} = statusPost;

    bandResults.meanMatrices.pre(:,:,iSubject) = embedSquareMatrix(subjectPre.thresholdedMatrix,availableElectrodes,numNodes);
    bandResults.meanMatrices.post(:,:,iSubject) = embedSquareMatrix(subjectPost.thresholdedMatrix,availableElectrodes,numNodes);
    bandResults.meanMatrices.diff(:,:,iSubject) = bandResults.meanMatrices.post(:,:,iSubject) - bandResults.meanMatrices.pre(:,:,iSubject);

    bandResults.partitionLabels.pre{iSubject} = embedVector(subjectPre.communityLabels,availableElectrodes,numNodes);
    bandResults.partitionLabels.post{iSubject} = embedVector(subjectPost.communityLabels,availableElectrodes,numNodes);

    for iMetric = 1:numel(globalSpecs)
        metricName = globalSpecs(iMetric).name;
        bandResults.global.(metricName).pre(iSubject) = subjectPre.global.(metricName);
        bandResults.global.(metricName).post(iSubject) = subjectPost.global.(metricName);
        bandResults.global.(metricName).diff(iSubject) = subjectPost.global.(metricName) - subjectPre.global.(metricName);
    end

    for iMetric = 1:numel(nodeSpecs)
        metricName = nodeSpecs(iMetric).name;
        bandResults.node.(metricName).pre(:,iSubject) = embedVector(subjectPre.node.(metricName),availableElectrodes,numNodes);
        bandResults.node.(metricName).post(:,iSubject) = embedVector(subjectPost.node.(metricName),availableElectrodes,numNodes);
        bandResults.node.(metricName).diff(:,iSubject) = bandResults.node.(metricName).post(:,iSubject) - bandResults.node.(metricName).pre(:,iSubject);
    end
end

bandResults.partitionSimilarity.pre = computePartitionSimilarityMatrix(bandResults.partitionLabels.pre);
bandResults.partitionSimilarity.post = computePartitionSimilarityMatrix(bandResults.partitionLabels.post);
bandResults.subjectMeanNMI.pre = computeSubjectMeanWithinGroupNMI(bandResults.partitionSimilarity.pre,bandResults.groupNames);
bandResults.subjectMeanNMI.post = computeSubjectMeanWithinGroupNMI(bandResults.partitionSimilarity.post,bandResults.groupNames);
bandResults.subjectMeanNMI.diff = bandResults.subjectMeanNMI.post - bandResults.subjectMeanNMI.pre;

bandResults.global.meanWithinGroupNMI.pre = bandResults.subjectMeanNMI.pre;
bandResults.global.meanWithinGroupNMI.post = bandResults.subjectMeanNMI.post;
bandResults.global.meanWithinGroupNMI.diff = bandResults.subjectMeanNMI.diff;
end

function [subjectMetrics,statusString] = analyzeSingleSubjectCommunityGraph(G,densityThreshold,gammaValue,numCommunityRuns)
statusString = 'ok';
W = thresholdDirectedGraphDensity(G,densityThreshold);
numPositiveEdges = nnz(W > 0);

subjectMetrics = initializeSubjectMetricStruct(size(W,1));
subjectMetrics.thresholdedMatrix = W;
subjectMetrics.numPositiveEdges = numPositiveEdges;
subjectMetrics.edgeDensityAchieved = numPositiveEdges / max(1,(numel(W) - size(W,1)));

if numPositiveEdges < 2
    statusString = 'too_few_edges_after_threshold';
    return;
end

[communityLabels,modularityQ] = detectBestDirectedCommunities(W,gammaValue,numCommunityRuns);
if isempty(communityLabels) || ~all(isfinite(communityLabels))
    statusString = 'community_detection_failed';
    return;
end

participationOut = participation_coef(W,communityLabels,1);
participationIn = participation_coef(W,communityLabels,2);
communityCounts = accumarray(communityLabels(:),1);

subjectMetrics.communityLabels = communityLabels(:);
subjectMetrics.node.participationOut = participationOut(:);
subjectMetrics.node.participationIn = participationIn(:);
subjectMetrics.global.numCommunities = numel(unique(communityLabels));
subjectMetrics.global.modularityQ = modularityQ;
subjectMetrics.global.meanParticipationOut = mean(participationOut,'omitnan');
subjectMetrics.global.meanParticipationIn = mean(participationIn,'omitnan');
subjectMetrics.global.maxCommunityFraction = max(communityCounts) / numel(communityLabels);
end

function subjectMetrics = initializeSubjectMetricStruct(numNodes)
subjectMetrics = struct;
subjectMetrics.thresholdedMatrix = zeros(numNodes,numNodes);
subjectMetrics.communityLabels = nan(numNodes,1);
subjectMetrics.numPositiveEdges = 0;
subjectMetrics.edgeDensityAchieved = NaN;
subjectMetrics.node.participationOut = nan(numNodes,1);
subjectMetrics.node.participationIn = nan(numNodes,1);
subjectMetrics.global.numCommunities = NaN;
subjectMetrics.global.modularityQ = NaN;
subjectMetrics.global.meanParticipationOut = NaN;
subjectMetrics.global.meanParticipationIn = NaN;
subjectMetrics.global.maxCommunityFraction = NaN;
subjectMetrics.global.meanWithinGroupNMI = NaN;
end

function W = thresholdDirectedGraphDensity(G,densityThreshold)
G = double(G);
G = removeSelfConnections(G);
G(~isfinite(G)) = 0;
G(G < 0) = 0;

offDiagMask = ~eye(size(G,1));
candidateMask = offDiagMask & G > 0;
candidateValues = G(candidateMask);

W = zeros(size(G));
if isempty(candidateValues)
    return;
end

nKeep = max(1,round(densityThreshold * nnz(offDiagMask)));
nKeep = min(nKeep,numel(candidateValues));

[sortedValues,sortedOrder] = sort(candidateValues,'descend');
keepLogical = false(size(candidateValues));
keepLogical(sortedOrder(1:nKeep)) = true;

linearIndices = find(candidateMask);
W(linearIndices(keepLogical)) = candidateValues(keepLogical);
end

function [bestCommunityLabels,bestQ] = detectBestDirectedCommunities(W,gammaValue,numCommunityRuns)
bestQ = -inf;
bestCommunityLabels = [];

for iRun = 1:numCommunityRuns
    try
        [communityLabels,Q] = community_louvain(W,gammaValue,[],'modularity');
    catch
        try
            [communityLabels,Q] = modularity_dir(W,gammaValue);
        catch
            communityLabels = [];
            Q = NaN;
        end
    end

    if isempty(communityLabels) || ~isfinite(Q)
        continue;
    end

    if Q > bestQ
        bestQ = Q;
        bestCommunityLabels = communityLabels(:);
    end
end

if isempty(bestCommunityLabels)
    bestCommunityLabels = nan(size(W,1),1);
    bestQ = NaN;
end
end

function partitionSimilarity = computePartitionSimilarityMatrix(partitionLabels)
numSubjects = numel(partitionLabels);
partitionSimilarity = nan(numSubjects,numSubjects);

for iSubject = 1:numSubjects
    partitionSimilarity(iSubject,iSubject) = 1;
    for jSubject = iSubject+1:numSubjects
        CiFull = partitionLabels{iSubject};
        CjFull = partitionLabels{jSubject};
        if isempty(CiFull) || isempty(CjFull)
            continue;
        end

        overlapMask = isfinite(CiFull) & isfinite(CjFull);
        if nnz(overlapMask) < 8
            continue;
        end

        Ci = CiFull(overlapMask);
        Cj = CjFull(overlapMask);
        [~,nmiVal] = partition_distance(Ci,Cj);
        partitionSimilarity(iSubject,jSubject) = nmiVal;
        partitionSimilarity(jSubject,iSubject) = nmiVal;
    end
end
end

function fullMatrix = embedSquareMatrix(localMatrix,availableElectrodes,numNodes)
fullMatrix = nan(numNodes,numNodes);
if isempty(localMatrix) || isempty(availableElectrodes)
    return;
end
fullMatrix(availableElectrodes,availableElectrodes) = localMatrix;
end

function fullVector = embedVector(localVector,availableElectrodes,numNodes)
fullVector = nan(numNodes,1);
if isempty(localVector) || isempty(availableElectrodes)
    return;
end
fullVector(availableElectrodes) = localVector(:);
end

function subjectMeanNMI = computeSubjectMeanWithinGroupNMI(partitionSimilarity,groupNames)
numSubjects = size(partitionSimilarity,1);
subjectMeanNMI = nan(numSubjects,1);

for iSubject = 1:numSubjects
    sameGroupMask = strcmp(groupNames,groupNames{iSubject});
    sameGroupMask(iSubject) = false;
    vals = partitionSimilarity(iSubject,sameGroupMask);
    vals = vals(isfinite(vals));
    if ~isempty(vals)
        subjectMeanNMI(iSubject) = mean(vals,'omitnan');
    end
end
end

function matrixOut = selectConditionMatrix(meanMatrices,indices,conditionType)
switch conditionType
    case 'pre'
        matrixOut = meanMatrices.pre(:,:,indices);
    case 'post'
        matrixOut = meanMatrices.post(:,:,indices);
    case 'diff'
        matrixOut = meanMatrices.diff(:,:,indices);
    otherwise
        error('Unknown conditionType: %s',conditionType);
end
end

function partitionSummary = summarizePartitionMatrices(partitionSimilarity,medPositions,ctrlPositions,conditionType)
switch conditionType
    case 'pre'
        matrixToUse = partitionSimilarity.pre;
    case 'post'
        matrixToUse = partitionSimilarity.post;
    case 'diff'
        matrixToUse = partitionSimilarity.post - partitionSimilarity.pre;
    otherwise
        error('Unknown conditionType: %s',conditionType);
end

partitionSummary = struct;
partitionSummary.conditionType = conditionType;
partitionSummary.meditatorMeanNMI = meanFiniteUpperTriangle(getSubMatrix(matrixToUse,medPositions,medPositions));
partitionSummary.controlMeanNMI = meanFiniteUpperTriangle(getSubMatrix(matrixToUse,ctrlPositions,ctrlPositions));
partitionSummary.crossGroupMeanNMI = meanFiniteValues(getSubMatrix(matrixToUse,medPositions,ctrlPositions));
end

function subMatrix = getSubMatrix(M,rowIdx,colIdx)
subMatrix = M(rowIdx,colIdx);
end

function meanVal = meanFiniteUpperTriangle(M)
if isempty(M)
    meanVal = NaN;
    return;
end

mask = triu(true(size(M)),1);
vals = M(mask);
vals = vals(isfinite(vals));
if isempty(vals)
    meanVal = NaN;
else
    meanVal = mean(vals,'omitnan');
end
end

function meanVal = meanFiniteValues(M)
vals = M(isfinite(M));
if isempty(vals)
    meanVal = NaN;
else
    meanVal = mean(vals,'omitnan');
end
end

function hFig = createCommunityFigure(bandStats,montageChanlocs,electrodeLabels,comparisonLabel,numMeditators,numControls,pairedDataFlag)
hFig = figure('Color','w','Name',[bandStats.name ' communities']);
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
plotCommunitySummaryText(gca,bandStats,comparisonLabel,numMeditators,numControls);

subplot(3,4,5);
plotTopologyTopoplot(gca,bandStats.nodeStats.participationOut,montageChanlocs,symmetricCLim(bandStats.nodeStats.participationOut.groupDifference));

subplot(3,4,6);
plotTopologyTopoplot(gca,bandStats.nodeStats.participationIn,montageChanlocs,symmetricCLim(bandStats.nodeStats.participationIn.groupDifference));

plotGlobalMetricDistribution(gcaForSubplot(3,4,7),bandStats.globalStats(1),pairedDataFlag);
plotGlobalMetricDistribution(gcaForSubplot(3,4,8),bandStats.globalStats(2),pairedDataFlag);
plotGlobalMetricDistribution(gcaForSubplot(3,4,9),bandStats.globalStats(3),pairedDataFlag);
plotGlobalMetricDistribution(gcaForSubplot(3,4,10),bandStats.globalStats(4),pairedDataFlag);
plotGlobalMetricDistribution(gcaForSubplot(3,4,11),bandStats.globalStats(5),pairedDataFlag);

subplot(3,4,12);
plotPartitionSummaryPanel(gca,bandStats);

if exist('sgtitle','file') == 2
    sgtitle(sprintf('%s band (%g-%g Hz) | %s | %s communities | density %.0f%%', ...
        bandStats.name,bandStats.range(1),bandStats.range(2),comparisonLabel,bandStats.conditionLabel,100*bandStats.densityThreshold),'Interpreter','none');
end
end

function plotCommunitySummaryText(ax,bandStats,comparisonLabel,numMeditators,numControls)
axis(ax,'off');

summaryLines = {
    sprintf('Comparison: %s',comparisonLabel)
    sprintf('Condition: %s',bandStats.conditionLabel)
    sprintf('Meditators: %d',numMeditators)
    sprintf('Controls: %d',numControls)
    sprintf('Montage electrodes: %d',bandStats.analysis.numNodes)
    'Subject electrodes vary by bad-channel mask'
    sprintf('Density kept: %.0f%%',100*bandStats.densityThreshold)
    sprintf('Gamma: %.2f',bandStats.gammaValue)
    sprintf('Community runs/subject: %d',bandStats.numCommunityRuns)
    ' '
    'Global metrics below:'
    'k, Q, mean Pout, mean Pin, max module frac'
    ' '
    sprintf('Within-group mean NMI (%s):',bandStats.conditionLabel)
    sprintf('Meditator: %.3f',bandStats.partitionSummary.meditatorMeanNMI)
    sprintf('Control: %.3f',bandStats.partitionSummary.controlMeanNMI)
    sprintf('Cross-group: %.3f',bandStats.partitionSummary.crossGroupMeanNMI)
    };

summaryText = sprintf('%s\n',summaryLines{:});
text(0,1,summaryText,'Parent',ax,'VerticalAlignment','top','FontName','Helvetica','FontSize',11);
end

function plotPartitionSummaryPanel(ax,bandStats)
axis(ax,'off');
metricNames = {bandStats.globalStats.name};
stabilityIdx = find(strcmp(metricNames,'meanWithinGroupNMI'),1);

if isempty(stabilityIdx)
    return;
end

metricStats = bandStats.globalStats(stabilityIdx);
summaryLines = {
    sprintf('%s',metricStats.label)
    sprintf('Meditator mean: %.3f',metricStats.meditatorMean)
    sprintf('Control mean: %.3f',metricStats.controlMean)
    sprintf('Group diff: %.3f',metricStats.groupDifference)
    sprintf('p = %.4f',metricStats.pValue)
    sprintf('q = %.4f',metricStats.qValue)
    sprintf('Effect size: %.3f',metricStats.effectSize)
    };
text(0,1,sprintf('%s\n',summaryLines{:}),'Parent',ax,'VerticalAlignment','top','FontName','Helvetica','FontSize',11);
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
        theta = 0;
        radius = 0;
        if isfield(montageChanlocs,'theta'); theta = montageChanlocs(idx).theta; end
        if isfield(montageChanlocs,'radius'); radius = montageChanlocs(idx).radius; end
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

function cLim = symmetricCLim(x)
vals = x(isfinite(x));
if isempty(vals)
    cLim = [];
    return;
end
maxAbs = max(abs(vals));
if maxAbs == 0
    cLim = [-1 1];
else
    cLim = [-maxAbs maxAbs];
end
end

function val = maxFiniteValue(x)
vals = x(isfinite(x));
if isempty(vals)
    val = 0;
else
    val = max(vals);
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

function nodeSpecs = getNodeMetricSpecs()
nodeSpecs(1).name = 'participationOut';
nodeSpecs(1).label = 'Participation Out';
nodeSpecs(1).description = 'Directed out-participation coefficient';

nodeSpecs(2).name = 'participationIn';
nodeSpecs(2).label = 'Participation In';
nodeSpecs(2).description = 'Directed in-participation coefficient';
end

function globalSpecs = getGlobalMetricSpecs()
globalSpecs(1).name = 'numCommunities';
globalSpecs(1).label = 'Number of communities';
globalSpecs(1).description = 'Count of detected modules';

globalSpecs(2).name = 'modularityQ';
globalSpecs(2).label = 'Modularity Q';
globalSpecs(2).description = 'Directed weighted modularity';

globalSpecs(3).name = 'meanParticipationOut';
globalSpecs(3).label = 'Mean participation out';
globalSpecs(3).description = 'Average out-participation across electrodes';

globalSpecs(4).name = 'meanParticipationIn';
globalSpecs(4).label = 'Mean participation in';
globalSpecs(4).description = 'Average in-participation across electrodes';

globalSpecs(5).name = 'maxCommunityFraction';
globalSpecs(5).label = 'Largest community fraction';
globalSpecs(5).description = 'Fraction of electrodes in the largest module';

globalSpecs(6).name = 'meanWithinGroupNMI';
globalSpecs(6).label = 'Mean within-group NMI';
globalSpecs(6).description = 'Subject-wise mean partition similarity within group';
end

function [medIndices,ctrlIndices] = getValidUnpairedIndices(results)
validMask = results.validSubjectMask(:) & ~cellfun(@isempty,results.selectedElectrodeList(:));
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

    if isempty(results.selectedElectrodeList{medIndex}) || isempty(results.selectedElectrodeList{ctrlIndex})
        continue;
    end

    medIndices(end+1,1) = medIndex; %#ok<AGROW>
    ctrlIndices(end+1,1) = ctrlIndex; %#ok<AGROW>
    pairNames{end+1,1} = [medName '_' ctrlName]; %#ok<AGROW>
end
end

function conditionType = validateConditionType(conditionType)
conditionType = lower(char(conditionType));
validConditionList = {'pre','post','diff'};
if ~ismember(conditionType,validConditionList)
    error('conditionType must be one of: %s',strjoin(validConditionList,', '));
end
end

function validateDensityThreshold(densityThreshold)
if ~isscalar(densityThreshold) || ~isfinite(densityThreshold) || densityThreshold <= 0 || densityThreshold > 1
    error('densityThreshold must be a scalar in the interval (0, 1].');
end
end

function label = getConditionLabel(conditionType)
switch conditionType
    case 'pre'
        label = 'Pre';
    case 'post'
        label = 'Post';
    case 'diff'
        label = 'Post - Pre';
    otherwise
        label = conditionType;
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

function bandTag = getBandInfoTag(bandInfo)
tagParts = cell(1,numel(bandInfo));
for iBand = 1:numel(bandInfo)
    safeName = regexprep(lower(bandInfo(iBand).name),'[^a-z0-9]+','');
    rangeVals = round(100*bandInfo(iBand).range);
    tagParts{iBand} = sprintf('%s_%g_%g',safeName,rangeVals(1),rangeVals(2));
end
bandTag = strjoin(tagParts,'__');
end

function projectRoot = getProjectRoot()
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

function ensureBCTOnPath(projectRoot)
fieldTripPath = fullfile(projectRoot,'fieldtrip-20260211');
bctPath = fullfile(fieldTripPath,'external','bct');
if exist('community_louvain','file') ~= 2
    addpath(fieldTripPath);
    if exist('ft_defaults','file') == 2
        ft_defaults;
    end
    addpath(bctPath);
end
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

function G = removeSelfConnections(G)
n = size(G,1);
G(1:n+1:end) = 0;
end

function ax = gcaForSubplot(m,n,p)
ax = subplot(m,n,p);
end

function ensureFolder(folderName)
if ~exist(folderName,'dir')
    mkdir(folderName);
end
end
