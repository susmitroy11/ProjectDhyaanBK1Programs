function stats = compareGrangerMeditatorsControlsBK1(protocolName,badEyeCondition,badTrialVersion,pairedDataFlag,saveOutputFlag,numPermutations,bandInfo)
% compareGrangerMeditatorsControlsBK1 Compare BK1 ROI-level Granger by group.
%
% Example:
%   stats = compareGrangerMeditatorsControlsBK1('M1');
%   stats = compareGrangerMeditatorsControlsBK1('M1','ep','v8',1);
%   bandInfo(1).name = 'custom1'; bandInfo(1).range = [25 45];
%   stats = compareGrangerMeditatorsControlsBK1('M1','ep','v8',0,1,5000,bandInfo);
%
% Notes:
%   This script works on ROI-level pre, post, and post-pre Granger matrices
%   from extractGrangerBandDataAllSubjectsBK1 and saves figures/statistics
%   in savedDataGranger/groupStats.

if ~exist('protocolName','var') || isempty(protocolName);       protocolName = 'M1';    end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition); badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion); badTrialVersion = 'v8'; end
if ~exist('pairedDataFlag','var') || isempty(pairedDataFlag);   pairedDataFlag = 0;     end
if ~exist('saveOutputFlag','var') || isempty(saveOutputFlag);   saveOutputFlag = 1;     end
if ~exist('numPermutations','var') || isempty(numPermutations); numPermutations = 5000; end
if ~exist('bandInfo','var') || isempty(bandInfo);               bandInfo = getDefaultBands(); end

grangerFolder = fullfile(fileparts(mfilename('fullpath')),'savedDataGranger');
summaryFolder = fullfile(grangerFolder,'summary');
statsFolder = fullfile(grangerFolder,'groupStats');
figureFolder = fullfile(statsFolder,'figures');

ensureFolder(summaryFolder);
ensureFolder(statsFolder);
ensureFolder(figureFolder);

bandInfo = validateBandInfo(bandInfo);
bandTag = getBandInfoTag(bandInfo);

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

conditionSpecs = getConditionSpecs();
numBands = numel(results.bandData);
numROIs = numel(results.roiLabels);

stats = struct;
stats.protocolName = protocolName;
stats.badEyeCondition = badEyeCondition;
stats.badTrialVersion = badTrialVersion;
stats.pairedDataFlag = pairedDataFlag;
stats.comparisonLabel = comparisonLabel;
stats.numPermutations = numPermutations;
stats.summaryFileName = summaryFileName;
stats.bandTag = bandTag;
stats.outputFileName = fullfile(statsFolder,[protocolName '_' badEyeCondition '_' badTrialVersion '_' bandTag '_' comparisonLabel '_grangerGroupStats.mat']);
stats.meditatorIndices = medIndices(:);
stats.controlIndices = ctrlIndices(:);
stats.meditatorSubjects = results.subjectNameList(medIndices);
stats.controlSubjects = results.subjectNameList(ctrlIndices);
stats.pairNames = pairNames;
stats.roiLabels = results.roiLabels;
stats.bandInfo = results.bandInfo;
stats.conditionSpecs = conditionSpecs;
stats.numMeditators = numel(medIndices);
stats.numControls = numel(ctrlIndices);

for iBand = 1:numBands
    stats.bandStats(iBand).name = results.bandData(iBand).name;
    stats.bandStats(iBand).range = results.bandData(iBand).range;
    stats.bandStats(iBand).figureFileName = fullfile(figureFolder,[protocolName '_' results.bandData(iBand).name '_' bandTag '_' comparisonLabel '_groupStats.png']);

    for iCondition = 1:numel(conditionSpecs)
        conditionStats = computeConditionStats(results.bandData(iBand),conditionSpecs(iCondition),medIndices,ctrlIndices,pairedDataFlag,numPermutations,numROIs,results.roiLabels);
        stats.bandStats(iBand).(conditionSpecs(iCondition).name) = conditionStats;
    end

    if saveOutputFlag
        hFig = createBandFigure(stats.bandStats(iBand),results.roiLabels,comparisonLabel,stats.numMeditators,stats.numControls,conditionSpecs);
        saveas(hFig,stats.bandStats(iBand).figureFileName);
        close(hFig);
    end
end

if saveOutputFlag
    save(stats.outputFileName,'stats','-v7.3');
end

disp(sprintf('Saved %s group stats for %s using %d meditators and %d controls.',comparisonLabel,protocolName,stats.numMeditators,stats.numControls)); %#ok<DSPS>
end

function conditionStats = computeConditionStats(bandData,conditionSpec,medIndices,ctrlIndices,pairedDataFlag,numPermutations,numROIs,roiLabels)
roiData = bandData.(conditionSpec.roiField);
netData = bandData.(conditionSpec.netField);

roiMed = roiData(:,:,medIndices);
roiCtrl = roiData(:,:,ctrlIndices);

conditionStats.name = conditionSpec.name;
conditionStats.label = conditionSpec.label;
conditionStats.roiField = conditionSpec.roiField;
conditionStats.netField = conditionSpec.netField;
conditionStats.meditatorMeanROI = mean(roiMed,3,'omitnan');
conditionStats.controlMeanROI = mean(roiCtrl,3,'omitnan');
conditionStats.groupDifferenceROI = conditionStats.meditatorMeanROI - conditionStats.controlMeanROI;
conditionStats.roiSampleCount = nan(numROIs,numROIs);
conditionStats.roiPValue = nan(numROIs,numROIs);
conditionStats.roiQValue = nan(numROIs,numROIs);
conditionStats.roiEffectSize = nan(numROIs,numROIs);

for iROI = 1:numROIs
    for jROI = 1:numROIs
        x = squeeze(roiMed(iROI,jROI,:));
        y = squeeze(roiCtrl(iROI,jROI,:));

        if pairedDataFlag
            [pVal,effectSize,nUsed] = pairedPermutationTest(x,y,numPermutations);
        else
            [pVal,effectSize,nUsed] = unpairedPermutationTest(x,y,numPermutations);
        end

        conditionStats.roiPValue(iROI,jROI) = pVal;
        conditionStats.roiEffectSize(iROI,jROI) = effectSize;
        conditionStats.roiSampleCount(iROI,jROI) = nUsed;
    end
end

[conditionStats.roiQValue,conditionStats.roiSignificantMask] = bhFDR(conditionStats.roiPValue,0.05);
conditionStats.significantPairs = getSignificantPairs(conditionStats,roiLabels);

medNet = netData(:,medIndices);
ctrlNet = netData(:,ctrlIndices);
conditionStats.meditatorMeanNet = mean(medNet,2,'omitnan');
conditionStats.controlMeanNet = mean(ctrlNet,2,'omitnan');
conditionStats.netDifference = conditionStats.meditatorMeanNet - conditionStats.controlMeanNet;
conditionStats.netSampleCount = nan(numROIs,1);
conditionStats.netPValue = nan(numROIs,1);
conditionStats.netEffectSize = nan(numROIs,1);

for iROI = 1:numROIs
    x = medNet(iROI,:)';
    y = ctrlNet(iROI,:)';
    if pairedDataFlag
        [pVal,effectSize,nUsed] = pairedPermutationTest(x,y,numPermutations);
    else
        [pVal,effectSize,nUsed] = unpairedPermutationTest(x,y,numPermutations);
    end
    conditionStats.netPValue(iROI) = pVal;
    conditionStats.netEffectSize(iROI) = effectSize;
    conditionStats.netSampleCount(iROI) = nUsed;
end

[conditionStats.netQValue,conditionStats.netSignificantMask] = bhFDR(conditionStats.netPValue,0.05);
end

function [medIndices,ctrlIndices] = getValidUnpairedIndices(results)
validMask = results.validSubjectMask(:);
groupList = results.groupNameList(:);

medIndices = find(validMask & strcmp(groupList,'Meditator'));
ctrlIndices = find(validMask & strcmp(groupList,'Control'));
end

function conditionSpecs = getConditionSpecs()
conditionSpecs(1).name = 'pre';
conditionSpecs(1).label = 'Pre';
conditionSpecs(1).roiField = 'roiPre';
conditionSpecs(1).netField = 'roiNetPre';

conditionSpecs(2).name = 'post';
conditionSpecs(2).label = 'Post';
conditionSpecs(2).roiField = 'roiPost';
conditionSpecs(2).netField = 'roiNetPost';

conditionSpecs(3).name = 'diff';
conditionSpecs(3).label = 'Post - Pre';
conditionSpecs(3).roiField = 'roiDiff';
conditionSpecs(3).netField = 'roiNetDiff';
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

function hFig = createBandFigure(bandStats,roiLabels,comparisonLabel,numMeditators,numControls,conditionSpecs)
hFig = figure('Color','w','Name',[bandStats.name ' group stats']);
hFig.Position(3:4) = [1800 1200];

for iCondition = 1:numel(conditionSpecs)
    conditionName = conditionSpecs(iCondition).name;
    cStats = bandStats.(conditionName);
    rowOffset = (iCondition-1) * 4;

    groupCLim = getConditionGroupCLim(cStats,conditionName);
    diffCLim = symmetricCLim(cStats.groupDifferenceROI);

    subplot(3,4,rowOffset+1);
    plotHeatmap(gca,cStats.meditatorMeanROI,roiLabels,groupCLim);
    title(sprintf('%s: Meditator (%d)',cStats.label,numMeditators),'Interpreter','none');

    subplot(3,4,rowOffset+2);
    plotHeatmap(gca,cStats.controlMeanROI,roiLabels,groupCLim);
    title(sprintf('%s: Control (%d)',cStats.label,numControls),'Interpreter','none');

    subplot(3,4,rowOffset+3);
    plotHeatmap(gca,cStats.groupDifferenceROI,roiLabels,diffCLim);
    title(sprintf('%s: Meditator - Control',cStats.label),'Interpreter','none');

    subplot(3,4,rowOffset+4);
    plotConditionNetSummary(gca,cStats,roiLabels);
    title(sprintf('%s: ROI Net',cStats.label),'Interpreter','none');
end

if exist('sgtitle','file') == 2
    sgtitle(sprintf('%s band (%g-%g Hz) | %s comparison',bandStats.name,bandStats.range(1),bandStats.range(2),comparisonLabel),'Interpreter','none');
end
end

function plotConditionNetSummary(ax,conditionStats,roiLabels)
barData = [conditionStats.meditatorMeanNet conditionStats.controlMeanNet];
bar(ax,barData,'grouped');
set(ax,'XTick',1:numel(roiLabels),'XTickLabel',roiLabels,'TickDir','out');
xtickangle(ax,45);
ylabel(ax,'Net Flow');
legend(ax,{'Meditator','Control'},'Location','best');
grid(ax,'on');

sigCountROI = sum(conditionStats.roiSignificantMask(:));
sigCountNet = sum(conditionStats.netSignificantMask(:));
text(0.02,0.98,sprintf('ROI q<0.05: %d\nNet q<0.05: %d',sigCountROI,sigCountNet), ...
    'Units','normalized','VerticalAlignment','top','Parent',ax);
end

function cLim = getConditionGroupCLim(conditionStats,conditionName)
if strcmp(conditionName,'diff')
    cLim = symmetricCLim([conditionStats.meditatorMeanROI(:); conditionStats.controlMeanROI(:)]);
else
    cLim = [0 max([maxFiniteValue(conditionStats.meditatorMeanROI) maxFiniteValue(conditionStats.controlMeanROI)])];
end
end

function plotHeatmap(ax,matrix,tickLabels,cLims)
imagesc(ax,matrix,'AlphaData',isfinite(matrix));
set(ax,'YDir','normal');
axis(ax,'square');
colorbar(ax);
colormap(ax,parula);
set(ax,'XTick',1:numel(tickLabels),'XTickLabel',tickLabels,'YTick',1:numel(tickLabels),'YTickLabel',tickLabels,'TickDir','out');
xtickangle(ax,45);
if ~isempty(cLims) && all(isfinite(cLims))
    caxis(ax,cLims);
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

function significantPairs = getSignificantPairs(conditionStats,roiLabels)
significantPairs = {};
[rowIdx,colIdx] = find(conditionStats.roiSignificantMask);
for iPair = 1:numel(rowIdx)
    significantPairs{iPair,1} = roiLabels{rowIdx(iPair)}; %#ok<AGROW>
    significantPairs{iPair,2} = roiLabels{colIdx(iPair)}; %#ok<AGROW>
    significantPairs{iPair,3} = conditionStats.groupDifferenceROI(rowIdx(iPair),colIdx(iPair)); %#ok<AGROW>
    significantPairs{iPair,4} = conditionStats.roiPValue(rowIdx(iPair),colIdx(iPair)); %#ok<AGROW>
    significantPairs{iPair,5} = conditionStats.roiQValue(rowIdx(iPair),colIdx(iPair)); %#ok<AGROW>
    significantPairs{iPair,6} = conditionStats.roiEffectSize(rowIdx(iPair),colIdx(iPair)); %#ok<AGROW>
end
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

function ensureFolder(folderName)
if ~exist(folderName,'dir')
    mkdir(folderName);
end
end
