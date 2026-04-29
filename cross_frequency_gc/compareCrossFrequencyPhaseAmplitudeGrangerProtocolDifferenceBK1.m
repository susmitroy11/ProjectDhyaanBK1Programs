function analysisResult = compareCrossFrequencyPhaseAmplitudeGrangerProtocolDifferenceBK1(protocolNameA,protocolNameB,pairSelection,badEyeCondition,badTrialVersion,cfgIn)
% compareCrossFrequencyPhaseAmplitudeGrangerProtocolDifferenceBK1
%
% Compare two protocols directly using matched-pair CF-GC outputs.
% This function reuses the multi-pair comparison workflow, but restricts
% the analysis to two named protocols and returns their difference map for
% every ROI and band-pair cell.
%
% Subtraction order is:
%   protocolNameA - protocolNameB
%
% Example:
%   compareCrossFrequencyPhaseAmplitudeGrangerProtocolDifferenceBK1( ...
%       'G2','G1','all','ep','v8')

if ~exist('protocolNameA','var') || isempty(protocolNameA)
    error('protocolNameA is required. Example: ''G2''.');
end
if ~exist('protocolNameB','var') || isempty(protocolNameB)
    error('protocolNameB is required. Example: ''G1''.');
end
if ~exist('pairSelection','var') || isempty(pairSelection);       pairSelection = 'all'; end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition);   badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion);   badTrialVersion = 'v8'; end
if ~exist('cfgIn','var') || isempty(cfgIn);                       cfgIn = struct(); end

cfg = applyDefaultCfg(cfgIn);
protocolNameA = validateSingleProtocolName(protocolNameA);
protocolNameB = validateSingleProtocolName(protocolNameB);
if strcmp(protocolNameA,protocolNameB)
    error('The two protocol names must be different.');
end

thisFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisFolder);
[pairList,selectionLabel] = cfGCUtils.resolvePairSelectionBK1(pairSelection,projectRoot);
familyDefs = getFamilyDefinitions();
differenceLabel = sprintf('%s - %s',protocolNameA,protocolNameB);
differenceToken = sprintf('%s_minus_%s',protocolNameA,protocolNameB);

optionsText = sprintf(['protocol_difference|protocolA=%s|protocolB=%s|pairs=%s|eye=%s|trial=%s|' ...
    'perm=%d|minobs=%d|minpairs=%d|fdr=%.4f|overview=%d|save=%d'], ...
    protocolNameA,protocolNameB,strjoin({pairList.pairLabel},','),badEyeCondition,badTrialVersion, ...
    cfg.numPermutations,cfg.minObservationCount,cfg.minValidSubjectsPerCell,cfg.fdrAlpha, ...
    cfg.makeOverviewFigure,cfg.saveOutputFlag);
analysisTag = cfGCUtils.makeAnalysisTag([selectionLabel '_' differenceToken],optionsText);
analysisRoot = buildAnalysisRoot(cfg.analysisFolder,cfg.analysisSubfolder,analysisTag);
summaryFileName = fullfile(analysisRoot,[differenceToken '_multiple_pairs_protocol_difference_summary.mat']);

[cachedResult,wasLoaded] = cfGCUtils.tryLoadCachedAnalysis(summaryFileName,cfg.saveOutputFlag && ~cfg.forceRebuild);
if wasLoaded && ~cfg.forceRebuild
    analysisResult = cachedResult;
    fprintf('Loaded cached protocol-difference comparison:\n%s\n',summaryFileName);
    return;
end

protocolAData = prepareProtocolArrays(pairList,protocolNameA,badEyeCondition,badTrialVersion,cfg);
protocolBData = prepareProtocolArrays(pairList,protocolNameB,badEyeCondition,badTrialVersion,cfg);

cfGCUtils.validateResultCompatibility(protocolAData.templateResults,protocolBData.templateResults);
commonValidPairMask = protocolAData.validPairMask & protocolBData.validPairMask;
commonIndices = find(commonValidPairMask);
if isempty(commonIndices)
    error(['No matched pairs had complete PRE/POST CF-GC files in both protocols %s and %s ' ...
        'for the requested selection.'],protocolNameA,protocolNameB);
end

analysisResult = struct();
analysisResult.mainOption = 'comparison';
analysisResult.subOption = 'protocol_difference';
analysisResult.protocolNameA = protocolNameA;
analysisResult.protocolNameB = protocolNameB;
analysisResult.differenceLabel = differenceLabel;
analysisResult.badEyeCondition = badEyeCondition;
analysisResult.badTrialVersion = badTrialVersion;
analysisResult.requestedPairList = pairList;
analysisResult.commonPairList = pairList(commonIndices);
analysisResult.protocolAValidPairMask = protocolAData.validPairMask;
analysisResult.protocolBValidPairMask = protocolBData.validPairMask;
analysisResult.commonValidPairMask = commonValidPairMask;
analysisResult.protocolAPairAvailability = protocolAData.pairAvailability;
analysisResult.protocolBPairAvailability = protocolBData.pairAvailability;
analysisResult.numRequestedPairs = numel(pairList);
analysisResult.numCommonPairs = numel(commonIndices);
analysisResult.analysisFolder = analysisRoot;
analysisResult.summaryFileName = summaryFileName;
analysisResult.roiLabels = protocolAData.roiLabels;
analysisResult.pairInfo = protocolAData.pairInfo;
analysisResult.phaseBandList = protocolAData.phaseBandList;
analysisResult.ampBandList = protocolAData.ampBandList;
analysisResult.familyDefinitions = familyDefs;

analysisResult.rawDifference = struct();
analysisResult.normalizedDifference = struct();
analysisResult.groupMeans = struct('raw',struct(),'normalized',struct());

for iFamily = 1:numel(familyDefs)
    familyName = familyDefs(iFamily).name;
    rawA = extractFamilyMetric(protocolAData.raw,familyName);
    rawB = extractFamilyMetric(protocolBData.raw,familyName);
    normA = extractFamilyMetric(protocolAData.normalized,familyName);
    normB = extractFamilyMetric(protocolBData.normalized,familyName);

    analysisResult.rawDifference.(familyName) = rawA(:,:,:,:,commonIndices) - rawB(:,:,:,:,commonIndices);
    analysisResult.normalizedDifference.(familyName) = normA(:,:,:,:,commonIndices) - normB(:,:,:,:,commonIndices);
    analysisResult.groupMeans.raw.(familyName) = mean(analysisResult.rawDifference.(familyName),5,'omitnan');
    analysisResult.groupMeans.normalized.(familyName) = mean(analysisResult.normalizedDifference.(familyName),5,'omitnan');
end

commonSampleCount = min( ...
    protocolAData.sampleCount.minAcrossAll(:,:,:,:,commonIndices), ...
    protocolBData.sampleCount.minAcrossAll(:,:,:,:,commonIndices));
analysisResult.sampleCount.minAcrossProtocols = commonSampleCount;
analysisResult.inference = computeProtocolDifferenceInference(protocolAData,protocolBData,commonIndices,familyDefs,cfg);
analysisResult.summaryOverview = buildSummaryOverview(analysisResult.groupMeans.raw,analysisResult.inference,commonSampleCount,familyDefs);
analysisResult.qValueSummaryTable = makeQSummaryTable(familyDefs,analysisResult.pairInfo,analysisResult.inference);

figureFiles = struct();
outputManifest = {};
if cfg.saveOutputFlag
    cfg.figureFolder = fullfile(analysisRoot,'figures');
    cfGCUtils.ensureFolder(cfg.figureFolder);

    allMeanMaps = cell(1,numel(familyDefs));
    for iFamily = 1:numel(familyDefs)
        allMeanMaps{iFamily} = analysisResult.groupMeans.raw.(familyDefs(iFamily).name);
    end
    diffCLim = cfGCUtils.resolveDiffCLim(allMeanMaps,cfg.diffColorLimits);

    for iFamily = 1:numel(familyDefs)
        familyName = familyDefs(iFamily).name;
        figureFiles.(familyDefs(iFamily).fileStem) = fullfile(cfg.figureFolder,[differenceToken '_' familyDefs(iFamily).fileStem '.png']);
        hFig = cfGCUtils.createBandPairGridFigure( ...
            analysisResult.groupMeans.raw.(familyName), ...
            analysisResult.pairInfo, ...
            analysisResult.roiLabels, ...
            diffCLim, ...
            sprintf('%s | %s',differenceLabel,familyDefs(iFamily).figureTitle), ...
            'cfGCUtils.blueWhiteRed', ...
            analysisResult.inference.(familyName).sig, ...
            analysisResult.inference.(familyName).q);
        cfGCUtils.saveFigureQuietly(hFig,figureFiles.(familyDefs(iFamily).fileStem));
        outputManifest{end+1,1} = figureFiles.(familyDefs(iFamily).fileStem); %#ok<AGROW>
    end

    qSummaryFileName = fullfile(analysisRoot,[differenceToken '_q_value_summary.csv']);
    writetable(analysisResult.qValueSummaryTable,qSummaryFileName);
    outputManifest{end+1,1} = qSummaryFileName; %#ok<AGROW>
    analysisResult.qValueSummaryFileName = qSummaryFileName;

    if cfg.makeOverviewFigure
        figureFiles.overview = fullfile(cfg.figureFolder,[differenceToken '_overviewOnly.png']);
        cfGCUtils.saveFigureQuietly(createOverviewFigure(analysisResult,differenceLabel),figureFiles.overview);
        outputManifest{end+1,1} = figureFiles.overview; %#ok<AGROW>
    end
end

analysisResult.figureFiles = figureFiles;
analysisResult.outputManifest = outputManifest;
analysisResult.notes = ['Protocol difference uses only matched meditator-control pairs with complete data in both protocols. ', ...
    'Displayed heatmaps are raw mean differences (' differenceLabel '). ', ...
    'Inference and q-values are computed from subject-normalized matched-pair protocol label swaps, ', ...
    'with BH-FDR applied separately within each band pair for every plot family.'];

cfGCUtils.saveAnalysisSummary(summaryFileName,analysisResult);
analysisResult.outputManifest{end+1,1} = summaryFileName;
cfGCUtils.saveAnalysisSummary(summaryFileName,analysisResult);
end

function cfg = applyDefaultCfg(cfgIn)
thisFolder = fileparts(mfilename('fullpath'));

cfg = cfgIn;
cfg = setDefault(cfg,'dataFolder',fullfile(thisFolder,'savedDataCrossFreqGranger'));
cfg = setDefault(cfg,'analysisFolder',fullfile(thisFolder,'analysis'));
cfg = setDefault(cfg,'analysisSubfolder',fullfile('comparison','protocol_difference'));
cfg = setDefault(cfg,'numPermutations',5000);
cfg = setDefault(cfg,'minObservationCount',20);
cfg = setDefault(cfg,'minValidSubjectsPerCell',3);
cfg = setDefault(cfg,'fdrAlpha',0.05);
cfg = setDefault(cfg,'diffColorLimits',[]);
cfg = setDefault(cfg,'makeOverviewFigure',1);
cfg = setDefault(cfg,'saveOutputFlag',1);
cfg = setDefault(cfg,'forceRebuild',0);
end

function protocolName = validateSingleProtocolName(protocolNameIn)
protocolList = cfGCUtils.resolveProtocolNameList(protocolNameIn);
if numel(protocolList) ~= 1
    error('Each protocol input must be a single protocol name.');
end
protocolName = protocolList{1};
end

function protocolData = prepareProtocolArrays(pairList,protocolName,badEyeCondition,badTrialVersion,cfg)
numPairsRequested = numel(pairList);
validPairMask = false(numPairsRequested,1);
pairAvailability = repmat({'missing'},numPairsRequested,1);
templateResults = [];

for iPair = 1:numPairsRequested
    pairInfoThis = pairList(iPair);
    try
        medPreResults = cfGCUtils.loadConditionResults(pairInfoThis.meditator,protocolName,'pre',badEyeCondition,badTrialVersion,cfg.dataFolder);
        medPostResults = cfGCUtils.loadConditionResults(pairInfoThis.meditator,protocolName,'post',badEyeCondition,badTrialVersion,cfg.dataFolder);
        ctrlPreResults = cfGCUtils.loadConditionResults(pairInfoThis.control,protocolName,'pre',badEyeCondition,badTrialVersion,cfg.dataFolder);
        ctrlPostResults = cfGCUtils.loadConditionResults(pairInfoThis.control,protocolName,'post',badEyeCondition,badTrialVersion,cfg.dataFolder);
    catch
        continue;
    end

    cfGCUtils.validateResultCompatibility(medPreResults,medPostResults);
    cfGCUtils.validateResultCompatibility(ctrlPreResults,ctrlPostResults);
    cfGCUtils.validateResultCompatibility(medPreResults,ctrlPreResults);
    cfGCUtils.validateResultCompatibility(medPostResults,ctrlPostResults);

    if isempty(templateResults)
        templateResults = medPreResults;
        [~,~,pairInfo,phaseBandList,ampBandList] = cfGCUtils.unpackBandPairGrid(medPreResults);
        numROIs = numel(medPreResults.roiLabels);
        numPhaseBands = numel(phaseBandList);
        numAmpBands = numel(ampBandList);
        [raw,normalized,sampleCount] = initializeProtocolArrays(numROIs,numPhaseBands,numAmpBands,numPairsRequested);
    end

    [medPreGrid,medPreSample] = cfGCUtils.unpackBandPairGrid(medPreResults);
    [medPostGrid,medPostSample] = cfGCUtils.unpackBandPairGrid(medPostResults);
    [ctrlPreGrid,ctrlPreSample] = cfGCUtils.unpackBandPairGrid(ctrlPreResults);
    [ctrlPostGrid,ctrlPostSample] = cfGCUtils.unpackBandPairGrid(ctrlPostResults);

    lowSampleMask = medPreSample < cfg.minObservationCount | medPostSample < cfg.minObservationCount | ...
        ctrlPreSample < cfg.minObservationCount | ctrlPostSample < cfg.minObservationCount;

    medPreGrid(lowSampleMask) = NaN;
    medPostGrid(lowSampleMask) = NaN;
    ctrlPreGrid(lowSampleMask) = NaN;
    ctrlPostGrid(lowSampleMask) = NaN;

    medPreNorm = cfGCUtils.normalizeSubjectGrid(medPreGrid);
    medPostNorm = cfGCUtils.normalizeSubjectGrid(medPostGrid);
    ctrlPreNorm = cfGCUtils.normalizeSubjectGrid(ctrlPreGrid);
    ctrlPostNorm = cfGCUtils.normalizeSubjectGrid(ctrlPostGrid);

    raw.preMed(:,:,:,:,iPair) = medPreGrid;
    raw.postMed(:,:,:,:,iPair) = medPostGrid;
    raw.preCtrl(:,:,:,:,iPair) = ctrlPreGrid;
    raw.postCtrl(:,:,:,:,iPair) = ctrlPostGrid;
    raw.deltaMed(:,:,:,:,iPair) = medPostGrid - medPreGrid;
    raw.deltaCtrl(:,:,:,:,iPair) = ctrlPostGrid - ctrlPreGrid;

    normalized.preMed(:,:,:,:,iPair) = medPreNorm;
    normalized.postMed(:,:,:,:,iPair) = medPostNorm;
    normalized.preCtrl(:,:,:,:,iPair) = ctrlPreNorm;
    normalized.postCtrl(:,:,:,:,iPair) = ctrlPostNorm;
    normalized.deltaMed(:,:,:,:,iPair) = medPostNorm - medPreNorm;
    normalized.deltaCtrl(:,:,:,:,iPair) = ctrlPostNorm - ctrlPreNorm;

    sampleCount.minAcrossAll(:,:,:,:,iPair) = min(cat(5,medPreSample,medPostSample,ctrlPreSample,ctrlPostSample),[],5);
    validPairMask(iPair) = true;
    pairAvailability{iPair} = 'ok';
end

if isempty(templateResults)
    error('No valid pairs with complete PRE/POST CF-GC data were found for protocol %s.',protocolName);
end

protocolData = struct();
protocolData.protocolName = protocolName;
protocolData.templateResults = templateResults;
protocolData.roiLabels = templateResults.roiLabels;
protocolData.pairInfo = pairInfo;
protocolData.phaseBandList = phaseBandList;
protocolData.ampBandList = ampBandList;
protocolData.validPairMask = validPairMask;
protocolData.pairAvailability = pairAvailability;
protocolData.raw = raw;
protocolData.normalized = normalized;
protocolData.sampleCount = sampleCount;
end

function [raw,normalized,sampleCount] = initializeProtocolArrays(numROIs,numPhaseBands,numAmpBands,numPairsRequested)
raw = struct();
raw.preMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
raw.postMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
raw.preCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
raw.postCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
raw.deltaMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
raw.deltaCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);

normalized = struct();
normalized.preMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
normalized.postMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
normalized.preCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
normalized.postCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
normalized.deltaMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
normalized.deltaCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);

sampleCount = struct();
sampleCount.minAcrossAll = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
end

function metric = extractFamilyMetric(dataStruct,familyName)
switch familyName
    case 'meditatorPre'
        metric = dataStruct.preMed;
    case 'meditatorPost'
        metric = dataStruct.postMed;
    case 'controlPre'
        metric = dataStruct.preCtrl;
    case 'controlPost'
        metric = dataStruct.postCtrl;
    case 'meditatorDelta'
        metric = dataStruct.deltaMed;
    case 'controlDelta'
        metric = dataStruct.deltaCtrl;
    case 'baselinePreDifference'
        metric = dataStruct.preMed - dataStruct.preCtrl;
    case 'postDifference'
        metric = dataStruct.postMed - dataStruct.postCtrl;
    case 'interaction'
        metric = dataStruct.deltaMed - dataStruct.deltaCtrl;
    otherwise
        error('Unknown family name: %s',familyName);
end
end

function inference = computeProtocolDifferenceInference(protocolAData,protocolBData,commonIndices,familyDefs,cfg)
dataSize = size(protocolAData.normalized.preMed);
numROIs = dataSize(1);
numPhase = dataSize(3);
numAmp = dataSize(4);

inference = struct();
for iFamily = 1:numel(familyDefs)
    familyName = familyDefs(iFamily).name;
    inference.(familyName) = struct( ...
        'p',nan(numROIs,numROIs,numPhase,numAmp), ...
        'effectDz',nan(numROIs,numROIs,numPhase,numAmp), ...
        'n',nan(numROIs,numROIs,numPhase,numAmp), ...
        'q',nan(numROIs,numROIs,numPhase,numAmp), ...
        'sig',false(numROIs,numROIs,numPhase,numAmp));
end

for iFamily = 1:numel(familyDefs)
    familyName = familyDefs(iFamily).name;
    metricA = extractFamilyMetric(protocolAData.normalized,familyName);
    metricB = extractFamilyMetric(protocolBData.normalized,familyName);

    for iROI = 1:numROIs
        for jROI = 1:numROIs
            for iPhase = 1:numPhase
                for iAmp = 1:numAmp
                    valuesA = squeeze(metricA(iROI,jROI,iPhase,iAmp,commonIndices));
                    valuesB = squeeze(metricB(iROI,jROI,iPhase,iAmp,commonIndices));
                    [pVal,effectDz,nUsed] = cfGCUtils.pairedLabelSwapPrePost(valuesB,valuesA,cfg.numPermutations);
                    inference.(familyName).p(iROI,jROI,iPhase,iAmp) = pVal;
                    inference.(familyName).effectDz(iROI,jROI,iPhase,iAmp) = effectDz;
                    inference.(familyName).n(iROI,jROI,iPhase,iAmp) = nUsed;
                end
            end
        end
    end

    invalidMask = inference.(familyName).n < cfg.minValidSubjectsPerCell;
    inference.(familyName).p(invalidMask) = NaN;
    inference.(familyName).effectDz(invalidMask) = NaN;
    [inference.(familyName).q,inference.(familyName).sig] = cfGCUtils.bhFDRByPair(inference.(familyName).p,cfg.fdrAlpha);
end
end

function summaryOverview = buildSummaryOverview(groupMeans,inference,commonSampleCount,familyDefs)
summaryOverview = struct();
summaryOverview.collapsed = struct();
summaryOverview.minQ = struct();

for iFamily = 1:numel(familyDefs)
    familyName = familyDefs(iFamily).name;
    summaryOverview.collapsed.(familyName) = cfGCUtils.meanAcrossPairs(groupMeans.(familyName));
    summaryOverview.minQ.(familyName) = cfGCUtils.minAcrossPairs(inference.(familyName).q);
end

summaryOverview.minSampleCount = cfGCUtils.meanAcrossPairs(mean(commonSampleCount,5,'omitnan'));
summaryOverview.sigCountInteraction = cfGCUtils.sumAcrossPairs(double(inference.interaction.sig));
summaryOverview.netFlow = struct();
summaryOverview.netFlow.meditatorPre = cfGCUtils.computeNetFlow(summaryOverview.collapsed.meditatorPre);
summaryOverview.netFlow.meditatorPost = cfGCUtils.computeNetFlow(summaryOverview.collapsed.meditatorPost);
summaryOverview.netFlow.controlPre = cfGCUtils.computeNetFlow(summaryOverview.collapsed.controlPre);
summaryOverview.netFlow.controlPost = cfGCUtils.computeNetFlow(summaryOverview.collapsed.controlPost);
summaryOverview.netFlow.meditatorDelta = cfGCUtils.computeNetFlow(summaryOverview.collapsed.meditatorDelta);
summaryOverview.netFlow.controlDelta = cfGCUtils.computeNetFlow(summaryOverview.collapsed.controlDelta);
summaryOverview.netFlow.interaction = cfGCUtils.computeNetFlow(summaryOverview.collapsed.interaction);
end

function qTable = makeQSummaryTable(familyDefs,pairInfo,inference)
numFamilies = numel(familyDefs);
numPhase = size(pairInfo,1);
numAmp = size(pairInfo,2);
numRows = numFamilies * numPhase * numAmp;

familyColumn = cell(numRows,1);
phaseColumn = nan(numRows,1);
ampColumn = nan(numRows,1);
bandLabelColumn = cell(numRows,1);
minQColumn = nan(numRows,1);
sigCountColumn = nan(numRows,1);

rowIndex = 0;
for iFamily = 1:numFamilies
    familyName = familyDefs(iFamily).name;
    for iPhase = 1:numPhase
        for iAmp = 1:numAmp
            rowIndex = rowIndex + 1;
            qThis = inference.(familyName).q(:,:,iPhase,iAmp);
            sigThis = inference.(familyName).sig(:,:,iPhase,iAmp);
            familyColumn{rowIndex,1} = familyName;
            phaseColumn(rowIndex,1) = iPhase;
            ampColumn(rowIndex,1) = iAmp;
            bandLabelColumn{rowIndex,1} = pairInfo(iPhase,iAmp).label;
            minQColumn(rowIndex,1) = cfGCUtils.minFiniteValue(qThis);
            sigCountColumn(rowIndex,1) = sum(sigThis(:));
        end
    end
end

qTable = table(familyColumn,phaseColumn,ampColumn,bandLabelColumn,minQColumn,sigCountColumn, ...
    'VariableNames',{'Family','PhaseBandIndex','AmplitudeBandIndex','BandLabel','MinQ','NumSignificantCells'});
end

function hFig = createOverviewFigure(analysisResult,differenceLabel)
hFig = cfGCUtils.createFigure([differenceLabel ' overview'],[100 100 1700 980]);
tiled = tiledlayout(hFig,3,4,'TileSpacing','compact','Padding','compact');
overview = analysisResult.summaryOverview;

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsed.meditatorPre,analysisResult.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (meditator PRE difference)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsed.meditatorPost,analysisResult.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (meditator POST difference)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsed.meditatorDelta,analysisResult.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (meditator delta difference)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.minQ.meditatorDelta,analysisResult.roiLabels,[0 0.1],'parula');
title(ax,'Overview only: min meditator-delta q');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsed.controlPre,analysisResult.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (control PRE difference)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsed.controlPost,analysisResult.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (control POST difference)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsed.controlDelta,analysisResult.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (control delta difference)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.minQ.controlDelta,analysisResult.roiLabels,[0 0.1],'parula');
title(ax,'Overview only: min control-delta q');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsed.baselinePreDifference,analysisResult.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (PRE group-difference change)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsed.postDifference,analysisResult.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (POST group-difference change)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsed.interaction,analysisResult.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (delta group-difference change)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.minQ.interaction,analysisResult.roiLabels,[0 0.1],'parula');
title(ax,'Overview only: min interaction q');

cfGCUtils.addFigureTitle(hFig,sprintf('%s | Protocol-difference supplementary overview',differenceLabel));
end

function familyDefs = getFamilyDefinitions()
familyDefs = struct( ...
    'name', { ...
        'meditatorPre', ...
        'meditatorPost', ...
        'controlPre', ...
        'controlPost', ...
        'meditatorDelta', ...
        'controlDelta', ...
        'baselinePreDifference', ...
        'postDifference', ...
        'interaction'}, ...
    'fileStem', { ...
        'meditatorPre', ...
        'meditatorPost', ...
        'controlPre', ...
        'controlPost', ...
        'meditatorDelta', ...
        'controlDelta', ...
        'preDifference', ...
        'postDifference', ...
        'deltaDifference'}, ...
    'figureTitle', { ...
        'Selected meditator PRE mean difference', ...
        'Selected meditator POST mean difference', ...
        'Selected control PRE mean difference', ...
        'Selected control POST mean difference', ...
        'Selected meditator POST - PRE mean difference', ...
        'Selected control POST - PRE mean difference', ...
        'PRE Meditator - Control difference', ...
        'POST Meditator - Control difference', ...
        'Delta Meditator - Control difference'});
end

function analysisRoot = buildAnalysisRoot(baseFolder,subFolder,analysisTag)
if isempty(subFolder)
    analysisRoot = fullfile(baseFolder,analysisTag);
else
    analysisRoot = fullfile(baseFolder,subFolder,analysisTag);
end
end

function s = setDefault(s,fieldName,defaultValue)
if ~isfield(s,fieldName) || isempty(s.(fieldName))
    s.(fieldName) = defaultValue;
end
end
