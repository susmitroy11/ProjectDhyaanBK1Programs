function netFlowResult = runNetFlowAnalysisCrossFrequencyPhaseAmplitudeGrangerBK1(pairSelection,protocolNameList,badEyeCondition,badTrialVersion,cfgIn)
% runNetFlowAnalysisCrossFrequencyPhaseAmplitudeGrangerBK1
%
% Run ROI-wise and band-wise net-flow analysis on saved CF-GC outputs.
% This script recomputes directed net flow from the saved PRE/POST CF-GC
% grids, performs matched-pair inference, and writes per-protocol plus
% cross-protocol figures to analysis_net_flow.
%
% Example:
%   runNetFlowAnalysisCrossFrequencyPhaseAmplitudeGrangerBK1('all','all','ep','v8')
%
% Outputs:
%   <cf_gc>/analysis_net_flow/
%       net_flow_summary.mat
%       figures/
%           overview_protocol_by_roi_interaction.png
%           overview_protocol_by_band_interaction.png
%           overview_frontal_dominance.png
%           <protocol>/
%               collapsed_connectivity_matrices.png
%               roi_net_flow_summary.png
%               band_resolved_interaction.png
%               frontal_dominance.png

if ~exist('pairSelection','var') || isempty(pairSelection);       pairSelection = 'all'; end
if ~exist('protocolNameList','var') || isempty(protocolNameList); protocolNameList = 'all'; end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition);   badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion);   badTrialVersion = 'v8'; end
if ~exist('cfgIn','var') || isempty(cfgIn);                       cfgIn = struct; end

cfg = applyDefaultCfg(cfgIn);
rng(cfg.randomSeed,'twister');

thisFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(thisFolder);
protocolNameList = cfGCUtils.resolveProtocolNameList(protocolNameList);
[pairList,selectionLabel] = cfGCUtils.resolvePairSelectionBK1(pairSelection,projectRoot);

cfGCUtils.ensureFolder(cfg.outputFolder);
cfGCUtils.ensureFolder(cfg.figureFolder);
cfGCUtils.ensureFolder(cfg.checkpointFolder);

preparedData = prepareNetFlowInputData(pairList,protocolNameList,badEyeCondition,badTrialVersion,cfg);
netFlowSummary = analyzeAllProtocols(preparedData,cfg);
overviewFigureFiles = createCrossProtocolFigures(netFlowSummary,cfg);
saveRunCheckpoint(netFlowSummary,overviewFigureFiles,cfg,'completed');

summaryFile = fullfile(cfg.outputFolder,'net_flow_summary.mat');
save(summaryFile,'netFlowSummary','-v7.3');

netFlowResult = struct;
netFlowResult.outputFolder = cfg.outputFolder;
netFlowResult.figureFolder = cfg.figureFolder;
netFlowResult.summaryFile = summaryFile;
netFlowResult.overviewFigureFiles = overviewFigureFiles;
netFlowResult.selectionLabel = selectionLabel;
netFlowResult.protocolNameList = protocolNameList;
netFlowResult.badEyeCondition = badEyeCondition;
netFlowResult.badTrialVersion = badTrialVersion;
netFlowResult.numRequestedPairs = numel(pairList);
netFlowResult.config = cfg;

fprintf('Net-flow analysis completed. Results saved under:\n%s\n',cfg.outputFolder);
end

function cfg = applyDefaultCfg(cfgIn)
thisFolder = fileparts(mfilename('fullpath'));

cfg = cfgIn;
cfg = setDefault(cfg,'dataFolder',fullfile(thisFolder,'savedDataCrossFreqGranger'));
cfg = setDefault(cfg,'outputFolder',fullfile(thisFolder,'analysis_net_flow'));
cfg = setDefault(cfg,'figureFolder',fullfile(thisFolder,'analysis_net_flow','figures'));
cfg = setDefault(cfg,'checkpointFolder',fullfile(thisFolder,'analysis_net_flow','checkpoints'));
cfg = setDefault(cfg,'numPermutations',5000);
cfg = setDefault(cfg,'minObservationCount',20);
cfg = setDefault(cfg,'minValidPairs',3);
cfg = setDefault(cfg,'fdrAlpha',0.05);
cfg = setDefault(cfg,'randomSeed',1);
cfg = setDefault(cfg,'maxROIsPerBandFigure',18);
cfg = setDefault(cfg,'figureFormat','png');
cfg = setDefault(cfg,'frontalLabelPattern','Frontal');
cfg = setDefault(cfg,'resumeIfAvailable',1);
end

function preparedData = prepareNetFlowInputData(pairList,protocolNameList,badEyeCondition,badTrialVersion,cfg)
preparedData = struct;
preparedData.createdOn = datetime('now');
preparedData.protocolNameList = protocolNameList;
preparedData.badEyeCondition = badEyeCondition;
preparedData.badTrialVersion = badTrialVersion;
preparedData.requestedPairList = pairList;
preparedData.protocolResults = [];

for iProtocol = 1:numel(protocolNameList)
    protocolData = prepareSingleProtocol(pairList,protocolNameList{iProtocol},badEyeCondition,badTrialVersion,cfg);
    if isempty(preparedData.protocolResults)
        preparedData.protocolResults = repmat(protocolData,1,numel(protocolNameList));
    end
    preparedData.protocolResults(iProtocol) = protocolData;
end

preparedData.roiLabels = preparedData.protocolResults(1).roiLabels;
preparedData.pairInfo = preparedData.protocolResults(1).pairInfo;
preparedData.phaseBandList = preparedData.protocolResults(1).phaseBandList;
preparedData.ampBandList = preparedData.protocolResults(1).ampBandList;
end

function protocolData = prepareSingleProtocol(pairList,protocolName,badEyeCondition,badTrialVersion,cfg)
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

        raw.preMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        raw.postMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        raw.preCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        raw.postCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        raw.deltaMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        raw.deltaCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        sampleCount.minAcrossAll = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
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

    raw.preMed(:,:,:,:,iPair) = medPreGrid;
    raw.postMed(:,:,:,:,iPair) = medPostGrid;
    raw.preCtrl(:,:,:,:,iPair) = ctrlPreGrid;
    raw.postCtrl(:,:,:,:,iPair) = ctrlPostGrid;
    raw.deltaMed(:,:,:,:,iPair) = medPostGrid - medPreGrid;
    raw.deltaCtrl(:,:,:,:,iPair) = ctrlPostGrid - ctrlPreGrid;
    sampleCount.minAcrossAll(:,:,:,:,iPair) = min(cat(5,medPreSample,medPostSample,ctrlPreSample,ctrlPostSample),[],5);

    validPairMask(iPair) = true;
    pairAvailability{iPair} = 'ok';
end

if isempty(templateResults)
    error('No valid pairs with complete PRE/POST CF-GC data were found for protocol %s.',protocolName);
end

validIndices = find(validPairMask);
protocolData = struct;
protocolData.protocolName = protocolName;
protocolData.numRequestedPairs = numPairsRequested;
protocolData.numValidPairs = numel(validIndices);
protocolData.validIndices = validIndices;
protocolData.validPairMask = validPairMask;
protocolData.requestedPairList = pairList;
protocolData.validPairList = pairList(validIndices);
protocolData.validPairLabels = {protocolData.validPairList.pairLabel};
protocolData.pairAvailability = pairAvailability;
protocolData.roiLabels = templateResults.roiLabels;
protocolData.pairInfo = pairInfo;
protocolData.phaseBandList = phaseBandList;
protocolData.ampBandList = ampBandList;
protocolData.raw = restrictToValidPairs(raw,validIndices);
protocolData.sampleCount = restrictToValidPairs(sampleCount,validIndices);
end

function out = restrictToValidPairs(in,validIndices)
fieldNames = fieldnames(in);
out = struct;
for iField = 1:numel(fieldNames)
    out.(fieldNames{iField}) = in.(fieldNames{iField})(:,:,:,:,validIndices);
end
end

function netFlowSummary = analyzeAllProtocols(preparedData,cfg)
numProtocols = numel(preparedData.protocolResults);
netFlowSummary = struct;
netFlowSummary.analysisType = 'cross_frequency_phase_amplitude_granger_net_flow';
netFlowSummary.createdOn = datetime('now');
netFlowSummary.config = cfg;
netFlowSummary.protocolNameList = preparedData.protocolNameList;
netFlowSummary.roiLabels = preparedData.roiLabels;
netFlowSummary.phaseBandList = preparedData.phaseBandList;
netFlowSummary.ampBandList = preparedData.ampBandList;
netFlowSummary.bandPairLabels = makeBandPairLabels(preparedData.pairInfo);
netFlowSummary.protocolResults = [];
netFlowSummary.completedProtocols = false(1,numProtocols);

for iProtocol = 1:numProtocols
    protocolData = preparedData.protocolResults(iProtocol);
    protocolResult = [];
    protocolCheckpointFile = getProtocolCheckpointFile(protocolData.protocolName,cfg);

    if cfg.resumeIfAvailable && exist(protocolCheckpointFile,'file')
        tmp = load(protocolCheckpointFile,'protocolResult');
        if isfield(tmp,'protocolResult')
            protocolResult = tmp.protocolResult;
        end
    end

    if isempty(protocolResult)
        protocolResult = analyzeSingleProtocol(protocolData,cfg);
        save(protocolCheckpointFile,'protocolResult','-v7.3');
    end

    if isempty(netFlowSummary.protocolResults)
        netFlowSummary.protocolResults = repmat(protocolResult,1,numProtocols);
    end
    netFlowSummary.protocolResults(iProtocol) = protocolResult;
    netFlowSummary.completedProtocols(iProtocol) = true;
    saveRunCheckpoint(netFlowSummary,struct,cfg,sprintf('completed_%s',protocolData.protocolName));
end
end

function protocolResult = analyzeSingleProtocol(protocolData,cfg)
roiLabels = protocolData.roiLabels;
numROIs = numel(roiLabels);
numPairs = protocolData.numValidPairs;
numPhaseBands = numel(protocolData.phaseBandList);
numAmpBands = numel(protocolData.ampBandList);

collapsedMatrices = struct;
collapsedMatrices.meditatorPre = nan(numROIs,numROIs,numPairs);
collapsedMatrices.meditatorPost = nan(numROIs,numROIs,numPairs);
collapsedMatrices.controlPre = nan(numROIs,numROIs,numPairs);
collapsedMatrices.controlPost = nan(numROIs,numROIs,numPairs);
collapsedMatrices.meditatorDelta = nan(numROIs,numROIs,numPairs);
collapsedMatrices.controlDelta = nan(numROIs,numROIs,numPairs);

roiNetFlow = struct;
roiNetFlow.meditatorPre = nan(numROIs,numPairs);
roiNetFlow.meditatorPost = nan(numROIs,numPairs);
roiNetFlow.controlPre = nan(numROIs,numPairs);
roiNetFlow.controlPost = nan(numROIs,numPairs);
roiNetFlow.meditatorDelta = nan(numROIs,numPairs);
roiNetFlow.controlDelta = nan(numROIs,numPairs);

bandNetFlow = struct;
bandNetFlow.meditatorPre = nan(numROIs,numPhaseBands,numAmpBands,numPairs);
bandNetFlow.meditatorPost = nan(numROIs,numPhaseBands,numAmpBands,numPairs);
bandNetFlow.controlPre = nan(numROIs,numPhaseBands,numAmpBands,numPairs);
bandNetFlow.controlPost = nan(numROIs,numPhaseBands,numAmpBands,numPairs);
bandNetFlow.meditatorDelta = nan(numROIs,numPhaseBands,numAmpBands,numPairs);
bandNetFlow.controlDelta = nan(numROIs,numPhaseBands,numAmpBands,numPairs);

for iPair = 1:numPairs
    collapsedMatrices.meditatorPre(:,:,iPair) = collapseBandGrid(protocolData.raw.preMed(:,:,:,:,iPair));
    collapsedMatrices.meditatorPost(:,:,iPair) = collapseBandGrid(protocolData.raw.postMed(:,:,:,:,iPair));
    collapsedMatrices.controlPre(:,:,iPair) = collapseBandGrid(protocolData.raw.preCtrl(:,:,:,:,iPair));
    collapsedMatrices.controlPost(:,:,iPair) = collapseBandGrid(protocolData.raw.postCtrl(:,:,:,:,iPair));
    collapsedMatrices.meditatorDelta(:,:,iPair) = collapseBandGrid(protocolData.raw.deltaMed(:,:,:,:,iPair));
    collapsedMatrices.controlDelta(:,:,iPair) = collapseBandGrid(protocolData.raw.deltaCtrl(:,:,:,:,iPair));

    roiNetFlow.meditatorPre(:,iPair) = computeNetFlowVector(collapsedMatrices.meditatorPre(:,:,iPair));
    roiNetFlow.meditatorPost(:,iPair) = computeNetFlowVector(collapsedMatrices.meditatorPost(:,:,iPair));
    roiNetFlow.controlPre(:,iPair) = computeNetFlowVector(collapsedMatrices.controlPre(:,:,iPair));
    roiNetFlow.controlPost(:,iPair) = computeNetFlowVector(collapsedMatrices.controlPost(:,:,iPair));
    roiNetFlow.meditatorDelta(:,iPair) = computeNetFlowVector(collapsedMatrices.meditatorDelta(:,:,iPair));
    roiNetFlow.controlDelta(:,iPair) = computeNetFlowVector(collapsedMatrices.controlDelta(:,:,iPair));

    bandNetFlow.meditatorPre(:,:,:,iPair) = computeBandResolvedNetFlow(protocolData.raw.preMed(:,:,:,:,iPair));
    bandNetFlow.meditatorPost(:,:,:,iPair) = computeBandResolvedNetFlow(protocolData.raw.postMed(:,:,:,:,iPair));
    bandNetFlow.controlPre(:,:,:,iPair) = computeBandResolvedNetFlow(protocolData.raw.preCtrl(:,:,:,:,iPair));
    bandNetFlow.controlPost(:,:,:,iPair) = computeBandResolvedNetFlow(protocolData.raw.postCtrl(:,:,:,:,iPair));
    bandNetFlow.meditatorDelta(:,:,:,iPair) = computeBandResolvedNetFlow(protocolData.raw.deltaMed(:,:,:,:,iPair));
    bandNetFlow.controlDelta(:,:,:,iPair) = computeBandResolvedNetFlow(protocolData.raw.deltaCtrl(:,:,:,:,iPair));
end

roiStats = computeROIStatistics(roiNetFlow,cfg);
bandStats = computeBandStatistics(bandNetFlow,cfg);
frontalStats = computeFrontalDominance(roiNetFlow,roiLabels,cfg);
meanMatrices = computeMeanMatrices(collapsedMatrices);
meanBandInteraction = mean(bandNetFlow.meditatorDelta - bandNetFlow.controlDelta,4,'omitnan');

protocolFigureFolder = fullfile(cfg.figureFolder,protocolData.protocolName);
cfGCUtils.ensureFolder(protocolFigureFolder);

figureFiles = struct;
figureFiles.collapsedMatrices = plotCollapsedMatrices(protocolData.protocolName,roiLabels,meanMatrices,protocolFigureFolder,cfg);
figureFiles.roiSummary = plotROINetFlowSummary(protocolData.protocolName,roiLabels,roiNetFlow,roiStats,protocolFigureFolder,cfg);
figureFiles.bandInteraction = plotBandResolvedInteraction(protocolData.protocolName,roiLabels,protocolData.pairInfo,meanBandInteraction,bandStats,protocolFigureFolder,cfg);
figureFiles.frontalDominance = plotFrontalDominance(protocolData.protocolName,frontalStats,protocolFigureFolder,cfg);

protocolResult = struct;
protocolResult.protocolName = protocolData.protocolName;
protocolResult.validPairLabels = protocolData.validPairLabels;
protocolResult.numValidPairs = protocolData.numValidPairs;
protocolResult.roiLabels = roiLabels;
protocolResult.phaseBandList = protocolData.phaseBandList;
protocolResult.ampBandList = protocolData.ampBandList;
protocolResult.bandPairLabels = makeBandPairLabels(protocolData.pairInfo);
protocolResult.collapsedMatrices = collapsedMatrices;
protocolResult.meanMatrices = meanMatrices;
protocolResult.roiNetFlow = roiNetFlow;
protocolResult.bandNetFlow = bandNetFlow;
protocolResult.roiStats = roiStats;
protocolResult.bandStats = bandStats;
protocolResult.frontalStats = frontalStats;
protocolResult.meanBandInteraction = meanBandInteraction;
protocolResult.figureFiles = figureFiles;
end

function meanMatrices = computeMeanMatrices(collapsedMatrices)
meanMatrices = struct;
fieldNames = fieldnames(collapsedMatrices);
for iField = 1:numel(fieldNames)
    meanMatrices.(fieldNames{iField}) = mean(collapsedMatrices.(fieldNames{iField}),3,'omitnan');
end
meanMatrices.preGroupDifference = meanMatrices.meditatorPre - meanMatrices.controlPre;
meanMatrices.postGroupDifference = meanMatrices.meditatorPost - meanMatrices.controlPost;
meanMatrices.interaction = meanMatrices.meditatorDelta - meanMatrices.controlDelta;
end

function roiStats = computeROIStatistics(roiNetFlow,cfg)
numROIs = size(roiNetFlow.meditatorPre,1);
families = { ...
    'meditatorChange','meditatorPre','meditatorPost'; ...
    'controlChange','controlPre','controlPost'; ...
    'baselinePreGroup','meditatorPre','controlPre'; ...
    'postGroup','meditatorPost','controlPost'; ...
    'interaction','meditatorDelta','controlDelta'};

roiStats = struct;
for iFamily = 1:size(families,1)
    familyName = families{iFamily,1};
    roiStats.(familyName).p = nan(numROIs,1);
    roiStats.(familyName).effect = nan(numROIs,1);
    roiStats.(familyName).n = nan(numROIs,1);
end

for iROI = 1:numROIs
    [roiStats.meditatorChange.p(iROI),roiStats.meditatorChange.effect(iROI),roiStats.meditatorChange.n(iROI)] = ...
        cfGCUtils.pairedLabelSwapPrePost(roiNetFlow.meditatorPre(iROI,:),roiNetFlow.meditatorPost(iROI,:),cfg.numPermutations);
    [roiStats.controlChange.p(iROI),roiStats.controlChange.effect(iROI),roiStats.controlChange.n(iROI)] = ...
        cfGCUtils.pairedLabelSwapPrePost(roiNetFlow.controlPre(iROI,:),roiNetFlow.controlPost(iROI,:),cfg.numPermutations);
    [roiStats.baselinePreGroup.p(iROI),roiStats.baselinePreGroup.effect(iROI),roiStats.baselinePreGroup.n(iROI)] = ...
        cfGCUtils.pairedLabelSwapMatchedGroups(roiNetFlow.meditatorPre(iROI,:),roiNetFlow.controlPre(iROI,:),cfg.numPermutations);
    [roiStats.postGroup.p(iROI),roiStats.postGroup.effect(iROI),roiStats.postGroup.n(iROI)] = ...
        cfGCUtils.pairedLabelSwapMatchedGroups(roiNetFlow.meditatorPost(iROI,:),roiNetFlow.controlPost(iROI,:),cfg.numPermutations);
    [roiStats.interaction.p(iROI),roiStats.interaction.effect(iROI),roiStats.interaction.n(iROI)] = ...
        cfGCUtils.pairedLabelSwapMatchedGroups(roiNetFlow.meditatorDelta(iROI,:),roiNetFlow.controlDelta(iROI,:),cfg.numPermutations);
end

familyNames = fieldnames(roiStats);
for iFamily = 1:numel(familyNames)
    invalidMask = roiStats.(familyNames{iFamily}).n < cfg.minValidPairs;
    roiStats.(familyNames{iFamily}).p(invalidMask) = NaN;
    roiStats.(familyNames{iFamily}).effect(invalidMask) = NaN;
    [roiStats.(familyNames{iFamily}).q,roiStats.(familyNames{iFamily}).sig] = applyBHFDR(roiStats.(familyNames{iFamily}).p,cfg.fdrAlpha);
end

roiStats.summary = struct;
roiStats.summary.meanMeditatorPre = mean(roiNetFlow.meditatorPre,2,'omitnan');
roiStats.summary.meanMeditatorPost = mean(roiNetFlow.meditatorPost,2,'omitnan');
roiStats.summary.meanControlPre = mean(roiNetFlow.controlPre,2,'omitnan');
roiStats.summary.meanControlPost = mean(roiNetFlow.controlPost,2,'omitnan');
roiStats.summary.meanMeditatorDelta = mean(roiNetFlow.meditatorDelta,2,'omitnan');
roiStats.summary.meanControlDelta = mean(roiNetFlow.controlDelta,2,'omitnan');
roiStats.summary.meanInteraction = roiStats.summary.meanMeditatorDelta - roiStats.summary.meanControlDelta;
roiStats.summary.semMeditatorPre = computeSEM(roiNetFlow.meditatorPre,2);
roiStats.summary.semMeditatorPost = computeSEM(roiNetFlow.meditatorPost,2);
roiStats.summary.semControlPre = computeSEM(roiNetFlow.controlPre,2);
roiStats.summary.semControlPost = computeSEM(roiNetFlow.controlPost,2);
roiStats.summary.semMeditatorDelta = computeSEM(roiNetFlow.meditatorDelta,2);
roiStats.summary.semControlDelta = computeSEM(roiNetFlow.controlDelta,2);
roiStats.summary.semInteraction = computeSEM(roiNetFlow.meditatorDelta - roiNetFlow.controlDelta,2);
end

function bandStats = computeBandStatistics(bandNetFlow,cfg)
dataDeltaDiff = bandNetFlow.meditatorDelta - bandNetFlow.controlDelta;
[numROIs,numPhaseBands,numAmpBands,~] = size(dataDeltaDiff);

bandStats = struct;
bandStats.interactionMean = mean(dataDeltaDiff,4,'omitnan');
bandStats.interactionSEM = computeSEM(dataDeltaDiff,4);
bandStats.interactionP = nan(numROIs,numPhaseBands,numAmpBands);
bandStats.interactionEffect = nan(numROIs,numPhaseBands,numAmpBands);
bandStats.interactionN = nan(numROIs,numPhaseBands,numAmpBands);

for iROI = 1:numROIs
    for iPhase = 1:numPhaseBands
        for iAmp = 1:numAmpBands
            medDelta = squeeze(bandNetFlow.meditatorDelta(iROI,iPhase,iAmp,:));
            ctrlDelta = squeeze(bandNetFlow.controlDelta(iROI,iPhase,iAmp,:));
            [bandStats.interactionP(iROI,iPhase,iAmp), ...
                bandStats.interactionEffect(iROI,iPhase,iAmp), ...
                bandStats.interactionN(iROI,iPhase,iAmp)] = ...
                cfGCUtils.pairedLabelSwapMatchedGroups(medDelta(:)',ctrlDelta(:)',cfg.numPermutations);
        end
    end
end

invalidMask = bandStats.interactionN < cfg.minValidPairs;
bandStats.interactionP(invalidMask) = NaN;
bandStats.interactionEffect(invalidMask) = NaN;

bandStats.interactionQ = nan(size(bandStats.interactionP));
bandStats.interactionSig = false(size(bandStats.interactionP));
for iROI = 1:numROIs
    [qVals,sigMask] = applyBHFDR(reshape(bandStats.interactionP(iROI,:,:),[],1),cfg.fdrAlpha);
    bandStats.interactionQ(iROI,:,:) = reshape(qVals,1,numPhaseBands,numAmpBands);
    bandStats.interactionSig(iROI,:,:) = reshape(sigMask,1,numPhaseBands,numAmpBands);
end

overallBandInteraction = squeeze(mean(bandStats.interactionMean,1,'omitnan'));
bandStats.overallBandInteraction = overallBandInteraction;
end

function frontalStats = computeFrontalDominance(roiNetFlow,roiLabels,cfg)
frontalIndices = find(contains(roiLabels,cfg.frontalLabelPattern,'IgnoreCase',true));
nonFrontalIndices = setdiff(1:numel(roiLabels),frontalIndices);

frontalStats = struct;
frontalStats.frontalIndices = frontalIndices;
frontalStats.nonFrontalIndices = nonFrontalIndices;
frontalStats.available = ~isempty(frontalIndices) && ~isempty(nonFrontalIndices);

if ~frontalStats.available
    frontalStats.values = struct;
    frontalStats.tests = struct;
    return;
end

frontalStats.values.meditatorPre = mean(roiNetFlow.meditatorPre(frontalIndices,:),1,'omitnan') - mean(roiNetFlow.meditatorPre(nonFrontalIndices,:),1,'omitnan');
frontalStats.values.meditatorPost = mean(roiNetFlow.meditatorPost(frontalIndices,:),1,'omitnan') - mean(roiNetFlow.meditatorPost(nonFrontalIndices,:),1,'omitnan');
frontalStats.values.controlPre = mean(roiNetFlow.controlPre(frontalIndices,:),1,'omitnan') - mean(roiNetFlow.controlPre(nonFrontalIndices,:),1,'omitnan');
frontalStats.values.controlPost = mean(roiNetFlow.controlPost(frontalIndices,:),1,'omitnan') - mean(roiNetFlow.controlPost(nonFrontalIndices,:),1,'omitnan');
frontalStats.values.meditatorDelta = mean(roiNetFlow.meditatorDelta(frontalIndices,:),1,'omitnan') - mean(roiNetFlow.meditatorDelta(nonFrontalIndices,:),1,'omitnan');
frontalStats.values.controlDelta = mean(roiNetFlow.controlDelta(frontalIndices,:),1,'omitnan') - mean(roiNetFlow.controlDelta(nonFrontalIndices,:),1,'omitnan');

frontalStats.tests = struct;
 [frontalStats.tests.meditatorChangeP,frontalStats.tests.meditatorChangeEffect,frontalStats.tests.meditatorChangeN] = ...
    cfGCUtils.pairedLabelSwapPrePost(frontalStats.values.meditatorPre,frontalStats.values.meditatorPost,cfg.numPermutations);
 [frontalStats.tests.controlChangeP,frontalStats.tests.controlChangeEffect,frontalStats.tests.controlChangeN] = ...
    cfGCUtils.pairedLabelSwapPrePost(frontalStats.values.controlPre,frontalStats.values.controlPost,cfg.numPermutations);
 [frontalStats.tests.baselinePreGroupP,frontalStats.tests.baselinePreGroupEffect,frontalStats.tests.baselinePreGroupN] = ...
    cfGCUtils.pairedLabelSwapMatchedGroups(frontalStats.values.meditatorPre,frontalStats.values.controlPre,cfg.numPermutations);
 [frontalStats.tests.postGroupP,frontalStats.tests.postGroupEffect,frontalStats.tests.postGroupN] = ...
    cfGCUtils.pairedLabelSwapMatchedGroups(frontalStats.values.meditatorPost,frontalStats.values.controlPost,cfg.numPermutations);
 [frontalStats.tests.interactionP,frontalStats.tests.interactionEffect,frontalStats.tests.interactionN] = ...
    cfGCUtils.pairedLabelSwapMatchedGroups(frontalStats.values.meditatorDelta,frontalStats.values.controlDelta,cfg.numPermutations);
end

function figureFile = plotCollapsedMatrices(protocolName,roiLabels,meanMatrices,protocolFigureFolder,cfg)
hFig = cfGCUtils.createFigure(sprintf('%s collapsed connectivity matrices',protocolName),[100 100 1750 1050]);
t = tiledlayout(hFig,2,4,'TileSpacing','compact','Padding','compact');

panelTitles = { ...
    'Meditator PRE', 'Meditator POST', 'Control PRE', 'Control POST', ...
    'Meditator DELTA', 'Control DELTA', 'PRE Med-Control', 'Interaction DELTA diff'};
panelFields = { ...
    'meditatorPre','meditatorPost','controlPre','controlPost', ...
    'meditatorDelta','controlDelta','preGroupDifference','interaction'};

absoluteMax = 0;
deltaMax = 0;
for iField = 1:numel(panelFields)
    dataMatrix = meanMatrices.(panelFields{iField});
    if iField <= 4
        absoluteMax = max([absoluteMax; abs(dataMatrix(:))],[],'omitnan');
    else
        deltaMax = max([deltaMax; abs(dataMatrix(:))],[],'omitnan');
    end
end
if ~(isfinite(absoluteMax) && absoluteMax > 0); absoluteMax = 1; end
if ~(isfinite(deltaMax) && deltaMax > 0); deltaMax = 1; end

for iPanel = 1:numel(panelFields)
    ax = nexttile(t,iPanel);
    dataMatrix = meanMatrices.(panelFields{iPanel});
    imagesc(ax,dataMatrix,'AlphaData',isfinite(dataMatrix));
    set(ax,'YDir','normal');
    axis(ax,'square');
    xticks(ax,1:numel(roiLabels));
    yticks(ax,1:numel(roiLabels));
    xticklabels(ax,roiLabels);
    yticklabels(ax,roiLabels);
    xtickangle(ax,45);
    title(ax,panelTitles{iPanel},'Interpreter','none');
    colormap(ax,cfGCUtils.blueWhiteRed(256));
    colorbar(ax);
    if iPanel <= 4
        caxis(ax,[-absoluteMax absoluteMax]);
    else
        caxis(ax,[-deltaMax deltaMax]);
    end
end

cfGCUtils.addFigureTitle(hFig,sprintf('%s collapsed ROI-to-ROI connectivity matrices',protocolName));
figureFile = fullfile(protocolFigureFolder,sprintf('collapsed_connectivity_matrices.%s',cfg.figureFormat));
cfGCUtils.saveFigureQuietly(hFig,figureFile);
end

function figureFile = plotROINetFlowSummary(protocolName,roiLabels,roiNetFlow,roiStats,protocolFigureFolder,cfg)
hFig = cfGCUtils.createFigure(sprintf('%s ROI net flow summary',protocolName),[100 100 1800 1100]);
t = tiledlayout(hFig,2,1,'TileSpacing','compact','Padding','compact');

x = 1:numel(roiLabels);
deltaDifference = roiNetFlow.meditatorDelta - roiNetFlow.controlDelta;

ax1 = nexttile(t,1);
hold(ax1,'on');
plotWithError(ax1,x,roiStats.summary.meanMeditatorPre,roiStats.summary.semMeditatorPre,[0.80 0.20 0.20],'Meditator PRE');
plotWithError(ax1,x,roiStats.summary.meanMeditatorPost,roiStats.summary.semMeditatorPost,[0.55 0.00 0.00],'Meditator POST');
plotWithError(ax1,x,roiStats.summary.meanControlPre,roiStats.summary.semControlPre,[0.20 0.35 0.80],'Control PRE');
plotWithError(ax1,x,roiStats.summary.meanControlPost,roiStats.summary.semControlPost,[0.00 0.10 0.45],'Control POST');
yline(ax1,0,'k:');
set(ax1,'XTick',x,'XTickLabel',roiLabels);
xtickangle(ax1,45);
xlim(ax1,[0.5 numel(roiLabels)+0.5]);
yLimTop = axisTightFromVectors([ ...
    roiStats.summary.meanMeditatorPre + roiStats.summary.semMeditatorPre; ...
    roiStats.summary.meanMeditatorPost + roiStats.summary.semMeditatorPost; ...
    roiStats.summary.meanControlPre + roiStats.summary.semControlPre; ...
    roiStats.summary.meanControlPost + roiStats.summary.semControlPost; ...
    roiStats.summary.meanMeditatorPre - roiStats.summary.semMeditatorPre; ...
    roiStats.summary.meanMeditatorPost - roiStats.summary.semMeditatorPost; ...
    roiStats.summary.meanControlPre - roiStats.summary.semControlPre; ...
    roiStats.summary.meanControlPost - roiStats.summary.semControlPost]);
ylim(ax1,yLimTop);
legend(ax1,'Location','eastoutside');
title(ax1,'ROI net flow means +/- SEM');
ylabel(ax1,'Net flow');

sigIndicesTop = find(roiStats.interaction.sig);
if ~isempty(sigIndicesTop)
    yStar = yLimTop(2) - 0.04 * range(yLimTop);
    text(ax1,sigIndicesTop,repmat(yStar,size(sigIndicesTop)),'*','HorizontalAlignment','center', ...
        'Color',[0.60 0 0],'FontSize',12,'FontWeight','bold');
end

ax2 = nexttile(t,2);
hold(ax2,'on');
plotWithError(ax2,x,roiStats.summary.meanMeditatorDelta,roiStats.summary.semMeditatorDelta,[0.80 0.20 0.20],'Meditator DELTA');
plotWithError(ax2,x,roiStats.summary.meanControlDelta,roiStats.summary.semControlDelta,[0.20 0.35 0.80],'Control DELTA');
plotWithError(ax2,x,roiStats.summary.meanInteraction,roiStats.summary.semInteraction,[0.10 0.10 0.10],'Interaction');
yline(ax2,0,'k:');
set(ax2,'XTick',x,'XTickLabel',roiLabels);
xtickangle(ax2,45);
xlim(ax2,[0.5 numel(roiLabels)+0.5]);
yLimBottom = axisTightFromVectors([ ...
    roiStats.summary.meanMeditatorDelta + roiStats.summary.semMeditatorDelta; ...
    roiStats.summary.meanControlDelta + roiStats.summary.semControlDelta; ...
    roiStats.summary.meanInteraction + roiStats.summary.semInteraction; ...
    roiStats.summary.meanMeditatorDelta - roiStats.summary.semMeditatorDelta; ...
    roiStats.summary.meanControlDelta - roiStats.summary.semControlDelta; ...
    roiStats.summary.meanInteraction - roiStats.summary.semInteraction]);
ylim(ax2,yLimBottom);
legend(ax2,'Location','eastoutside');
title(ax2,'Pre/post change and interaction');
ylabel(ax2,'Net flow change');

sigIndicesBottom = find(roiStats.interaction.sig);
if ~isempty(sigIndicesBottom)
    yStar = yLimBottom(2) - 0.04 * range(yLimBottom);
    text(ax2,sigIndicesBottom,repmat(yStar,size(sigIndicesBottom)),'*','HorizontalAlignment','center', ...
        'Color',[0.60 0 0],'FontSize',12,'FontWeight','bold');
end

cfGCUtils.addFigureTitle(hFig,sprintf('%s ROI net flow summary (%d valid pairs)',protocolName,size(deltaDifference,2)));
figureFile = fullfile(protocolFigureFolder,sprintf('roi_net_flow_summary.%s',cfg.figureFormat));
cfGCUtils.saveFigureQuietly(hFig,figureFile);
end

function figureFile = plotBandResolvedInteraction(protocolName,roiLabels,pairInfo,meanBandInteraction,bandStats,protocolFigureFolder,cfg)
bandPairLabels = makeBandPairLabels(pairInfo);
numROIs = numel(roiLabels);
maxRows = min(cfg.maxROIsPerBandFigure,numROIs);

plotOrder = 1:maxRows;
if isfield(bandStats,'interactionMean') && ~isempty(bandStats.interactionMean)
    roiStrength = squeeze(mean(mean(abs(bandStats.interactionMean),3,'omitnan'),2,'omitnan'));
    roiStrength(~isfinite(roiStrength)) = -inf;
    [~,sortIdx] = sort(roiStrength,'descend');
    plotOrder = sortIdx(1:maxRows);
end

effectMatrix = flattenBandCube(meanBandInteraction(plotOrder,:,:));
sigMatrix = flattenBandCube(double(bandStats.interactionSig(plotOrder,:,:)));

hFig = cfGCUtils.createFigure(sprintf('%s band-resolved interaction',protocolName),[100 100 1800 900]);
t = tiledlayout(hFig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(t,1);
imagesc(ax1,effectMatrix,'AlphaData',isfinite(effectMatrix));
set(ax1,'YDir','normal');
xticks(ax1,1:numel(bandPairLabels));
yticks(ax1,1:numel(plotOrder));
xticklabels(ax1,bandPairLabels);
yticklabels(ax1,roiLabels(plotOrder));
xtickangle(ax1,45);
title(ax1,'Mean interaction net flow by ROI and band pair');
colorbar(ax1);
colormap(ax1,cfGCUtils.blueWhiteRed(256));
bandMax = max(abs(effectMatrix(:)),[],'omitnan');
if ~(isfinite(bandMax) && bandMax > 0); bandMax = 1; end
caxis(ax1,[-bandMax bandMax]);

sigLocations = find(sigMatrix > 0);
if ~isempty(sigLocations)
    [sigRows,sigCols] = ind2sub(size(sigMatrix),sigLocations);
    hold(ax1,'on');
    plot(ax1,sigCols,sigRows,'k.','MarkerSize',8);
end

ax2 = nexttile(t,2);
imagesc(ax2,sigMatrix,'AlphaData',isfinite(sigMatrix));
set(ax2,'YDir','normal');
xticks(ax2,1:numel(bandPairLabels));
yticks(ax2,1:numel(plotOrder));
xticklabels(ax2,bandPairLabels);
yticklabels(ax2,roiLabels(plotOrder));
xtickangle(ax2,45);
title(ax2,'FDR-significant interaction cells (1 = significant)');
colorbar(ax2);
colormap(ax2,parula(64));
caxis(ax2,[0 1]);

cfGCUtils.addFigureTitle(hFig,sprintf('%s band-resolved interaction net flow',protocolName));
figureFile = fullfile(protocolFigureFolder,sprintf('band_resolved_interaction.%s',cfg.figureFormat));
cfGCUtils.saveFigureQuietly(hFig,figureFile);
end

function figureFile = plotFrontalDominance(protocolName,frontalStats,protocolFigureFolder,cfg)
figureFile = fullfile(protocolFigureFolder,sprintf('frontal_dominance.%s',cfg.figureFormat));

if ~frontalStats.available
    hFig = cfGCUtils.createFigure(sprintf('%s frontal dominance unavailable',protocolName),[100 100 1000 450]);
    ax = axes(hFig);
    axis(ax,'off');
    text(ax,0.5,0.5,sprintf('No frontal/non-frontal split available for %s ROI labels.',protocolName), ...
        'HorizontalAlignment','center','FontSize',14);
    cfGCUtils.addFigureTitle(hFig,sprintf('%s frontal dominance summary',protocolName));
    cfGCUtils.saveFigureQuietly(hFig,figureFile);
    return;
end

vals = frontalStats.values;
tests = frontalStats.tests;

means = [ ...
    mean(vals.meditatorPre,'omitnan'), mean(vals.meditatorPost,'omitnan'); ...
    mean(vals.controlPre,'omitnan'), mean(vals.controlPost,'omitnan'); ...
    mean(vals.meditatorDelta,'omitnan'), mean(vals.controlDelta,'omitnan')];
sems = [ ...
    computeSEM(vals.meditatorPre,2), computeSEM(vals.meditatorPost,2); ...
    computeSEM(vals.controlPre,2), computeSEM(vals.controlPost,2); ...
    computeSEM(vals.meditatorDelta,2), computeSEM(vals.controlDelta,2)];

labels = {'Meditator dominance','Control dominance','Delta dominance'};
conditionLabels = {'State 1','State 2'};

hFig = cfGCUtils.createFigure(sprintf('%s frontal dominance',protocolName),[100 100 1500 500]);
t = tiledlayout(hFig,1,3,'TileSpacing','compact','Padding','compact');

for iPanel = 1:3
    ax = nexttile(t,iPanel);
    bar(ax,1:2,means(iPanel,:),'FaceColor',[0.75 0.75 0.75]);
    hold(ax,'on');
    errorbar(ax,1:2,means(iPanel,:),sems(iPanel,:),'k.','LineWidth',1.2);
    yline(ax,0,'k:');
    set(ax,'XTick',1:2,'XTickLabel',conditionLabels);
    title(ax,labels{iPanel});
    ylabel(ax,'Frontal minus non-frontal net flow');
end

annotationText = sprintf(['Meditator pre/post p = %.4g | Control pre/post p = %.4g | ' ...
    'PRE group p = %.4g | POST group p = %.4g | Interaction p = %.4g'], ...
    tests.meditatorChangeP,tests.controlChangeP,tests.baselinePreGroupP,tests.postGroupP,tests.interactionP);
annotation(hFig,'textbox',[0.03 0.01 0.94 0.07],'String',annotationText,'EdgeColor','none', ...
    'HorizontalAlignment','center','Interpreter','none');
cfGCUtils.addFigureTitle(hFig,sprintf('%s frontal dominance summary',protocolName));
cfGCUtils.saveFigureQuietly(hFig,figureFile);
end

function overviewFigureFiles = createCrossProtocolFigures(netFlowSummary,cfg)
overviewFigureFiles = struct;
overviewFigureFiles.protocolByROI = plotCrossProtocolROIOverview(netFlowSummary,cfg);
overviewFigureFiles.protocolByBand = plotCrossProtocolBandOverview(netFlowSummary,cfg);
overviewFigureFiles.frontalDominance = plotCrossProtocolFrontalOverview(netFlowSummary,cfg);
end

function figureFile = plotCrossProtocolROIOverview(netFlowSummary,cfg)
protocolResults = netFlowSummary.protocolResults;
numProtocols = numel(protocolResults);
numROIs = numel(netFlowSummary.roiLabels);

effectMatrix = nan(numProtocols,numROIs);
qMatrix = nan(numProtocols,numROIs);
for iProtocol = 1:numProtocols
    effectMatrix(iProtocol,:) = protocolResults(iProtocol).roiStats.interaction.effect;
    qMatrix(iProtocol,:) = protocolResults(iProtocol).roiStats.interaction.q;
end

hFig = cfGCUtils.createFigure('Cross-protocol ROI interaction overview',[100 100 1600 700]);
t = tiledlayout(hFig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(t,1);
imagesc(ax1,effectMatrix,'AlphaData',isfinite(effectMatrix));
set(ax1,'YDir','normal');
xticks(ax1,1:numROIs);
yticks(ax1,1:numProtocols);
xticklabels(ax1,netFlowSummary.roiLabels);
yticklabels(ax1,netFlowSummary.protocolNameList);
xtickangle(ax1,45);
title(ax1,'Interaction effect size by protocol and ROI');
colorbar(ax1);
colormap(ax1,cfGCUtils.blueWhiteRed(256));
maxVal = max(abs(effectMatrix(:)),[],'omitnan');
if ~(isfinite(maxVal) && maxVal > 0); maxVal = 1; end
caxis(ax1,[-maxVal maxVal]);

sigMask = qMatrix <= cfg.fdrAlpha;
sigLocations = find(sigMask);
if ~isempty(sigLocations)
    [sigRows,sigCols] = ind2sub(size(sigMask),sigLocations);
    hold(ax1,'on');
    plot(ax1,sigCols,sigRows,'k.','MarkerSize',10);
end

ax2 = nexttile(t,2);
imagesc(ax2,-log10(qMatrix),'AlphaData',isfinite(qMatrix));
set(ax2,'YDir','normal');
xticks(ax2,1:numROIs);
yticks(ax2,1:numProtocols);
xticklabels(ax2,netFlowSummary.roiLabels);
yticklabels(ax2,netFlowSummary.protocolNameList);
xtickangle(ax2,45);
title(ax2,'-log10(FDR q) for ROI interaction');
colorbar(ax2);

cfGCUtils.addFigureTitle(hFig,'Cross-protocol ROI interaction overview');
figureFile = fullfile(cfg.figureFolder,sprintf('overview_protocol_by_roi_interaction.%s',cfg.figureFormat));
cfGCUtils.saveFigureQuietly(hFig,figureFile);
end

function figureFile = plotCrossProtocolBandOverview(netFlowSummary,cfg)
protocolResults = netFlowSummary.protocolResults;
numProtocols = numel(protocolResults);
bandPairLabels = netFlowSummary.bandPairLabels;
numBandPairs = numel(bandPairLabels);

effectMatrix = nan(numProtocols,numBandPairs);
for iProtocol = 1:numProtocols
    overallInteraction = protocolResults(iProtocol).bandStats.overallBandInteraction;
    effectMatrix(iProtocol,:) = flattenBandCube(overallInteraction);
end

hFig = cfGCUtils.createFigure('Cross-protocol band interaction overview',[100 100 1650 600]);
ax = axes(hFig);
imagesc(ax,effectMatrix,'AlphaData',isfinite(effectMatrix));
set(ax,'YDir','normal');
xticks(ax,1:numBandPairs);
yticks(ax,1:numProtocols);
xticklabels(ax,bandPairLabels);
yticklabels(ax,netFlowSummary.protocolNameList);
xtickangle(ax,45);
title(ax,'Mean interaction net flow averaged across ROIs');
colorbar(ax);
colormap(ax,cfGCUtils.blueWhiteRed(256));
maxVal = max(abs(effectMatrix(:)),[],'omitnan');
if ~(isfinite(maxVal) && maxVal > 0); maxVal = 1; end
caxis(ax,[-maxVal maxVal]);

cfGCUtils.addFigureTitle(hFig,'Cross-protocol band-pair interaction overview');
figureFile = fullfile(cfg.figureFolder,sprintf('overview_protocol_by_band_interaction.%s',cfg.figureFormat));
cfGCUtils.saveFigureQuietly(hFig,figureFile);
end

function figureFile = plotCrossProtocolFrontalOverview(netFlowSummary,cfg)
protocolResults = netFlowSummary.protocolResults;
protocolNames = netFlowSummary.protocolNameList;
numProtocols = numel(protocolResults);
interactionValues = nan(numProtocols,1);
interactionP = nan(numProtocols,1);

for iProtocol = 1:numProtocols
    if protocolResults(iProtocol).frontalStats.available
        interactionValues(iProtocol) = mean(protocolResults(iProtocol).frontalStats.values.meditatorDelta - ...
            protocolResults(iProtocol).frontalStats.values.controlDelta,'omitnan');
        interactionP(iProtocol) = protocolResults(iProtocol).frontalStats.tests.interactionP;
    end
end

hFig = cfGCUtils.createFigure('Cross-protocol frontal dominance overview',[100 100 1300 500]);
ax = axes(hFig);
bar(ax,interactionValues,'FaceColor',[0.55 0.55 0.55]);
hold(ax,'on');
yline(ax,0,'k:');
set(ax,'XTick',1:numProtocols,'XTickLabel',protocolNames);
xtickangle(ax,45);
ylabel(ax,'Mean frontal dominance interaction');
title(ax,'Frontal minus non-frontal dominance interaction by protocol');

sigIdx = find(interactionP <= cfg.fdrAlpha);
if ~isempty(sigIdx)
    yLim = ylim(ax);
    yText = yLim(2) - 0.05 * range(yLim);
    text(ax,sigIdx,repmat(yText,size(sigIdx)),'*','HorizontalAlignment','center', ...
        'FontSize',12,'FontWeight','bold','Color',[0.6 0 0]);
end

cfGCUtils.addFigureTitle(hFig,'Cross-protocol frontal dominance overview');
figureFile = fullfile(cfg.figureFolder,sprintf('overview_frontal_dominance.%s',cfg.figureFormat));
cfGCUtils.saveFigureQuietly(hFig,figureFile);
end

function collapsed = collapseBandGrid(gridData)
collapsed = squeeze(mean(mean(gridData,4,'omitnan'),3,'omitnan'));
end

function bandNetFlow = computeBandResolvedNetFlow(gridData)
[numROIs,~,numPhaseBands,numAmpBands] = size(gridData);
bandNetFlow = nan(numROIs,numPhaseBands,numAmpBands);
for iPhase = 1:numPhaseBands
    for iAmp = 1:numAmpBands
        bandNetFlow(:,iPhase,iAmp) = computeNetFlowVector(gridData(:,:,iPhase,iAmp));
    end
end
end

function netFlowVector = computeNetFlowVector(matrixData)
outFlow = mean(matrixData,2,'omitnan');
inFlow = mean(matrixData,1,'omitnan')';
netFlowVector = outFlow - inFlow;
end

function semVal = computeSEM(dataArray,dim)
if ~exist('dim','var') || isempty(dim)
    dim = 1;
end
nUsed = sum(isfinite(dataArray),dim);
semVal = std(dataArray,0,dim,'omitnan') ./ sqrt(max(nUsed,1));
semVal(nUsed == 0) = NaN;
end

function plotWithError(ax,x,y,semVec,colorValue,labelText)
errorbar(ax,x,y,semVec,'-o', ...
    'Color',colorValue, ...
    'MarkerFaceColor',colorValue, ...
    'MarkerEdgeColor',colorValue, ...
    'LineWidth',1.2, ...
    'MarkerSize',4, ...
    'DisplayName',labelText);
end

function yLimits = axisTightFromVectors(values)
finiteVals = values(isfinite(values));
if isempty(finiteVals)
    yLimits = [-1 1];
    return;
end
minVal = min(finiteVals);
maxVal = max(finiteVals);
if minVal == maxVal
    pad = max(1e-3,abs(minVal) * 0.2 + 1e-3);
else
    pad = 0.12 * (maxVal - minVal);
end
yLimits = [minVal - pad, maxVal + pad];
end

function bandPairLabels = makeBandPairLabels(pairInfo)
numPhaseBands = size(pairInfo,1);
numAmpBands = size(pairInfo,2);
bandPairLabels = cell(1,numPhaseBands * numAmpBands);
idx = 0;
for iPhase = 1:numPhaseBands
    for iAmp = 1:numAmpBands
        idx = idx + 1;
        phaseRange = pairInfo(iPhase,iAmp).phaseRange;
        ampRange = pairInfo(iPhase,iAmp).ampRange;
        bandPairLabels{idx} = sprintf('P %.1f-%.1f | A %.1f-%.1f',phaseRange(1),phaseRange(2),ampRange(1),ampRange(2));
    end
end
end

function flatMatrix = flattenBandCube(dataCube)
if ismatrix(dataCube)
    flatMatrix = reshape(dataCube,1,[]);
    return;
end
numRows = size(dataCube,1);
flatMatrix = reshape(dataCube,numRows,[]);
end

function [qVals,sigMask] = applyBHFDR(pVals,alpha)
pVals = pVals(:);
qVals = nan(size(pVals));
sigMask = false(size(pVals));

validMask = isfinite(pVals);
if ~any(validMask)
    return;
end

validP = pVals(validMask);
[sortedP,sortIdx] = sort(validP);
m = numel(sortedP);
qSorted = sortedP .* m ./ (1:m)';
qSorted = flipud(cummin(flipud(qSorted)));
qSorted = min(qSorted,1);

sigSorted = sortedP <= ((1:m)'/m) * alpha;
if any(sigSorted)
    lastSig = find(sigSorted,1,'last');
    threshold = sortedP(lastSig);
    sigSorted = sortedP <= threshold;
else
    sigSorted = false(size(sortedP));
end

restoreIdx = zeros(size(sortIdx));
restoreIdx(sortIdx) = 1:m;
qVals(validMask) = qSorted(restoreIdx);
sigMask(validMask) = sigSorted(restoreIdx);
end

function checkpointFile = getProtocolCheckpointFile(protocolName,cfg)
checkpointFile = fullfile(cfg.checkpointFolder,sprintf('protocol_%s_checkpoint.mat',sanitizeToken(protocolName)));
end

function saveRunCheckpoint(netFlowSummary,overviewFigureFiles,cfg,statusText)
checkpointData = struct;
checkpointData.savedOn = datetime('now');
checkpointData.status = statusText;
checkpointData.netFlowSummary = netFlowSummary;
checkpointData.overviewFigureFiles = overviewFigureFiles;
checkpointData.config = cfg;
save(fullfile(cfg.checkpointFolder,'net_flow_run_checkpoint.mat'),'checkpointData','-v7.3');
end

function token = sanitizeToken(textValue)
token = regexprep(char(string(textValue)),'[^a-zA-Z0-9_-]','_');
end

function s = setDefault(s,fieldName,defaultValue)
if ~isfield(s,fieldName) || isempty(s.(fieldName))
    s.(fieldName) = defaultValue;
end
end
