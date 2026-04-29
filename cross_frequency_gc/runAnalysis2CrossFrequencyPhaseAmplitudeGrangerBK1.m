function analysis2Result = runAnalysis2CrossFrequencyPhaseAmplitudeGrangerBK1(pairSelection,protocolNameList,badEyeCondition,badTrialVersion,cfgIn)
% runAnalysis2CrossFrequencyPhaseAmplitudeGrangerBK1
%
% Run secondary summary analyses on saved CF-GC comparison outputs.
% This layer creates compact MAT summaries without generating the large
% figure sets used by the primary visualization/comparison workflows.
%
% It performs:
%   1) Protocol consistency summaries across the 8 protocols
%   2) Effect-size aggregation summaries
%   3) Raw-vs-normalized inference comparison
%   4) Pair-level heterogeneity summaries
%   5) Global net-flow summaries
%
% Results are saved as five MAT files under <cf_gc>/analysis_2/.
%
% Example:
%   runAnalysis2CrossFrequencyPhaseAmplitudeGrangerBK1('all','all','ep','v8')

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

cfGCUtils.ensureFolder(cfg.analysisFolder);
cfGCUtils.ensureFolder(cfg.checkpointFolder);
preparedData = prepareProtocolData(pairList,protocolNameList,badEyeCondition,badTrialVersion,cfg);

normalizedInference = computeInferenceAcrossProtocols(preparedData,'normalized',cfg);
rawInference = computeInferenceAcrossProtocols(preparedData,'raw',cfg);

protocolConsistency = buildProtocolConsistencySummary(preparedData,normalizedInference,cfg);
effectSizeSummary = buildEffectSizeSummary(preparedData,normalizedInference,rawInference,cfg);
rawVsNormalized = buildRawVsNormalizedSummary(preparedData,normalizedInference,rawInference,cfg);
pairHeterogeneity = buildPairHeterogeneitySummary(preparedData,cfg);
globalNetFlow = buildGlobalNetFlowSummary(preparedData,cfg);

protocolConsistencyFile = fullfile(cfg.analysisFolder,'protocol_consistency_summary.mat');
effectSizeFile = fullfile(cfg.analysisFolder,'effect_size_summary.mat');
rawVsNormalizedFile = fullfile(cfg.analysisFolder,'raw_vs_normalized_inference.mat');
pairHeterogeneityFile = fullfile(cfg.analysisFolder,'pair_heterogeneity_summary.mat');
globalNetFlowFile = fullfile(cfg.analysisFolder,'global_net_flow_summary.mat');

save(protocolConsistencyFile,'protocolConsistency','-v7.3');
save(effectSizeFile,'effectSizeSummary','-v7.3');
save(rawVsNormalizedFile,'rawVsNormalized','-v7.3');
save(pairHeterogeneityFile,'pairHeterogeneity','-v7.3');
save(globalNetFlowFile,'globalNetFlow','-v7.3');

analysis2Result = struct;
analysis2Result.analysisFolder = cfg.analysisFolder;
analysis2Result.selectionLabel = selectionLabel;
analysis2Result.badEyeCondition = badEyeCondition;
analysis2Result.badTrialVersion = badTrialVersion;
analysis2Result.protocolNameList = protocolNameList;
analysis2Result.numRequestedPairs = numel(pairList);
analysis2Result.outputFiles = {protocolConsistencyFile; effectSizeFile; rawVsNormalizedFile; pairHeterogeneityFile; globalNetFlowFile};
analysis2Result.config = cfg;

fprintf('Analysis-2 completed. Saved 5 MAT files under:\n%s\n',cfg.analysisFolder);
end

function cfg = applyDefaultCfg(cfgIn)
thisFolder = fileparts(mfilename('fullpath'));

cfg = cfgIn;
cfg = setDefault(cfg,'dataFolder',fullfile(thisFolder,'savedDataCrossFreqGranger'));
cfg = setDefault(cfg,'analysisFolder',fullfile(thisFolder,'analysis_2'));
cfg = setDefault(cfg,'checkpointFolder',fullfile(thisFolder,'analysis_2','checkpoints'));
cfg = setDefault(cfg,'numPermutations',5000);
cfg = setDefault(cfg,'minObservationCount',20);
cfg = setDefault(cfg,'minValidSubjectsPerCell',3);
cfg = setDefault(cfg,'fdrAlpha',0.05);
cfg = setDefault(cfg,'consistencyPThreshold',0.05);
cfg = setDefault(cfg,'moderateEffectThreshold',0.50);
cfg = setDefault(cfg,'strongEffectThreshold',0.80);
cfg = setDefault(cfg,'maxRankedCells',50);
cfg = setDefault(cfg,'randomSeed',1);
cfg = setDefault(cfg,'resumeIfAvailable',1);
end

function preparedData = prepareProtocolData(pairList,protocolNameList,badEyeCondition,badTrialVersion,cfg)
preparedData = struct;
preparedData.createdOn = datetime('now');
preparedData.protocolNameList = protocolNameList;
preparedData.badEyeCondition = badEyeCondition;
preparedData.badTrialVersion = badTrialVersion;
preparedData.requestedPairList = pairList;
preparedData.protocolResults = [];

for iProtocol = 1:numel(protocolNameList)
    protocolName = protocolNameList{iProtocol};
    checkpointFile = fullfile(cfg.checkpointFolder,sprintf('prepared_protocol_%s.mat',protocolName));
    if cfg.resumeIfAvailable && exist(checkpointFile,'file')
        tmp = load(checkpointFile,'protocolData');
        protocolData = tmp.protocolData;
        fprintf('Loaded preparation checkpoint for protocol %s\n',protocolName);
    else
        protocolData = prepareSingleProtocol(pairList,protocolName,badEyeCondition,badTrialVersion,cfg);
        save(checkpointFile,'protocolData','-v7.3');
        fprintf('Saved preparation checkpoint for protocol %s\n',protocolName);
    end
    if isempty(preparedData.protocolResults)
        preparedData.protocolResults = repmat(protocolData,1,numel(protocolNameList));
    end
    preparedData.protocolResults(iProtocol) = protocolData;
end

preparedData.roiLabels = preparedData.protocolResults(1).roiLabels;
preparedData.pairInfo = preparedData.protocolResults(1).pairInfo;
preparedData.phaseBandList = preparedData.protocolResults(1).phaseBandList;
preparedData.ampBandList = preparedData.protocolResults(1).ampBandList;
preparedData.familyDefinitions = getFamilyDefinitions();
preparedData.normalizationInfo = struct( ...
    'functionName','cfGCUtils.normalizeSubjectGrid', ...
    'note','Analysis-2 reuses the same normalization helper as the main comparison pipeline to keep normalized inference directly comparable.');
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

        normalized.preMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        normalized.postMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        normalized.preCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        normalized.postCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        normalized.deltaMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        normalized.deltaCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);

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

validIndices = find(validPairMask);
protocolData = struct;
protocolData.protocolName = protocolName;
protocolData.requestedPairList = pairList;
protocolData.validPairMask = validPairMask;
protocolData.validPairList = pairList(validIndices);
protocolData.pairAvailability = pairAvailability;
protocolData.numRequestedPairs = numPairsRequested;
protocolData.numValidPairs = numel(validIndices);
protocolData.validIndices = validIndices;
protocolData.roiLabels = templateResults.roiLabels;
protocolData.pairInfo = pairInfo;
protocolData.phaseBandList = phaseBandList;
protocolData.ampBandList = ampBandList;
protocolData.raw = restrictToValidPairs(raw,validIndices);
protocolData.normalized = restrictToValidPairs(normalized,validIndices);
protocolData.sampleCount = restrictToValidPairs(sampleCount,validIndices);
protocolData.validPairLabels = {protocolData.validPairList.pairLabel};
protocolData.validMeditatorNames = {protocolData.validPairList.meditator};
protocolData.validControlNames = {protocolData.validPairList.control};
end

function out = restrictToValidPairs(in,validIndices)
fieldNames = fieldnames(in);
out = struct;
for iField = 1:numel(fieldNames)
    out.(fieldNames{iField}) = in.(fieldNames{iField})(:,:,:,:,validIndices);
end
end

function inferenceSummary = computeInferenceAcrossProtocols(preparedData,spaceName,cfg)
protocolResults = preparedData.protocolResults;
summaryCheckpointFile = fullfile(cfg.checkpointFolder,sprintf('%s_inference_summary.mat',spaceName));
if cfg.resumeIfAvailable && exist(summaryCheckpointFile,'file')
    tmp = load(summaryCheckpointFile,'inferenceSummary');
    inferenceSummary = tmp.inferenceSummary;
    if numel(inferenceSummary.protocolResults) == numel(protocolResults)
        fprintf('Loaded full %s inference summary checkpoint\n',spaceName);
        return;
    end
end

inferenceSummary = struct;
inferenceSummary.spaceName = spaceName;
inferenceSummary.protocolNameList = preparedData.protocolNameList;
inferenceSummary.familyDefinitions = preparedData.familyDefinitions;
inferenceSummary.protocolResults = [];

for iProtocol = 1:numel(protocolResults)
    protocolData = protocolResults(iProtocol);
    checkpointFile = fullfile(cfg.checkpointFolder,sprintf('%s_inference_%s.mat',spaceName,protocolData.protocolName));
    if cfg.resumeIfAvailable && exist(checkpointFile,'file')
        tmp = load(checkpointFile,'result');
        result = tmp.result;
        fprintf('Loaded %s inference checkpoint for protocol %s\n',spaceName,protocolData.protocolName);
    else
        inference = computeCellwiseInference(protocolData.(spaceName),cfg);

        result = struct;
        result.protocolName = protocolData.protocolName;
        result.numValidPairs = protocolData.numValidPairs;
        result.validPairLabels = protocolData.validPairLabels;
        result.inference = inference;
        save(checkpointFile,'result','-v7.3');
        fprintf('Saved %s inference checkpoint for protocol %s\n',spaceName,protocolData.protocolName);
    end
    if isempty(inferenceSummary.protocolResults)
        inferenceSummary.protocolResults = repmat(result,1,numel(protocolResults));
    end
    inferenceSummary.protocolResults(iProtocol) = result;
    save(summaryCheckpointFile,'inferenceSummary','-v7.3');
end
end

function inference = computeCellwiseInference(dataStruct,cfg)
dataSize = size(dataStruct.preMed);
numROIs = dataSize(1);
numPhase = dataSize(3);
numAmp = dataSize(4);

familyDefs = getFamilyDefinitions();
inference = struct;

for iFamily = 1:numel(familyDefs)
    inference.(familyDefs(iFamily).pField) = nan(numROIs,numROIs,numPhase,numAmp);
    inference.(familyDefs(iFamily).effectField) = nan(numROIs,numROIs,numPhase,numAmp);
    inference.(familyDefs(iFamily).nField) = nan(numROIs,numROIs,numPhase,numAmp);
end

for iROI = 1:numROIs
    for jROI = 1:numROIs
        for iPhase = 1:numPhase
            for iAmp = 1:numAmp
                preMed = squeeze(dataStruct.preMed(iROI,jROI,iPhase,iAmp,:));
                postMed = squeeze(dataStruct.postMed(iROI,jROI,iPhase,iAmp,:));
                preCtrl = squeeze(dataStruct.preCtrl(iROI,jROI,iPhase,iAmp,:));
                postCtrl = squeeze(dataStruct.postCtrl(iROI,jROI,iPhase,iAmp,:));
                deltaMed = postMed - preMed;
                deltaCtrl = postCtrl - preCtrl;

                variableSet = struct;
                variableSet.preMed = preMed;
                variableSet.postMed = postMed;
                variableSet.preCtrl = preCtrl;
                variableSet.postCtrl = postCtrl;
                variableSet.deltaMed = deltaMed;
                variableSet.deltaCtrl = deltaCtrl;

                for iFamily = 1:numel(familyDefs)
                    [pVal,effectVal,nUsed] = runFamilyTest(variableSet,familyDefs(iFamily),cfg.numPermutations);
                    inference.(familyDefs(iFamily).pField)(iROI,jROI,iPhase,iAmp) = pVal;
                    inference.(familyDefs(iFamily).effectField)(iROI,jROI,iPhase,iAmp) = effectVal;
                    inference.(familyDefs(iFamily).nField)(iROI,jROI,iPhase,iAmp) = nUsed;
                end
            end
        end
    end
end

inference = validateInferenceCompleteness(inference,familyDefs);
inference = applyMinimumNMask(inference,cfg.minValidSubjectsPerCell);
for iFamily = 1:numel(familyDefs)
    [inference.(familyDefs(iFamily).qField),inference.(familyDefs(iFamily).sigField)] = ...
        cfGCUtils.bhFDRByPair(inference.(familyDefs(iFamily).pField),cfg.fdrAlpha);
end
end

function inference = applyMinimumNMask(inference,minN)
familyDefs = getFamilyDefinitions();
for iFamily = 1:numel(familyDefs)
    if ~isfield(inference,familyDefs(iFamily).nField)
        error('Inference field missing required N-field: %s',familyDefs(iFamily).nField);
    end
    invalidMask = inference.(familyDefs(iFamily).nField) < minN;
    inference.(familyDefs(iFamily).pField)(invalidMask) = NaN;
    inference.(familyDefs(iFamily).effectField)(invalidMask) = NaN;
end
end

function [pVal,effectVal,nUsed] = runFamilyTest(variableSet,familyDef,numPermutations)
switch familyDef.testKind
    case 'paired_prepost'
        xValues = variableSet.(familyDef.xField);
        yValues = variableSet.(familyDef.yField);
        [pVal,effectVal,nUsed] = cfGCUtils.pairedLabelSwapPrePost(xValues,yValues,numPermutations);

    case 'matched_groups'
        xValues = variableSet.(familyDef.xField);
        yValues = variableSet.(familyDef.yField);
        [pVal,effectVal,nUsed] = cfGCUtils.pairedLabelSwapMatchedGroups(xValues,yValues,numPermutations);

    otherwise
        error('Unknown family test kind: %s',familyDef.testKind);
end
end

function inference = validateInferenceCompleteness(inference,familyDefs)
for iFamily = 1:numel(familyDefs)
    requiredFields = {familyDefs(iFamily).pField, familyDefs(iFamily).effectField, familyDefs(iFamily).nField};
    for iField = 1:numel(requiredFields)
        if ~isfield(inference,requiredFields{iField})
            error('Inference is missing required field: %s',requiredFields{iField});
        end
    end

    nValues = inference.(familyDefs(iFamily).nField);
    pValues = inference.(familyDefs(iFamily).pField);
    effectValues = inference.(familyDefs(iFamily).effectField);
    hasUsableSamples = any(nValues(:) >= 2);
    hasAnyFiniteResult = any(isfinite(pValues(:))) || any(isfinite(effectValues(:)));

    if hasUsableSamples && ~hasAnyFiniteResult
        error('Inference family %s appears initialized but never populated.',familyDefs(iFamily).name);
    end
end
end

function protocolConsistency = buildProtocolConsistencySummary(preparedData,normalizedInference,cfg)
familyDefs = normalizedInference.familyDefinitions;
protocolConsistency = struct;
protocolConsistency.analysisType = 'protocol_consistency';
protocolConsistency.protocolNameList = preparedData.protocolNameList;
protocolConsistency.roiLabels = preparedData.roiLabels;
protocolConsistency.pairInfo = preparedData.pairInfo;
protocolConsistency.familySummaries = [];

for iFamily = 1:numel(familyDefs)
    pStack = stackInferenceField(normalizedInference.protocolResults,familyDefs(iFamily).pField);
    qStack = stackInferenceField(normalizedInference.protocolResults,familyDefs(iFamily).qField);
    effectStack = stackInferenceField(normalizedInference.protocolResults,familyDefs(iFamily).effectField);

    effectSign = sign(effectStack);
    meanEffect = mean(effectStack,5,'omitnan');
    meanAbsEffect = mean(abs(effectStack),5,'omitnan');
    medianEffect = median(effectStack,5,'omitnan');
    pCount = sum(pStack < cfg.consistencyPThreshold,5,'omitnan');
    qCount = sum(qStack < cfg.fdrAlpha,5,'omitnan');
    moderateCount = sum(abs(effectStack) >= cfg.moderateEffectThreshold,5,'omitnan');
    strongCount = sum(abs(effectStack) >= cfg.strongEffectThreshold,5,'omitnan');
    signAgreementFraction = computeSignAgreementFraction(effectSign);

    familySummary = struct;
    familySummary.familyName = familyDefs(iFamily).name;
    familySummary.pField = familyDefs(iFamily).pField;
    familySummary.effectField = familyDefs(iFamily).effectField;
    familySummary.pStack = pStack;
    familySummary.qStack = qStack;
    familySummary.effectStack = effectStack;
    familySummary.protocolCountP05 = pCount;
    familySummary.protocolCountQ05 = qCount;
    familySummary.protocolCountModerateEffect = moderateCount;
    familySummary.protocolCountStrongEffect = strongCount;
    familySummary.meanEffect = meanEffect;
    familySummary.medianEffect = medianEffect;
    familySummary.meanAbsEffect = meanAbsEffect;
    familySummary.signAgreementFraction = signAgreementFraction;
    familySummary.topRankedCells = rankCellsByConsistency(preparedData,pStack,qStack,effectStack,meanEffect,meanAbsEffect,signAgreementFraction,cfg);

    if isempty(protocolConsistency.familySummaries)
        protocolConsistency.familySummaries = repmat(familySummary,1,numel(familyDefs));
    end
    protocolConsistency.familySummaries(iFamily) = familySummary;
end
end

function effectSizeSummary = buildEffectSizeSummary(preparedData,normalizedInference,rawInference,cfg)
familyDefs = normalizedInference.familyDefinitions;
effectSizeSummary = struct;
effectSizeSummary.analysisType = 'effect_size_summary';
effectSizeSummary.protocolNameList = preparedData.protocolNameList;
effectSizeSummary.roiLabels = preparedData.roiLabels;
effectSizeSummary.pairInfo = preparedData.pairInfo;
effectSizeSummary.normalized = summarizeEffectStacks(preparedData,normalizedInference,familyDefs,cfg);
effectSizeSummary.raw = summarizeEffectStacks(preparedData,rawInference,familyDefs,cfg);
end

function summary = summarizeEffectStacks(preparedData,inferenceSummary,familyDefs,cfg)
summary = struct;
summary.spaceName = inferenceSummary.spaceName;
summary.familySummaries = [];

for iFamily = 1:numel(familyDefs)
    effectStack = stackInferenceField(inferenceSummary.protocolResults,familyDefs(iFamily).effectField);
    meanEffect = mean(effectStack,5,'omitnan');
    medianEffect = median(effectStack,5,'omitnan');
    meanAbsEffect = mean(abs(effectStack),5,'omitnan');
    signAgreementFraction = computeSignAgreementFraction(sign(effectStack));
    moderateCount = sum(abs(effectStack) >= cfg.moderateEffectThreshold,5,'omitnan');
    strongCount = sum(abs(effectStack) >= cfg.strongEffectThreshold,5,'omitnan');
    consistencyWeightedEffect = meanAbsEffect .* signAgreementFraction;

    familySummary = struct;
    familySummary.familyName = familyDefs(iFamily).name;
    familySummary.effectField = familyDefs(iFamily).effectField;
    familySummary.effectStack = effectStack;
    familySummary.meanEffect = meanEffect;
    familySummary.medianEffect = medianEffect;
    familySummary.meanAbsEffect = meanAbsEffect;
    familySummary.signAgreementFraction = signAgreementFraction;
    familySummary.protocolCountModerateEffect = moderateCount;
    familySummary.protocolCountStrongEffect = strongCount;
    familySummary.consistencyWeightedEffect = consistencyWeightedEffect;
    familySummary.topRankedCells = rankCellsByEffect(preparedData,effectStack,consistencyWeightedEffect,meanEffect,meanAbsEffect,signAgreementFraction,cfg);

    if isempty(summary.familySummaries)
        summary.familySummaries = repmat(familySummary,1,numel(familyDefs));
    end
    summary.familySummaries(iFamily) = familySummary;
end
end

function rawVsNormalized = buildRawVsNormalizedSummary(preparedData,normalizedInference,rawInference,cfg)
familyDefs = normalizedInference.familyDefinitions;
rawVsNormalized = struct;
rawVsNormalized.analysisType = 'raw_vs_normalized';
rawVsNormalized.protocolNameList = preparedData.protocolNameList;
rawVsNormalized.roiLabels = preparedData.roiLabels;
rawVsNormalized.pairInfo = preparedData.pairInfo;
rawVsNormalized.familySummaries = [];

for iFamily = 1:numel(familyDefs)
    normalizedEffect = stackInferenceField(normalizedInference.protocolResults,familyDefs(iFamily).effectField);
    rawEffect = stackInferenceField(rawInference.protocolResults,familyDefs(iFamily).effectField);
    normalizedP = stackInferenceField(normalizedInference.protocolResults,familyDefs(iFamily).pField);
    rawP = stackInferenceField(rawInference.protocolResults,familyDefs(iFamily).pField);
    normalizedQ = stackInferenceField(normalizedInference.protocolResults,familyDefs(iFamily).qField);
    rawQ = stackInferenceField(rawInference.protocolResults,familyDefs(iFamily).qField);

    familySummary = struct;
    familySummary.familyName = familyDefs(iFamily).name;
    familySummary.effectDifference = rawEffect - normalizedEffect;
    familySummary.pDifference = rawP - normalizedP;
    familySummary.qDifference = rawQ - normalizedQ;
    familySummary.protocolComparisons = [];

    for iProtocol = 1:numel(preparedData.protocolNameList)
        rawEffectThis = rawEffect(:,:,:,:,iProtocol);
        normEffectThis = normalizedEffect(:,:,:,:,iProtocol);
        rawPThis = rawP(:,:,:,:,iProtocol);
        normPThis = normalizedP(:,:,:,:,iProtocol);
        rawQThis = rawQ(:,:,:,:,iProtocol);
        normQThis = normalizedQ(:,:,:,:,iProtocol);

        protocolComparison = struct;
        protocolComparison.protocolName = preparedData.protocolNameList{iProtocol};
        protocolComparison.effectCorrelation = safeVectorCorrelation(rawEffectThis,normEffectThis);
        protocolComparison.pCorrelation = safeVectorCorrelation(rawPThis,normPThis);
        protocolComparison.qCorrelation = safeVectorCorrelation(rawQThis,normQThis);
        protocolComparison.signAgreementFraction = computeFiniteSignAgreement(rawEffectThis,normEffectThis);
        protocolComparison.rawP05Count = sum(rawPThis(:) < cfg.consistencyPThreshold,'omitnan');
        protocolComparison.normalizedP05Count = sum(normPThis(:) < cfg.consistencyPThreshold,'omitnan');
        protocolComparison.rawQ05Count = sum(rawQThis(:) < cfg.fdrAlpha,'omitnan');
        protocolComparison.normalizedQ05Count = sum(normQThis(:) < cfg.fdrAlpha,'omitnan');
        protocolComparison.rawOnlyP05Count = sum(rawPThis(:) < cfg.consistencyPThreshold & ~(normPThis(:) < cfg.consistencyPThreshold),'omitnan');
        protocolComparison.normalizedOnlyP05Count = sum(normPThis(:) < cfg.consistencyPThreshold & ~(rawPThis(:) < cfg.consistencyPThreshold),'omitnan');
        protocolComparison.bothP05Count = sum(rawPThis(:) < cfg.consistencyPThreshold & normPThis(:) < cfg.consistencyPThreshold,'omitnan');
        protocolComparison.effectShiftMean = mean((rawEffectThis(:) - normEffectThis(:)),'omitnan');
        protocolComparison.effectAbsShiftMean = mean(abs(rawEffectThis(:) - normEffectThis(:)),'omitnan');

        if isempty(familySummary.protocolComparisons)
            familySummary.protocolComparisons = repmat(protocolComparison,1,numel(preparedData.protocolNameList));
        end
        familySummary.protocolComparisons(iProtocol) = protocolComparison;
    end

    if isempty(rawVsNormalized.familySummaries)
        rawVsNormalized.familySummaries = repmat(familySummary,1,numel(familyDefs));
    end
    rawVsNormalized.familySummaries(iFamily) = familySummary;
end
end

function pairHeterogeneity = buildPairHeterogeneitySummary(preparedData,cfg)
familyDefs = getFamilyDefinitions();
spaceNames = {'normalized','raw'};

pairHeterogeneity = struct;
pairHeterogeneity.analysisType = 'pair_heterogeneity';
pairHeterogeneity.protocolNameList = preparedData.protocolNameList;
pairHeterogeneity.roiLabels = preparedData.roiLabels;
pairHeterogeneity.pairInfo = preparedData.pairInfo;
pairHeterogeneity.spaceSummaries = [];

for iSpace = 1:numel(spaceNames)
    spaceSummary = struct;
    spaceSummary.spaceName = spaceNames{iSpace};
    spaceSummary.protocolResults = [];

    for iProtocol = 1:numel(preparedData.protocolResults)
        protocolData = preparedData.protocolResults(iProtocol);
        protocolSummary = struct;
        protocolSummary.protocolName = protocolData.protocolName;
        protocolSummary.validPairLabels = protocolData.validPairLabels;
        protocolSummary.familySummaries = [];

        for iFamily = 1:numel(familyDefs)
            pairwiseValues = computePairwiseContributionArray(protocolData.(spaceNames{iSpace}),familyDefs(iFamily).name);
            heterogeneityMetrics = summarizePairwiseContributions(pairwiseValues);

            familySummary = struct;
            familySummary.familyName = familyDefs(iFamily).name;
            familySummary.pairwiseValues = pairwiseValues;
            familySummary.metrics = heterogeneityMetrics;
            familySummary.topDistributedCells = rankPairwiseCells(preparedData,heterogeneityMetrics,cfg);
            if isempty(protocolSummary.familySummaries)
                protocolSummary.familySummaries = repmat(familySummary,1,numel(familyDefs));
            end
            protocolSummary.familySummaries(iFamily) = familySummary;
        end

        if isempty(spaceSummary.protocolResults)
            spaceSummary.protocolResults = repmat(protocolSummary,1,numel(preparedData.protocolResults));
        end
        spaceSummary.protocolResults(iProtocol) = protocolSummary;
    end

    if isempty(pairHeterogeneity.spaceSummaries)
        pairHeterogeneity.spaceSummaries = repmat(spaceSummary,1,numel(spaceNames));
    end
    pairHeterogeneity.spaceSummaries(iSpace) = spaceSummary;
end
end

function globalNetFlow = buildGlobalNetFlowSummary(preparedData,cfg)
globalNetFlow = struct;
globalNetFlow.analysisType = 'global_net_flow';
globalNetFlow.protocolNameList = preparedData.protocolNameList;
globalNetFlow.roiLabels = preparedData.roiLabels;
globalNetFlow.protocolResults = [];

frontalIndices = find(contains(preparedData.roiLabels,'Frontal','IgnoreCase',true));
nonFrontalIndices = setdiff(1:numel(preparedData.roiLabels),frontalIndices);

for iProtocol = 1:numel(preparedData.protocolResults)
    protocolData = preparedData.protocolResults(iProtocol);
    [netFlowData,collapsedMatrices] = computeNetFlowData(protocolData);
    netFlowInference = computeNetFlowInference(netFlowData,cfg,frontalIndices,nonFrontalIndices);

    protocolResult = struct;
    protocolResult.protocolName = protocolData.protocolName;
    protocolResult.validPairLabels = protocolData.validPairLabels;
    protocolResult.collapsedMatrices = collapsedMatrices;
    protocolResult.netFlowData = netFlowData;
    protocolResult.inference = netFlowInference;
    if isempty(globalNetFlow.protocolResults)
        globalNetFlow.protocolResults = repmat(protocolResult,1,numel(preparedData.protocolResults));
    end
    globalNetFlow.protocolResults(iProtocol) = protocolResult;
end
end

function [netFlowData,collapsedMatrices] = computeNetFlowData(protocolData)
numPairs = protocolData.numValidPairs;
numROIs = numel(protocolData.roiLabels);

collapsedMatrices = struct;
collapsedMatrices.meditatorPre = nan(numROIs,numROIs,numPairs);
collapsedMatrices.meditatorPost = nan(numROIs,numROIs,numPairs);
collapsedMatrices.controlPre = nan(numROIs,numROIs,numPairs);
collapsedMatrices.controlPost = nan(numROIs,numROIs,numPairs);
collapsedMatrices.meditatorDelta = nan(numROIs,numROIs,numPairs);
collapsedMatrices.controlDelta = nan(numROIs,numROIs,numPairs);

netFlowData = struct;
netFlowData.meditatorPre = nan(numROIs,numPairs);
netFlowData.meditatorPost = nan(numROIs,numPairs);
netFlowData.controlPre = nan(numROIs,numPairs);
netFlowData.controlPost = nan(numROIs,numPairs);
netFlowData.meditatorDelta = nan(numROIs,numPairs);
netFlowData.controlDelta = nan(numROIs,numPairs);

for iPair = 1:numPairs
    collapsedMatrices.meditatorPre(:,:,iPair) = collapseBandPairGrid(protocolData.raw.preMed(:,:,:,:,iPair));
    collapsedMatrices.meditatorPost(:,:,iPair) = collapseBandPairGrid(protocolData.raw.postMed(:,:,:,:,iPair));
    collapsedMatrices.controlPre(:,:,iPair) = collapseBandPairGrid(protocolData.raw.preCtrl(:,:,:,:,iPair));
    collapsedMatrices.controlPost(:,:,iPair) = collapseBandPairGrid(protocolData.raw.postCtrl(:,:,:,:,iPair));
    collapsedMatrices.meditatorDelta(:,:,iPair) = collapseBandPairGrid(protocolData.raw.deltaMed(:,:,:,:,iPair));
    collapsedMatrices.controlDelta(:,:,iPair) = collapseBandPairGrid(protocolData.raw.deltaCtrl(:,:,:,:,iPair));

    netFlowData.meditatorPre(:,iPair) = cfGCUtils.computeNetFlow(collapsedMatrices.meditatorPre(:,:,iPair));
    netFlowData.meditatorPost(:,iPair) = cfGCUtils.computeNetFlow(collapsedMatrices.meditatorPost(:,:,iPair));
    netFlowData.controlPre(:,iPair) = cfGCUtils.computeNetFlow(collapsedMatrices.controlPre(:,:,iPair));
    netFlowData.controlPost(:,iPair) = cfGCUtils.computeNetFlow(collapsedMatrices.controlPost(:,:,iPair));
    netFlowData.meditatorDelta(:,iPair) = cfGCUtils.computeNetFlow(collapsedMatrices.meditatorDelta(:,:,iPair));
    netFlowData.controlDelta(:,iPair) = cfGCUtils.computeNetFlow(collapsedMatrices.controlDelta(:,:,iPair));
end
end

function collapsed = collapseBandPairGrid(gridData)
collapsed = nan(size(gridData,1),size(gridData,2));
for iRow = 1:size(gridData,1)
    for iCol = 1:size(gridData,2)
        tmp = squeeze(gridData(iRow,iCol,:,:));
        collapsed(iRow,iCol) = mean(tmp(:),'omitnan');
    end
end
end

function inference = computeNetFlowInference(netFlowData,cfg,frontalIndices,nonFrontalIndices)
numROIs = size(netFlowData.meditatorPre,1);

inference = struct;
fieldNames = {'meditatorChangeP','meditatorChangeDz','meditatorChangeN', ...
    'controlChangeP','controlChangeDz','controlChangeN', ...
    'baselinePreGroupP','baselinePreGroupEffect','baselinePreGroupN', ...
    'postGroupP','postGroupEffect','postGroupN', ...
    'interactionP','interactionEffect','interactionN'};
for iField = 1:numel(fieldNames)
    inference.(fieldNames{iField}) = nan(numROIs,1);
end

for iROI = 1:numROIs
    [pVal,effectVal,nUsed] = cfGCUtils.pairedLabelSwapPrePost(netFlowData.meditatorPre(iROI,:),netFlowData.meditatorPost(iROI,:),cfg.numPermutations);
    inference.meditatorChangeP(iROI) = pVal;
    inference.meditatorChangeDz(iROI) = effectVal;
    inference.meditatorChangeN(iROI) = nUsed;

    [pVal,effectVal,nUsed] = cfGCUtils.pairedLabelSwapPrePost(netFlowData.controlPre(iROI,:),netFlowData.controlPost(iROI,:),cfg.numPermutations);
    inference.controlChangeP(iROI) = pVal;
    inference.controlChangeDz(iROI) = effectVal;
    inference.controlChangeN(iROI) = nUsed;

    [pVal,effectVal,nUsed] = cfGCUtils.pairedLabelSwapMatchedGroups(netFlowData.meditatorPre(iROI,:),netFlowData.controlPre(iROI,:),cfg.numPermutations);
    inference.baselinePreGroupP(iROI) = pVal;
    inference.baselinePreGroupEffect(iROI) = effectVal;
    inference.baselinePreGroupN(iROI) = nUsed;

    [pVal,effectVal,nUsed] = cfGCUtils.pairedLabelSwapMatchedGroups(netFlowData.meditatorPost(iROI,:),netFlowData.controlPost(iROI,:),cfg.numPermutations);
    inference.postGroupP(iROI) = pVal;
    inference.postGroupEffect(iROI) = effectVal;
    inference.postGroupN(iROI) = nUsed;

    [pVal,effectVal,nUsed] = cfGCUtils.pairedLabelSwapMatchedGroups(netFlowData.meditatorDelta(iROI,:),netFlowData.controlDelta(iROI,:),cfg.numPermutations);
    inference.interactionP(iROI) = pVal;
    inference.interactionEffect(iROI) = effectVal;
    inference.interactionN(iROI) = nUsed;
end

inference = applyNetFlowMinimumNMask(inference,cfg.minValidSubjectsPerCell);
familyDefs = getFamilyDefinitions();
for iFamily = 1:numel(familyDefs)
    [inference.(familyDefs(iFamily).qField),inference.(familyDefs(iFamily).sigField)] = ...
        bhFDRVector(inference.(familyDefs(iFamily).pField),cfg.fdrAlpha);
end

frontalDominance = struct;
frontalDominance.meditatorPre = mean(netFlowData.meditatorPre(frontalIndices,:),1,'omitnan') - mean(netFlowData.meditatorPre(nonFrontalIndices,:),1,'omitnan');
frontalDominance.meditatorPost = mean(netFlowData.meditatorPost(frontalIndices,:),1,'omitnan') - mean(netFlowData.meditatorPost(nonFrontalIndices,:),1,'omitnan');
frontalDominance.controlPre = mean(netFlowData.controlPre(frontalIndices,:),1,'omitnan') - mean(netFlowData.controlPre(nonFrontalIndices,:),1,'omitnan');
frontalDominance.controlPost = mean(netFlowData.controlPost(frontalIndices,:),1,'omitnan') - mean(netFlowData.controlPost(nonFrontalIndices,:),1,'omitnan');
frontalDominance.meditatorDelta = mean(netFlowData.meditatorDelta(frontalIndices,:),1,'omitnan') - mean(netFlowData.meditatorDelta(nonFrontalIndices,:),1,'omitnan');
frontalDominance.controlDelta = mean(netFlowData.controlDelta(frontalIndices,:),1,'omitnan') - mean(netFlowData.controlDelta(nonFrontalIndices,:),1,'omitnan');

frontalDominanceTests = struct;
[frontalDominanceTests.meditatorChangeP,frontalDominanceTests.meditatorChangeDz,frontalDominanceTests.meditatorChangeN] = ...
    cfGCUtils.pairedLabelSwapPrePost(frontalDominance.meditatorPre,frontalDominance.meditatorPost,cfg.numPermutations);
[frontalDominanceTests.controlChangeP,frontalDominanceTests.controlChangeDz,frontalDominanceTests.controlChangeN] = ...
    cfGCUtils.pairedLabelSwapPrePost(frontalDominance.controlPre,frontalDominance.controlPost,cfg.numPermutations);
[frontalDominanceTests.baselinePreGroupP,frontalDominanceTests.baselinePreGroupEffect,frontalDominanceTests.baselinePreGroupN] = ...
    cfGCUtils.pairedLabelSwapMatchedGroups(frontalDominance.meditatorPre,frontalDominance.controlPre,cfg.numPermutations);
[frontalDominanceTests.postGroupP,frontalDominanceTests.postGroupEffect,frontalDominanceTests.postGroupN] = ...
    cfGCUtils.pairedLabelSwapMatchedGroups(frontalDominance.meditatorPost,frontalDominance.controlPost,cfg.numPermutations);
[frontalDominanceTests.interactionP,frontalDominanceTests.interactionEffect,frontalDominanceTests.interactionN] = ...
    cfGCUtils.pairedLabelSwapMatchedGroups(frontalDominance.meditatorDelta,frontalDominance.controlDelta,cfg.numPermutations);

inference.frontalDominance = frontalDominance;
inference.frontalDominanceTests = frontalDominanceTests;
end

function inference = applyNetFlowMinimumNMask(inference,minN)
familyDefs = getFamilyDefinitions();
for iFamily = 1:numel(familyDefs)
    invalidMask = inference.(familyDefs(iFamily).nField) < minN;
    inference.(familyDefs(iFamily).pField)(invalidMask) = NaN;
    inference.(familyDefs(iFamily).effectField)(invalidMask) = NaN;
end
end

function pairwiseValues = computePairwiseContributionArray(dataStruct,familyName)
switch familyName
    case 'meditatorChange'
        pairwiseValues = dataStruct.postMed - dataStruct.preMed;
    case 'controlChange'
        pairwiseValues = dataStruct.postCtrl - dataStruct.preCtrl;
    case 'baselinePreGroup'
        pairwiseValues = dataStruct.preMed - dataStruct.preCtrl;
    case 'postGroup'
        pairwiseValues = dataStruct.postMed - dataStruct.postCtrl;
    case 'interaction'
        pairwiseValues = (dataStruct.postMed - dataStruct.preMed) - (dataStruct.postCtrl - dataStruct.preCtrl);
    otherwise
        error('Unknown family name: %s',familyName);
end
end

function metrics = summarizePairwiseContributions(pairwiseValues)
dataSize = size(pairwiseValues);
metrics = struct;
metrics.mean = mean(pairwiseValues,5,'omitnan');
metrics.median = median(pairwiseValues,5,'omitnan');
metrics.std = std(pairwiseValues,0,5,'omitnan');
metrics.iqr = computeIQR(pairwiseValues);
metrics.numPairsUsed = sum(isfinite(pairwiseValues),5);
metrics.signAgreementFraction = computeSignAgreementFraction(sign(pairwiseValues));
metrics.maxPairContributionFraction = computeMaxContributionFraction(pairwiseValues);
metrics.effectivePairCount = computeEffectivePairCount(pairwiseValues);
metrics.leaveOneOutSignAgreementFraction = computeLeaveOneOutSignAgreement(pairwiseValues);
metrics.shape = dataSize;
end

function protocolResults = stackInferenceField(protocolResultsIn,fieldName)
numProtocols = numel(protocolResultsIn);
firstField = protocolResultsIn(1).inference.(fieldName);
protocolResults = nan([size(firstField) numProtocols]);
for iProtocol = 1:numProtocols
    protocolResults(:,:,:,:,iProtocol) = protocolResultsIn(iProtocol).inference.(fieldName);
end
end

function signAgreementFraction = computeSignAgreementFraction(signStack)
positiveCount = sum(signStack > 0,5,'omitnan');
negativeCount = sum(signStack < 0,5,'omitnan');
finiteCount = sum(isfinite(signStack) & signStack ~= 0,5);
signAgreementFraction = nan(size(positiveCount));
validMask = finiteCount > 0;
signAgreementFraction(validMask) = max(positiveCount(validMask),negativeCount(validMask)) ./ finiteCount(validMask);
end

function rankedCells = rankCellsByConsistency(preparedData,pStack,qStack,effectStack,meanEffect,meanAbsEffect,signAgreementFraction,cfg)
score = sum(pStack < cfg.consistencyPThreshold,5,'omitnan') + ...
    0.01 * meanAbsEffect + 0.001 * signAgreementFraction;
rankedCells = makeRankedCellStruct(preparedData,score,pStack,qStack,effectStack,meanEffect,meanAbsEffect,signAgreementFraction,cfg.maxRankedCells);
end

function rankedCells = rankCellsByEffect(preparedData,effectStack,consistencyWeightedEffect,meanEffect,meanAbsEffect,signAgreementFraction,cfg)
dummyP = nan(size(effectStack));
dummyQ = nan(size(effectStack));
rankedCells = makeRankedCellStruct(preparedData,consistencyWeightedEffect,dummyP,dummyQ,effectStack,meanEffect,meanAbsEffect,signAgreementFraction,cfg.maxRankedCells);
end

function rankedCells = rankPairwiseCells(preparedData,metrics,cfg)
score = metrics.signAgreementFraction .* metrics.effectivePairCount;
dummyStack = nan([size(score) 1]);
rankedCells = makeRankedCellStruct(preparedData,score,dummyStack,dummyStack,dummyStack,metrics.mean,abs(metrics.mean),metrics.signAgreementFraction,cfg.maxRankedCells);
end

function rankedCells = makeRankedCellStruct(preparedData,score,pStack,qStack,effectStack,meanEffect,meanAbsEffect,signAgreementFraction,maxRankedCells)
scoreVector = score(:);
finiteMask = isfinite(scoreVector);
validIndices = find(finiteMask);
if isempty(validIndices)
    rankedCells = struct([]);
    return;
end

[~,order] = sort(scoreVector(validIndices),'descend');
topIndices = validIndices(order(1:min(maxRankedCells,numel(order))));
rankedCellBuffer = cell(numel(topIndices),1);

for iCell = 1:numel(topIndices)
    linearIndex = topIndices(iCell);
    [iROI,jROI,iPhase,iAmp] = ind2sub(size(score),linearIndex);
    protocolsP = {};
    protocolsQ = {};

    if ndims(pStack) == 5 && size(pStack,5) == numel(preparedData.protocolNameList)
        thisP = squeeze(pStack(iROI,jROI,iPhase,iAmp,:));
        thisQ = squeeze(qStack(iROI,jROI,iPhase,iAmp,:));
        thisEffect = squeeze(effectStack(iROI,jROI,iPhase,iAmp,:));
        protocolsP = preparedData.protocolNameList(thisP < 0.05);
        protocolsQ = preparedData.protocolNameList(thisQ < 0.05);
    else
        thisEffect = [];
    end

    entry = struct;
    entry.rank = iCell;
    entry.sourceROI = preparedData.roiLabels{iROI};
    entry.targetROI = preparedData.roiLabels{jROI};
    entry.phaseBandIndex = preparedData.phaseBandList(iPhase);
    entry.ampBandIndex = preparedData.ampBandList(iAmp);
    entry.phaseRange = preparedData.pairInfo(iPhase,iAmp).phaseRange;
    entry.ampRange = preparedData.pairInfo(iPhase,iAmp).ampRange;
    entry.bandLabel = preparedData.pairInfo(iPhase,iAmp).label;
    entry.score = score(iROI,jROI,iPhase,iAmp);
    entry.meanEffect = meanEffect(iROI,jROI,iPhase,iAmp);
    entry.meanAbsEffect = meanAbsEffect(iROI,jROI,iPhase,iAmp);
    entry.signAgreementFraction = signAgreementFraction(iROI,jROI,iPhase,iAmp);
    entry.protocolCountP05 = numel(protocolsP);
    entry.protocolCountQ05 = numel(protocolsQ);
    entry.protocolsP05 = protocolsP;
    entry.protocolsQ05 = protocolsQ;
    entry.protocolEffects = thisEffect;
    rankedCellBuffer{iCell,1} = entry;
end

if isempty(rankedCellBuffer)
    rankedCells = struct([]);
else
    rankedCells = vertcat(rankedCellBuffer{:});
end
end

function rho = safeVectorCorrelation(a,b)
a = a(:);
b = b(:);
validMask = isfinite(a) & isfinite(b);
if sum(validMask) < 3
    rho = NaN;
    return;
end
rho = corr(a(validMask),b(validMask));
end

function agreementFraction = computeFiniteSignAgreement(a,b)
a = sign(a(:));
b = sign(b(:));
validMask = isfinite(a) & isfinite(b) & a ~= 0 & b ~= 0;
if ~any(validMask)
    agreementFraction = NaN;
    return;
end
agreementFraction = mean(a(validMask) == b(validMask));
end

function qValues = localBHVector(pVector)
qValues = nan(size(pVector));
validMask = isfinite(pVector);
validP = pVector(validMask);
if isempty(validP)
    return;
end

[sortedP,order] = sort(validP(:));
n = numel(sortedP);
qSorted = nan(n,1);
qSorted(n) = sortedP(n);
for i = n-1:-1:1
    qSorted(i) = min((n/i) * sortedP(i),qSorted(i+1));
end
qSorted = min(qSorted,1);

tmp = nan(n,1);
tmp(order) = qSorted;
qValues(validMask) = tmp;
end

function [qValues,sigMask] = bhFDRVector(pValues,alpha)
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.05;
end
qValues = localBHVector(pValues(:));
qValues = reshape(qValues,size(pValues));
sigMask = qValues < alpha;
end

function iqrValues = computeIQR(values)
iqrValues = nan(size(values,1),size(values,2),size(values,3),size(values,4));
for iRow = 1:size(values,1)
    for iCol = 1:size(values,2)
        for iPhase = 1:size(values,3)
            for iAmp = 1:size(values,4)
                tmp = squeeze(values(iRow,iCol,iPhase,iAmp,:));
                tmp = tmp(isfinite(tmp));
                if ~isempty(tmp)
                    iqrValues(iRow,iCol,iPhase,iAmp) = iqr(tmp);
                end
            end
        end
    end
end
end

function fraction = computeMaxContributionFraction(values)
fraction = nan(size(values,1),size(values,2),size(values,3),size(values,4));
for iRow = 1:size(values,1)
    for iCol = 1:size(values,2)
        for iPhase = 1:size(values,3)
            for iAmp = 1:size(values,4)
                tmp = abs(squeeze(values(iRow,iCol,iPhase,iAmp,:)));
                tmp = tmp(isfinite(tmp));
                if isempty(tmp) || sum(tmp) <= eps
                    continue;
                end
                fraction(iRow,iCol,iPhase,iAmp) = max(tmp) / sum(tmp);
            end
        end
    end
end
end

function effectiveN = computeEffectivePairCount(values)
effectiveN = nan(size(values,1),size(values,2),size(values,3),size(values,4));
for iRow = 1:size(values,1)
    for iCol = 1:size(values,2)
        for iPhase = 1:size(values,3)
            for iAmp = 1:size(values,4)
                tmp = abs(squeeze(values(iRow,iCol,iPhase,iAmp,:)));
                tmp = tmp(isfinite(tmp));
                if isempty(tmp) || sum(tmp.^2) <= eps
                    continue;
                end
                effectiveN(iRow,iCol,iPhase,iAmp) = (sum(tmp).^2) / sum(tmp.^2);
            end
        end
    end
end
end

function agreement = computeLeaveOneOutSignAgreement(values)
agreement = nan(size(values,1),size(values,2),size(values,3),size(values,4));
for iRow = 1:size(values,1)
    for iCol = 1:size(values,2)
        for iPhase = 1:size(values,3)
            for iAmp = 1:size(values,4)
                tmp = squeeze(values(iRow,iCol,iPhase,iAmp,:));
                tmp = tmp(isfinite(tmp));
                if numel(tmp) < 3
                    continue;
                end
                fullMean = mean(tmp);
                if fullMean == 0
                    continue;
                end
                looSigns = false(numel(tmp),1);
                for iLeave = 1:numel(tmp)
                    looMean = mean(tmp(setdiff(1:numel(tmp),iLeave)));
                    looSigns(iLeave) = sign(looMean) == sign(fullMean);
                end
                agreement(iRow,iCol,iPhase,iAmp) = mean(looSigns);
            end
        end
    end
end
end

function familyDefs = getFamilyDefinitions()
familyDefs = [...
    struct('name','meditatorChange','testKind','paired_prepost','xField','preMed','yField','postMed','pField','meditatorChangeP','effectField','meditatorChangeDz','nField','meditatorChangeN','qField','meditatorChangeQ','sigField','meditatorChangeSig'), ...
    struct('name','controlChange','testKind','paired_prepost','xField','preCtrl','yField','postCtrl','pField','controlChangeP','effectField','controlChangeDz','nField','controlChangeN','qField','controlChangeQ','sigField','controlChangeSig'), ...
    struct('name','baselinePreGroup','testKind','matched_groups','xField','preMed','yField','preCtrl','pField','baselinePreGroupP','effectField','baselinePreGroupEffect','nField','baselinePreGroupN','qField','baselinePreGroupQ','sigField','baselinePreGroupSig'), ...
    struct('name','postGroup','testKind','matched_groups','xField','postMed','yField','postCtrl','pField','postGroupP','effectField','postGroupEffect','nField','postGroupN','qField','postGroupQ','sigField','postGroupSig'), ...
    struct('name','interaction','testKind','matched_groups','xField','deltaMed','yField','deltaCtrl','pField','interactionP','effectField','interactionEffect','nField','interactionN','qField','interactionQ','sigField','interactionSig')];
end

function s = setDefault(s,fieldName,defaultValue)
if ~isfield(s,fieldName) || isempty(s.(fieldName))
    s.(fieldName) = defaultValue;
end
end
