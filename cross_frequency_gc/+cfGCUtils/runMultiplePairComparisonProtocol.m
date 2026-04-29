function stats = runMultiplePairComparisonProtocol(pairList,protocolName,badEyeCondition,badTrialVersion,cfg)
% runMultiplePairComparisonProtocol Execute one protocol in the multi-pair pipeline.
numPairsRequested = numel(pairList);
stats = struct;
stats.protocolName = protocolName;
stats.badEyeCondition = badEyeCondition;
stats.badTrialVersion = badTrialVersion;
stats.comparisonMode = 'multiple_pairs';
stats.numPermutations = cfg.numPermutations;
stats.cfg = cfg;
stats.requestedPairList = pairList;
stats.validPairMask = false(numPairsRequested,1);
stats.pairAvailability = repmat({'missing'},numPairsRequested,1);

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
        [~,~,pairGridInfo,phaseBandList,ampBandList] = cfGCUtils.unpackBandPairGrid(medPreResults);
        numROIs = numel(medPreResults.roiLabels);
        numPhaseBands = numel(phaseBandList);
        numAmpBands = numel(ampBandList);

        stats.roiLabels = medPreResults.roiLabels;
        stats.pairInfo = pairGridInfo;
        stats.phaseBandList = phaseBandList;
        stats.ampBandList = ampBandList;
        stats.raw.preMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.raw.postMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.raw.preCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.raw.postCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.raw.deltaMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.raw.deltaCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.normalized.preMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.normalized.postMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.normalized.preCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.normalized.postCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.normalized.deltaMed = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.normalized.deltaCtrl = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
        stats.sampleCount.minAcrossAll = nan(numROIs,numROIs,numPhaseBands,numAmpBands,numPairsRequested);
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

    stats.raw.preMed(:,:,:,:,iPair) = medPreGrid;
    stats.raw.postMed(:,:,:,:,iPair) = medPostGrid;
    stats.raw.preCtrl(:,:,:,:,iPair) = ctrlPreGrid;
    stats.raw.postCtrl(:,:,:,:,iPair) = ctrlPostGrid;
    stats.raw.deltaMed(:,:,:,:,iPair) = medPostGrid - medPreGrid;
    stats.raw.deltaCtrl(:,:,:,:,iPair) = ctrlPostGrid - ctrlPreGrid;
    stats.normalized.preMed(:,:,:,:,iPair) = medPreNorm;
    stats.normalized.postMed(:,:,:,:,iPair) = medPostNorm;
    stats.normalized.preCtrl(:,:,:,:,iPair) = ctrlPreNorm;
    stats.normalized.postCtrl(:,:,:,:,iPair) = ctrlPostNorm;
    stats.normalized.deltaMed(:,:,:,:,iPair) = medPostNorm - medPreNorm;
    stats.normalized.deltaCtrl(:,:,:,:,iPair) = ctrlPostNorm - ctrlPreNorm;
    stats.sampleCount.minAcrossAll(:,:,:,:,iPair) = min(cat(5,medPreSample,medPostSample,ctrlPreSample,ctrlPostSample),[],5);
    stats.validPairMask(iPair) = true;
    stats.pairAvailability{iPair} = 'ok';
end

if isempty(templateResults)
    error('None of the requested matched pairs had complete PRE/POST CF-GC files for protocol %s.',protocolName);
end

validIndices = find(stats.validPairMask);
stats.validPairList = pairList(validIndices);
stats.numValidPairs = numel(validIndices);
stats.numRequestedPairs = numPairsRequested;

stats.groupMeans.raw.meditatorPre = mean(stats.raw.preMed(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.raw.meditatorPost = mean(stats.raw.postMed(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.raw.meditatorDelta = mean(stats.raw.deltaMed(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.raw.controlPre = mean(stats.raw.preCtrl(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.raw.controlPost = mean(stats.raw.postCtrl(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.raw.controlDelta = mean(stats.raw.deltaCtrl(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.raw.baselinePreDifference = stats.groupMeans.raw.meditatorPre - stats.groupMeans.raw.controlPre;
stats.groupMeans.raw.postDifference = stats.groupMeans.raw.meditatorPost - stats.groupMeans.raw.controlPost;
stats.groupMeans.raw.interaction = stats.groupMeans.raw.meditatorDelta - stats.groupMeans.raw.controlDelta;

stats.groupMeans.normalized.meditatorPre = mean(stats.normalized.preMed(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.normalized.meditatorPost = mean(stats.normalized.postMed(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.normalized.meditatorDelta = mean(stats.normalized.deltaMed(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.normalized.controlPre = mean(stats.normalized.preCtrl(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.normalized.controlPost = mean(stats.normalized.postCtrl(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.normalized.controlDelta = mean(stats.normalized.deltaCtrl(:,:,:,:,validIndices),5,'omitnan');
stats.groupMeans.normalized.baselinePreDifference = stats.groupMeans.normalized.meditatorPre - stats.groupMeans.normalized.controlPre;
stats.groupMeans.normalized.postDifference = stats.groupMeans.normalized.meditatorPost - stats.groupMeans.normalized.controlPost;
stats.groupMeans.normalized.interaction = stats.groupMeans.normalized.meditatorDelta - stats.groupMeans.normalized.controlDelta;

stats.inference = computeInference(stats.normalized,validIndices,cfg);
stats.summaryOverview = buildOverview(stats.groupMeans.raw,stats.inference,stats.sampleCount.minAcrossAll(:,:,:,:,validIndices));

figureFiles = struct;
outputManifest = {};
if cfg.saveOutputFlag
    prefix = fullfile(cfg.figureFolder,[protocolName '_' badEyeCondition '_' badTrialVersion '_multiplePairs_']);
    figureFiles.meditatorPre = [prefix 'meditatorPre.png'];
    figureFiles.meditatorPost = [prefix 'meditatorPost.png'];
    figureFiles.controlPre = [prefix 'controlPre.png'];
    figureFiles.controlPost = [prefix 'controlPost.png'];
    figureFiles.meditatorDelta = [prefix 'meditatorDelta.png'];
    figureFiles.controlDelta = [prefix 'controlDelta.png'];
    figureFiles.preDifference = [prefix 'preDifference.png'];
    figureFiles.postDifference = [prefix 'postDifference.png'];
    figureFiles.deltaDifference = [prefix 'deltaDifference.png'];

    commonCLim = cfGCUtils.resolveCommonCLim({stats.groupMeans.raw.meditatorPre,stats.groupMeans.raw.meditatorPost, ...
        stats.groupMeans.raw.controlPre,stats.groupMeans.raw.controlPost},cfg.commonColorLimits);
    rawDiffCLim = cfGCUtils.resolveDiffCLim({stats.groupMeans.raw.meditatorDelta,stats.groupMeans.raw.controlDelta, ...
        stats.groupMeans.raw.baselinePreDifference,stats.groupMeans.raw.postDifference,stats.groupMeans.raw.interaction},cfg.diffColorLimits);

    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(stats.groupMeans.raw.meditatorPre,stats.pairInfo,stats.roiLabels,commonCLim, ...
        sprintf('%s | Selected meditator PRE mean',protocolName),'parula'),figureFiles.meditatorPre);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(stats.groupMeans.raw.meditatorPost,stats.pairInfo,stats.roiLabels,commonCLim, ...
        sprintf('%s | Selected meditator POST mean',protocolName),'parula'),figureFiles.meditatorPost);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(stats.groupMeans.raw.controlPre,stats.pairInfo,stats.roiLabels,commonCLim, ...
        sprintf('%s | Selected control PRE mean',protocolName),'parula'),figureFiles.controlPre);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(stats.groupMeans.raw.controlPost,stats.pairInfo,stats.roiLabels,commonCLim, ...
        sprintf('%s | Selected control POST mean',protocolName),'parula'),figureFiles.controlPost);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(stats.groupMeans.raw.meditatorDelta,stats.pairInfo,stats.roiLabels,rawDiffCLim, ...
        sprintf('%s | Selected meditator POST - PRE mean',protocolName),'cfGCUtils.blueWhiteRed',stats.inference.meditatorChangeSig,stats.inference.meditatorChangeQ),figureFiles.meditatorDelta);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(stats.groupMeans.raw.controlDelta,stats.pairInfo,stats.roiLabels,rawDiffCLim, ...
        sprintf('%s | Selected control POST - PRE mean',protocolName),'cfGCUtils.blueWhiteRed',stats.inference.controlChangeSig,stats.inference.controlChangeQ),figureFiles.controlDelta);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(stats.groupMeans.raw.baselinePreDifference,stats.pairInfo,stats.roiLabels,rawDiffCLim, ...
        sprintf('%s | PRE Meditator - Control',protocolName),'cfGCUtils.blueWhiteRed',stats.inference.baselinePreGroupSig,stats.inference.baselinePreGroupQ),figureFiles.preDifference);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(stats.groupMeans.raw.postDifference,stats.pairInfo,stats.roiLabels,rawDiffCLim, ...
        sprintf('%s | POST Meditator - Control',protocolName),'cfGCUtils.blueWhiteRed',stats.inference.postGroupSig,stats.inference.postGroupQ),figureFiles.postDifference);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(stats.groupMeans.raw.interaction,stats.pairInfo,stats.roiLabels,rawDiffCLim, ...
        sprintf('%s | Delta Meditator - Control',protocolName),'cfGCUtils.blueWhiteRed',stats.inference.interactionSig,stats.inference.interactionQ),figureFiles.deltaDifference);

    outputManifest = [outputManifest; struct2cell(figureFiles)]; %#ok<AGROW>
    if cfg.makeOverviewFigure
        figureFiles.overview = [prefix 'overviewOnly.png'];
        cfGCUtils.saveFigureQuietly(createOverviewFigure(stats),figureFiles.overview);
        outputManifest{end+1,1} = figureFiles.overview; %#ok<AGROW>
    end
end

stats.figureFiles = figureFiles;
stats.outputManifest = outputManifest;
stats.notes = ['Multiple-pair comparison uses matched meditator-control pairs only. ', ...
    'Inference is based on subject-normalized matrices, within-subject pre/post label swaps, ', ...
    'matched-group label swaps, and FDR applied separately within each band pair.'];
end

function inference = computeInference(dataStruct,validIndices,cfg)
dataSize = size(dataStruct.preMed);
numROIs = dataSize(1);
numPhase = dataSize(3);
numAmp = dataSize(4);

inference = struct;
fieldList = {'meditatorChangeP','meditatorChangeDz','meditatorChangeN', ...
    'controlChangeP','controlChangeDz','controlChangeN', ...
    'baselinePreGroupP','baselinePreGroupEffect','baselinePreGroupN', ...
    'postGroupP','postGroupEffect','postGroupN', ...
    'interactionP','interactionEffect','interactionN'};

for iField = 1:numel(fieldList)
    inference.(fieldList{iField}) = nan(numROIs,numROIs,numPhase,numAmp);
end

for iROI = 1:numROIs
    for jROI = 1:numROIs
        for iPhase = 1:numPhase
            for iAmp = 1:numAmp
                preMed = squeeze(dataStruct.preMed(iROI,jROI,iPhase,iAmp,validIndices));
                postMed = squeeze(dataStruct.postMed(iROI,jROI,iPhase,iAmp,validIndices));
                preCtrl = squeeze(dataStruct.preCtrl(iROI,jROI,iPhase,iAmp,validIndices));
                postCtrl = squeeze(dataStruct.postCtrl(iROI,jROI,iPhase,iAmp,validIndices));

                [pVal,effectDz,nUsed] = cfGCUtils.pairedLabelSwapPrePost(preMed,postMed,cfg.numPermutations);
                inference.meditatorChangeP(iROI,jROI,iPhase,iAmp) = pVal;
                inference.meditatorChangeDz(iROI,jROI,iPhase,iAmp) = effectDz;
                inference.meditatorChangeN(iROI,jROI,iPhase,iAmp) = nUsed;

                [pVal,effectDz,nUsed] = cfGCUtils.pairedLabelSwapPrePost(preCtrl,postCtrl,cfg.numPermutations);
                inference.controlChangeP(iROI,jROI,iPhase,iAmp) = pVal;
                inference.controlChangeDz(iROI,jROI,iPhase,iAmp) = effectDz;
                inference.controlChangeN(iROI,jROI,iPhase,iAmp) = nUsed;

                [pVal,effectVal,nUsed] = cfGCUtils.pairedLabelSwapMatchedGroups(preMed,preCtrl,cfg.numPermutations);
                inference.baselinePreGroupP(iROI,jROI,iPhase,iAmp) = pVal;
                inference.baselinePreGroupEffect(iROI,jROI,iPhase,iAmp) = effectVal;
                inference.baselinePreGroupN(iROI,jROI,iPhase,iAmp) = nUsed;

                [pVal,effectVal,nUsed] = cfGCUtils.pairedLabelSwapMatchedGroups(postMed,postCtrl,cfg.numPermutations);
                inference.postGroupP(iROI,jROI,iPhase,iAmp) = pVal;
                inference.postGroupEffect(iROI,jROI,iPhase,iAmp) = effectVal;
                inference.postGroupN(iROI,jROI,iPhase,iAmp) = nUsed;

                deltaMed = postMed - preMed;
                deltaCtrl = postCtrl - preCtrl;
                [pVal,effectVal,nUsed] = cfGCUtils.pairedLabelSwapMatchedGroups(deltaMed,deltaCtrl,cfg.numPermutations);
                inference.interactionP(iROI,jROI,iPhase,iAmp) = pVal;
                inference.interactionEffect(iROI,jROI,iPhase,iAmp) = effectVal;
                inference.interactionN(iROI,jROI,iPhase,iAmp) = nUsed;
            end
        end
    end
end

inference = applyMinimumNMask(inference,cfg.minValidSubjectsPerCell);
[inference.meditatorChangeQ,inference.meditatorChangeSig] = cfGCUtils.bhFDRByPair(inference.meditatorChangeP,cfg.fdrAlpha);
[inference.controlChangeQ,inference.controlChangeSig] = cfGCUtils.bhFDRByPair(inference.controlChangeP,cfg.fdrAlpha);
[inference.baselinePreGroupQ,inference.baselinePreGroupSig] = cfGCUtils.bhFDRByPair(inference.baselinePreGroupP,cfg.fdrAlpha);
[inference.postGroupQ,inference.postGroupSig] = cfGCUtils.bhFDRByPair(inference.postGroupP,cfg.fdrAlpha);
[inference.interactionQ,inference.interactionSig] = cfGCUtils.bhFDRByPair(inference.interactionP,cfg.fdrAlpha);
end

function inference = applyMinimumNMask(inference,minN)
maskFields = {
    'meditatorChangeP','meditatorChangeDz','meditatorChangeN';
    'controlChangeP','controlChangeDz','controlChangeN';
    'baselinePreGroupP','baselinePreGroupEffect','baselinePreGroupN';
    'postGroupP','postGroupEffect','postGroupN';
    'interactionP','interactionEffect','interactionN'};

for iRow = 1:size(maskFields,1)
    nField = maskFields{iRow,3};
    invalidMask = inference.(nField) < minN;
    inference.(maskFields{iRow,1})(invalidMask) = NaN;
    inference.(maskFields{iRow,2})(invalidMask) = NaN;
end
end

function summaryOverview = buildOverview(groupMeans,inference,minSampleData)
summaryOverview.collapsedMeditatorDelta = cfGCUtils.meanAcrossPairs(groupMeans.meditatorDelta);
summaryOverview.collapsedControlDelta = cfGCUtils.meanAcrossPairs(groupMeans.controlDelta);
summaryOverview.collapsedInteraction = cfGCUtils.meanAcrossPairs(groupMeans.interaction);
summaryOverview.collapsedBaselinePreDifference = cfGCUtils.meanAcrossPairs(groupMeans.baselinePreDifference);
summaryOverview.collapsedPostDifference = cfGCUtils.meanAcrossPairs(groupMeans.postDifference);
summaryOverview.collapsedInteractionMinQ = cfGCUtils.minAcrossPairs(inference.interactionQ);
summaryOverview.collapsedBaselineMinQ = cfGCUtils.minAcrossPairs(inference.baselinePreGroupQ);
summaryOverview.collapsedPostMinQ = cfGCUtils.minAcrossPairs(inference.postGroupQ);
summaryOverview.collapsedInteractionSigCount = cfGCUtils.sumAcrossPairs(double(inference.interactionSig));
summaryOverview.collapsedMinSampleCount = cfGCUtils.meanAcrossPairs(mean(minSampleData,5,'omitnan'));
summaryOverview.meditatorNetDelta = cfGCUtils.computeNetFlow(summaryOverview.collapsedMeditatorDelta);
summaryOverview.controlNetDelta = cfGCUtils.computeNetFlow(summaryOverview.collapsedControlDelta);
end

function hFig = createOverviewFigure(stats)
hFig = cfGCUtils.createFigure([stats.protocolName ' overview'],[100 100 1600 950]);
tiled = tiledlayout(hFig,3,3,'TileSpacing','compact','Padding','compact');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,stats.summaryOverview.collapsedMeditatorDelta,stats.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (meditator delta)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,stats.summaryOverview.collapsedControlDelta,stats.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (control delta)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,stats.summaryOverview.collapsedPostMinQ,stats.roiLabels,[0 0.1],'parula');
title(ax,'Overview only: min POST-group q across band pairs');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,stats.summaryOverview.collapsedInteraction,stats.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (interaction)');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,stats.summaryOverview.collapsedInteractionMinQ,stats.roiLabels,[0 0.1],'parula');
title(ax,'Overview only: min interaction q across band pairs');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,stats.summaryOverview.collapsedBaselineMinQ,stats.roiLabels,[0 0.1],'parula');
title(ax,'Overview only: min baseline q across band pairs');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,stats.summaryOverview.collapsedPostDifference,stats.roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (POST difference)');

ax = nexttile(tiled);
bar(ax,[stats.summaryOverview.meditatorNetDelta stats.summaryOverview.controlNetDelta],'grouped');
set(ax,'XTick',1:numel(stats.roiLabels),'XTickLabel',stats.roiLabels);
xtickangle(ax,45);
grid(ax,'on');
legend(ax,{'Meditator','Control'},'Location','best');
ylabel(ax,'Overview-only net delta');
title(ax,'Overview only: collapsed net flow');

cfGCUtils.addFigureTitle(hFig,sprintf('%s | Overview-only supplementary summary',stats.protocolName));
end

function hFig = createBandPairGridFigure(gridData,pairInfo,roiLabels,cLims,figureTitle,cMap,sigMask,qGrid)
if ~exist('sigMask','var')
    sigMask = [];
end
if ~exist('qGrid','var')
    qGrid = [];
end

hFig = cfGCUtils.createFigure(figureTitle,[100 100 1550 980]);
tiled = tiledlayout(hFig,size(gridData,3),size(gridData,4),'TileSpacing','compact','Padding','compact');

for iPhase = 1:size(gridData,3)
    for iAmp = 1:size(gridData,4)
        ax = nexttile(tiled);
        cfGCUtils.plotHeatmap(ax,gridData(:,:,iPhase,iAmp),roiLabels,cLims,cMap);
        xlabel(ax,'Target ROI (Amplitude)');
        ylabel(ax,'Source ROI (Phase)');
        title(ax,pairInfo(iPhase,iAmp).label,'FontWeight','normal');

        if ~isempty(sigMask)
            hold(ax,'on');
            [rowIndex,colIndex] = find(sigMask(:,:,iPhase,iAmp));
            if ~isempty(rowIndex)
                plot(ax,colIndex,rowIndex,'k.','MarkerSize',14);
            end
            hold(ax,'off');
        end

        addQValueLabel(ax,qGrid,iPhase,iAmp);
    end
end

cfGCUtils.addFigureTitle(hFig,figureTitle);
end

function addQValueLabel(ax,qGrid,iPhase,iAmp)
qLabel = 'q: n/a';

if ~isempty(qGrid)
    qMin = cfGCUtils.minFiniteValue(qGrid(:,:,iPhase,iAmp));
    if isfinite(qMin)
        qLabel = sprintf('min q = %.3f',qMin);
    end
end

text(ax,0.02,0.98,qLabel, ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'Color','k', ...
    'FontSize',9, ...
    'BackgroundColor','w', ...
    'Margin',2);
end
