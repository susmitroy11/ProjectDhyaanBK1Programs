function result = runSinglePairComparisonProtocol(pairSelection,protocolName,badEyeCondition,badTrialVersion,cfg)
% runSinglePairComparisonProtocol Build descriptive outputs for one matched pair.
medPreResults = cfGCUtils.loadConditionResults(pairSelection.meditator,protocolName,'pre',badEyeCondition,badTrialVersion,cfg.dataFolder);
medPostResults = cfGCUtils.loadConditionResults(pairSelection.meditator,protocolName,'post',badEyeCondition,badTrialVersion,cfg.dataFolder);
ctrlPreResults = cfGCUtils.loadConditionResults(pairSelection.control,protocolName,'pre',badEyeCondition,badTrialVersion,cfg.dataFolder);
ctrlPostResults = cfGCUtils.loadConditionResults(pairSelection.control,protocolName,'post',badEyeCondition,badTrialVersion,cfg.dataFolder);

cfGCUtils.validateResultCompatibility(medPreResults,medPostResults);
cfGCUtils.validateResultCompatibility(ctrlPreResults,ctrlPostResults);
cfGCUtils.validateResultCompatibility(medPreResults,ctrlPreResults);
cfGCUtils.validateResultCompatibility(medPostResults,ctrlPostResults);

[medPreGrid,medPreSample,pairInfo,phaseBandList,ampBandList] = cfGCUtils.unpackBandPairGrid(medPreResults);
[medPostGrid,medPostSample] = cfGCUtils.unpackBandPairGrid(medPostResults);
[ctrlPreGrid,ctrlPreSample] = cfGCUtils.unpackBandPairGrid(ctrlPreResults);
[ctrlPostGrid,ctrlPostSample] = cfGCUtils.unpackBandPairGrid(ctrlPostResults);

lowSampleMask = medPreSample < cfg.minObservationCount | medPostSample < cfg.minObservationCount | ...
    ctrlPreSample < cfg.minObservationCount | ctrlPostSample < cfg.minObservationCount;
medPreGrid(lowSampleMask) = NaN;
medPostGrid(lowSampleMask) = NaN;
ctrlPreGrid(lowSampleMask) = NaN;
ctrlPostGrid(lowSampleMask) = NaN;

medDelta = medPostGrid - medPreGrid;
ctrlDelta = ctrlPostGrid - ctrlPreGrid;
baselinePreDifference = medPreGrid - ctrlPreGrid;
postDifference = medPostGrid - ctrlPostGrid;
deltaDifference = medDelta - ctrlDelta;
minSampleGrid = min(cat(5,medPreSample,medPostSample,ctrlPreSample,ctrlPostSample),[],5);

overview = struct;
overview.collapsedMeditatorDelta = cfGCUtils.meanAcrossPairs(medDelta);
overview.collapsedControlDelta = cfGCUtils.meanAcrossPairs(ctrlDelta);
overview.collapsedBaselinePreDifference = cfGCUtils.meanAcrossPairs(baselinePreDifference);
overview.collapsedPostDifference = cfGCUtils.meanAcrossPairs(postDifference);
overview.collapsedDeltaDifference = cfGCUtils.meanAcrossPairs(deltaDifference);
overview.collapsedMinSampleCount = cfGCUtils.meanAcrossPairs(minSampleGrid);
overview.meditatorNetDelta = cfGCUtils.computeNetFlow(overview.collapsedMeditatorDelta);
overview.controlNetDelta = cfGCUtils.computeNetFlow(overview.collapsedControlDelta);
overview.deltaDifferenceNetFlow = cfGCUtils.computeNetFlow(overview.collapsedDeltaDifference);

commonCLim = cfGCUtils.resolveCommonCLim({medPreGrid,medPostGrid,ctrlPreGrid,ctrlPostGrid},cfg.commonColorLimits);
diffCLim = cfGCUtils.resolveDiffCLim({baselinePreDifference,postDifference,medDelta,ctrlDelta,deltaDifference},cfg.diffColorLimits);
figureFiles = struct;
outputManifest = {};

if cfg.saveOutputFlag
    prefix = fullfile(cfg.figureFolder,[protocolName '_' pairSelection.pairLabel '_' badEyeCondition '_' badTrialVersion '_']);
    figureFiles.meditatorPre = [prefix 'meditatorPre.png'];
    figureFiles.meditatorPost = [prefix 'meditatorPost.png'];
    figureFiles.controlPre = [prefix 'controlPre.png'];
    figureFiles.controlPost = [prefix 'controlPost.png'];
    figureFiles.preDifference = [prefix 'preDifference.png'];
    figureFiles.postDifference = [prefix 'postDifference.png'];
    figureFiles.meditatorDelta = [prefix 'meditatorDelta.png'];
    figureFiles.controlDelta = [prefix 'controlDelta.png'];
    figureFiles.deltaDifference = [prefix 'deltaDifference.png'];

    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(medPreGrid,pairInfo,medPreResults.roiLabels,commonCLim, ...
        sprintf('%s | %s | PRE',pairSelection.meditator,protocolName),'parula'),figureFiles.meditatorPre);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(medPostGrid,pairInfo,medPostResults.roiLabels,commonCLim, ...
        sprintf('%s | %s | POST',pairSelection.meditator,protocolName),'parula'),figureFiles.meditatorPost);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(ctrlPreGrid,pairInfo,ctrlPreResults.roiLabels,commonCLim, ...
        sprintf('%s | %s | PRE',pairSelection.control,protocolName),'parula'),figureFiles.controlPre);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(ctrlPostGrid,pairInfo,ctrlPostResults.roiLabels,commonCLim, ...
        sprintf('%s | %s | POST',pairSelection.control,protocolName),'parula'),figureFiles.controlPost);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(baselinePreDifference,pairInfo,medPreResults.roiLabels,diffCLim, ...
        sprintf('%s | %s | PRE Meditator - Control',pairSelection.pairLabel,protocolName),'cfGCUtils.blueWhiteRed'),figureFiles.preDifference);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(postDifference,pairInfo,medPostResults.roiLabels,diffCLim, ...
        sprintf('%s | %s | POST difference (Meditator - Control)',pairSelection.pairLabel,protocolName),'cfGCUtils.blueWhiteRed'),figureFiles.postDifference);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(medDelta,pairInfo,medPreResults.roiLabels,diffCLim, ...
        sprintf('%s | %s | POST - PRE',pairSelection.meditator,protocolName),'cfGCUtils.blueWhiteRed'),figureFiles.meditatorDelta);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(ctrlDelta,pairInfo,ctrlPreResults.roiLabels,diffCLim, ...
        sprintf('%s | %s | POST - PRE',pairSelection.control,protocolName),'cfGCUtils.blueWhiteRed'),figureFiles.controlDelta);
    cfGCUtils.saveFigureQuietly(createBandPairGridFigure(deltaDifference,pairInfo,medPreResults.roiLabels,diffCLim, ...
        sprintf('%s | %s | Delta Meditator - Control',pairSelection.pairLabel,protocolName),'cfGCUtils.blueWhiteRed'),figureFiles.deltaDifference);

    outputManifest = [outputManifest; struct2cell(figureFiles)]; %#ok<AGROW>
    if cfg.makeOverviewFigure
        figureFiles.overview = [prefix 'overviewOnly.png'];
        cfGCUtils.saveFigureQuietly(createOverviewFigure(protocolName,pairSelection,medPreResults.roiLabels,overview),figureFiles.overview);
        outputManifest{end+1,1} = figureFiles.overview; %#ok<AGROW>
    end
end

result = struct;
result.protocolName = protocolName;
result.badEyeCondition = badEyeCondition;
result.badTrialVersion = badTrialVersion;
result.comparisonMode = 'single_pair';
result.meditatorSubject = pairSelection.meditator;
result.controlSubject = pairSelection.control;
result.pairIndex = pairSelection.pairIndex;
result.pairLabel = pairSelection.pairLabel;
result.roiLabels = medPreResults.roiLabels;
result.pairInfo = pairInfo;
result.phaseBandList = phaseBandList;
result.ampBandList = ampBandList;
result.minObservationCount = cfg.minObservationCount;
result.lowSampleMask = lowSampleMask;
result.maskedCellCount = sum(lowSampleMask(:));
result.meditator.pre = medPreGrid;
result.meditator.post = medPostGrid;
result.meditator.delta = medDelta;
result.control.pre = ctrlPreGrid;
result.control.post = ctrlPostGrid;
result.control.delta = ctrlDelta;
result.preDifference = baselinePreDifference;
result.postDifference = postDifference;
result.deltaDifference = deltaDifference;
result.sampleCount.meditatorPre = medPreSample;
result.sampleCount.meditatorPost = medPostSample;
result.sampleCount.controlPre = ctrlPreSample;
result.sampleCount.controlPost = ctrlPostSample;
result.sampleCount.minAcrossAll = minSampleGrid;
result.overview = overview;
result.figureFiles = figureFiles;
result.outputManifest = outputManifest;
result.notes = ['Single-pair comparison is descriptive only. ', ...
    'No inferential statistics are reported because one matched pair does not support stable population-level inference.'];
end

function hFig = createOverviewFigure(protocolName,pairSelection,roiLabels,overview)
hFig = cfGCUtils.createFigure([protocolName ' ' pairSelection.pairLabel ' overview'],[100 100 1600 950]);
tiled = tiledlayout(hFig,3,3,'TileSpacing','compact','Padding','compact');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsedBaselinePreDifference,roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean PRE difference');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsedPostDifference,roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean POST difference');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsedMeditatorDelta,roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean meditator delta');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsedControlDelta,roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean control delta');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsedDeltaDifference,roiLabels,[],'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean delta difference');

ax = nexttile(tiled);
bar(ax,[overview.meditatorNetDelta overview.controlNetDelta overview.deltaDifferenceNetFlow],'grouped');
set(ax,'XTick',1:numel(roiLabels),'XTickLabel',roiLabels);
xtickangle(ax,45);
grid(ax,'on');
legend(ax,{'Meditator delta','Control delta','Delta difference'},'Location','best');
ylabel(ax,'Overview-only net flow');
title(ax,'Overview only: collapsed net flow');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsedMinSampleCount,roiLabels,[],'parula');
title(ax,'Overview only: mean min sample count');

cfGCUtils.addFigureTitle(hFig,sprintf('%s | %s | Single-pair overview',protocolName,pairSelection.pairLabel));
end

function hFig = createBandPairGridFigure(gridData,pairInfo,roiLabels,cLims,figureTitle,cMap)
hFig = cfGCUtils.createFigure(figureTitle,[100 100 1550 980]);
tiled = tiledlayout(hFig,size(gridData,3),size(gridData,4),'TileSpacing','compact','Padding','compact');

for iPhase = 1:size(gridData,3)
    for iAmp = 1:size(gridData,4)
        ax = nexttile(tiled);
        cfGCUtils.plotHeatmap(ax,gridData(:,:,iPhase,iAmp),roiLabels,cLims,cMap);
        xlabel(ax,'Target ROI (Amplitude)');
        ylabel(ax,'Source ROI (Phase)');
        title(ax,pairInfo(iPhase,iAmp).label,'FontWeight','normal');
        text(ax,0.02,0.98,'q: n/a', ...
            'Units','normalized', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','top', ...
            'Color','k', ...
            'FontSize',9, ...
            'BackgroundColor','w', ...
            'Margin',2);
    end
end

cfGCUtils.addFigureTitle(hFig,figureTitle);
end
