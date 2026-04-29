function result = runSubjectVisualizationProtocol(subjectName,protocolName,badEyeCondition,badTrialVersion,cfg)
% runSubjectVisualizationProtocol Build figure-ready outputs for one subject.
preResults = cfGCUtils.loadConditionResults(subjectName,protocolName,'pre',badEyeCondition,badTrialVersion,cfg.dataFolder);
postResults = cfGCUtils.loadConditionResults(subjectName,protocolName,'post',badEyeCondition,badTrialVersion,cfg.dataFolder);
cfGCUtils.validateResultCompatibility(preResults,postResults);

[preGrid,preSampleGrid,pairInfo,phaseBandList,ampBandList] = cfGCUtils.unpackBandPairGrid(preResults);
[postGrid,postSampleGrid] = cfGCUtils.unpackBandPairGrid(postResults);

if cfg.sampleCountThreshold > 0
    preGrid(preSampleGrid < cfg.sampleCountThreshold) = NaN;
    postGrid(postSampleGrid < cfg.sampleCountThreshold) = NaN;
end

diffGrid = postGrid - preGrid;
minSampleGrid = min(preSampleGrid,postSampleGrid);

overview.collapsedPre = cfGCUtils.meanAcrossPairs(preGrid);
overview.collapsedPost = cfGCUtils.meanAcrossPairs(postGrid);
overview.collapsedDiff = overview.collapsedPost - overview.collapsedPre;
overview.collapsedMinSampleCount = cfGCUtils.meanAcrossPairs(minSampleGrid);
[strongestSrc,strongestDst,strongestProfile] = getStrongestPairProfile(overview.collapsedDiff,preGrid,postGrid,diffGrid,pairInfo);

figureFiles = struct;
outputManifest = {};
commonCLim = cfGCUtils.resolveCommonCLim({preGrid,postGrid},cfg.commonColorLimits);
diffCLim = cfGCUtils.resolveDiffCLim(diffGrid,cfg.diffColorLimits);

if cfg.saveOutputFlag
    figureFiles.pre = fullfile(cfg.figureFolder,[subjectName '_' protocolName '_' badEyeCondition '_' badTrialVersion '_cfGC_preGrid.png']);
    figureFiles.post = fullfile(cfg.figureFolder,[subjectName '_' protocolName '_' badEyeCondition '_' badTrialVersion '_cfGC_postGrid.png']);
    figureFiles.diff = fullfile(cfg.figureFolder,[subjectName '_' protocolName '_' badEyeCondition '_' badTrialVersion '_cfGC_diffGrid.png']);

    cfGCUtils.saveFigureQuietly(cfGCUtils.createBandPairGridFigure(preGrid,pairInfo,postResults.roiLabels,commonCLim, ...
        sprintf('%s | %s | PRE',subjectName,protocolName),'parula'),figureFiles.pre);
    cfGCUtils.saveFigureQuietly(cfGCUtils.createBandPairGridFigure(postGrid,pairInfo,postResults.roiLabels,commonCLim, ...
        sprintf('%s | %s | POST',subjectName,protocolName),'parula'),figureFiles.post);
    cfGCUtils.saveFigureQuietly(cfGCUtils.createBandPairGridFigure(diffGrid,pairInfo,postResults.roiLabels,diffCLim, ...
        sprintf('%s | %s | POST - PRE',subjectName,protocolName),'cfGCUtils.blueWhiteRed'),figureFiles.diff);
    outputManifest = [outputManifest; struct2cell(figureFiles)]; %#ok<AGROW>

    if cfg.showSampleCountFigure
        figureFiles.sampleCount = fullfile(cfg.figureFolder,[subjectName '_' protocolName '_' badEyeCondition '_' badTrialVersion '_cfGC_sampleCountGrid.png']);
        cfGCUtils.saveFigureQuietly(cfGCUtils.createBandPairGridFigure(minSampleGrid,pairInfo,postResults.roiLabels,[], ...
            sprintf('%s | %s | Min Sample Count (PRE, POST)',subjectName,protocolName),'parula'),figureFiles.sampleCount);
        outputManifest{end+1,1} = figureFiles.sampleCount; %#ok<AGROW>
    end

    if cfg.showOverviewFigure
        figureFiles.overview = fullfile(cfg.figureFolder,[subjectName '_' protocolName '_' badEyeCondition '_' badTrialVersion '_cfGC_overviewOnly.png']);
        cfGCUtils.saveFigureQuietly(createOverviewFigure(subjectName,protocolName,postResults.roiLabels,overview, ...
            strongestSrc,strongestDst,strongestProfile,diffCLim),figureFiles.overview);
        outputManifest{end+1,1} = figureFiles.overview; %#ok<AGROW>
    end
end

result = struct;
result.subjectName = subjectName;
result.protocolName = protocolName;
result.badEyeCondition = badEyeCondition;
result.badTrialVersion = badTrialVersion;
result.roiLabels = postResults.roiLabels;
result.phaseBandList = phaseBandList;
result.ampBandList = ampBandList;
result.pairInfo = pairInfo;
result.preGrid = preGrid;
result.postGrid = postGrid;
result.diffGrid = diffGrid;
result.preSampleGrid = preSampleGrid;
result.postSampleGrid = postSampleGrid;
result.minSampleGrid = minSampleGrid;
result.overview = overview;
result.strongestSourceROI = strongestSrc;
result.strongestTargetROI = strongestDst;
result.strongestPairProfile = strongestProfile;
result.figureFiles = figureFiles;
result.outputManifest = outputManifest;
result.notes = ['Single-subject CF-GC visualization is descriptive. ', ...
    'Primary interpretation should use the 9 individual band-pair heatmaps. ', ...
    'Collapsed matrices are overview-only and should not be treated as a primary biological result.'];
end

function hFig = createOverviewFigure(subjectName,protocolName,roiLabels,overview,strongestSrc,strongestDst,strongestProfile,diffCLim)
hFig = cfGCUtils.createFigure([subjectName ' ' protocolName ' overview'],[100 100 1500 900]);
tiled = tiledlayout(hFig,2,3,'TileSpacing','compact','Padding','compact');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsedPre,roiLabels,[],'parula');
title(ax,'Overview only: mean across band pairs (PRE)');
xlabel(ax,'Target ROI');
ylabel(ax,'Source ROI');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsedPost,roiLabels,[],'parula');
title(ax,'Overview only: mean across band pairs (POST)');
xlabel(ax,'Target ROI');
ylabel(ax,'Source ROI');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsedDiff,roiLabels,diffCLim,'cfGCUtils.blueWhiteRed');
title(ax,'Overview only: mean across band pairs (POST - PRE)');
xlabel(ax,'Target ROI');
ylabel(ax,'Source ROI');

ax = nexttile(tiled);
cfGCUtils.plotHeatmap(ax,overview.collapsedMinSampleCount,roiLabels,[],'parula');
title(ax,'Overview only: mean min sample count');
xlabel(ax,'Target ROI');
ylabel(ax,'Source ROI');

ax = nexttile(tiled,[1 2]);
if ~isempty(strongestProfile.pre)
    plot(ax,strongestProfile.pre,'-o','Color',[0.75 0.2 0.2],'LineWidth',1.5); hold(ax,'on');
    plot(ax,strongestProfile.post,'-o','Color',[0.2 0.35 0.8],'LineWidth',1.5);
    plot(ax,strongestProfile.diff,'-s','Color',[0.15 0.15 0.15],'LineWidth',1.5);
    set(ax,'XTick',1:numel(strongestProfile.pairLabels),'XTickLabel',strongestProfile.pairLabels);
    xtickangle(ax,45);
    xlabel(ax,'Band Pair');
    ylabel(ax,'CF-GC');
    legend(ax,{'PRE','POST','POST - PRE'},'Location','best');
    grid(ax,'on');
    title(ax,sprintf('Strongest overview change: %s -> %s',roiLabels{strongestSrc},roiLabels{strongestDst}),'Interpreter','none');
else
    axis(ax,'off');
    text(0.5,0.5,'No valid ROI pair found.','HorizontalAlignment','center','Parent',ax);
end

cfGCUtils.addFigureTitle(hFig,sprintf('%s | %s | Overview-only supplementary figure',subjectName,protocolName));
end

function [sourceIndex,targetIndex,profile] = getStrongestPairProfile(collapsedDiff,preGrid,postGrid,diffGrid,pairInfo)
absMatrix = abs(collapsedDiff);
absMatrix(~isfinite(absMatrix)) = -Inf;
[~,linearIndex] = max(absMatrix(:));

if isempty(linearIndex) || isinf(absMatrix(linearIndex))
    sourceIndex = NaN;
    targetIndex = NaN;
    profile = struct('pre',[],'post',[],'diff',[],'pairLabels',{{}});
    return;
end

[sourceIndex,targetIndex] = ind2sub(size(absMatrix),linearIndex);
profile.pre = [];
profile.post = [];
profile.diff = [];
profile.pairLabels = {};

for iPhase = 1:size(preGrid,3)
    for iAmp = 1:size(preGrid,4)
        profile.pre(end+1) = preGrid(sourceIndex,targetIndex,iPhase,iAmp); %#ok<AGROW>
        profile.post(end+1) = postGrid(sourceIndex,targetIndex,iPhase,iAmp); %#ok<AGROW>
        profile.diff(end+1) = diffGrid(sourceIndex,targetIndex,iPhase,iAmp); %#ok<AGROW>
        profile.pairLabels{end+1} = sprintf('%g-%g | %g-%g', ... %#ok<AGROW>
            pairInfo(iPhase,iAmp).phaseRange(1),pairInfo(iPhase,iAmp).phaseRange(2), ...
            pairInfo(iPhase,iAmp).ampRange(1),pairInfo(iPhase,iAmp).ampRange(2));
    end
end
end
