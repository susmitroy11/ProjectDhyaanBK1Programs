function viz = visualizeGrangerSingleSubjectBK1(subjectName,protocolName,bandSpec,badEyeCondition,badTrialVersion,sourceElectrode,targetElectrode,saveFigureFlag)
% visualizeGrangerSingleSubjectBK1 Visualize one BK1 subject's Granger result.
%
% Example:
%   visualizeGrangerSingleSubjectBK1('028HB','M1','alpha');
%   visualizeGrangerSingleSubjectBK1('028HB','M1',[20 32], 'ep', 'v8', 2, 10, 1);

if ~exist('subjectName','var') || isempty(subjectName);         subjectName = '028HB'; end
if ~exist('protocolName','var') || isempty(protocolName);       protocolName = 'M1';    end
if ~exist('bandSpec','var') || isempty(bandSpec);               bandSpec = 'alpha';     end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition); badEyeCondition = 'ep'; end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion); badTrialVersion = 'v8'; end
if ~exist('sourceElectrode','var');                             sourceElectrode = [];    end
if ~exist('targetElectrode','var');                             targetElectrode = [];    end
if ~exist('saveFigureFlag','var') || isempty(saveFigureFlag);   saveFigureFlag = 0;     end

grangerFolder = fullfile(fileparts(mfilename('fullpath')),'savedDataGranger');
dataFileName = fullfile(grangerFolder,subjectName,[protocolName '_' badEyeCondition '_' badTrialVersion '_granger.mat']);

if ~exist(dataFileName,'file')
    error('Could not find Granger file: %s',dataFileName);
end

x = load(dataFileName);
bandInfo = resolveBandSpec(bandSpec);

goodElectrodes = intersect(getSelectedElectrodes(x.metaPre),getSelectedElectrodes(x.metaPost),'stable');
if isempty(goodElectrodes)
    error('No common good electrodes found for %s %s.',subjectName,protocolName);
end

freqMaskPre = x.freqPre >= bandInfo.range(1) & x.freqPre <= bandInfo.range(2);
freqMaskPost = x.freqPost >= bandInfo.range(1) & x.freqPost <= bandInfo.range(2);
if ~any(freqMaskPre) || ~any(freqMaskPost)
    error('No frequencies found in requested band [%g %g] Hz.',bandInfo.range(1),bandInfo.range(2));
end

Gpre = averageBandMatrix(x.grangerPre,freqMaskPre);
Gpost = averageBandMatrix(x.grangerPost,freqMaskPost);
Gdiff = Gpost - Gpre;

Gpre = removeSelfConnections(Gpre);
Gpost = removeSelfConnections(Gpost);
Gdiff = removeSelfConnections(Gdiff);

[roiGroups,roiLabels] = getDefaultROIGroups();
roiPre = computeROIMatrix(Gpre,roiGroups);
roiPost = computeROIMatrix(Gpost,roiGroups);
roiDiff = computeROIMatrix(Gdiff,roiGroups);

if isempty(sourceElectrode) || isempty(targetElectrode)
    [sourceElectrode,targetElectrode] = findStrongestPair(Gdiff,goodElectrodes);
end

validatePair(sourceElectrode,targetElectrode,goodElectrodes);

pairSpectra.preForward = squeeze(x.grangerPre(sourceElectrode,targetElectrode,:));
pairSpectra.preReverse = squeeze(x.grangerPre(targetElectrode,sourceElectrode,:));
pairSpectra.postForward = squeeze(x.grangerPost(sourceElectrode,targetElectrode,:));
pairSpectra.postReverse = squeeze(x.grangerPost(targetElectrode,sourceElectrode,:));

netFlow = computeNetFlow(Gdiff);
plotElectrodes = chooseTopElectrodes(netFlow,goodElectrodes,15);

hFig = figure('Color','w','Name',[subjectName ' ' protocolName ' Granger']);
hFig.Position(3:4) = [1600 900];

tickPos = unique(round(linspace(1,numel(goodElectrodes),min(10,numel(goodElectrodes)))));
tickLabels = makeTickLabels(x.label(goodElectrodes(tickPos)),goodElectrodes(tickPos));

commonCLim = [0 max([maxFiniteValue(Gpre(goodElectrodes,goodElectrodes)) maxFiniteValue(Gpost(goodElectrodes,goodElectrodes))])];
diffCLim = symmetricCLim(Gdiff(goodElectrodes,goodElectrodes));

subplot(2,3,1);
plotHeatmap(gca,Gpre(goodElectrodes,goodElectrodes),tickPos,tickLabels,commonCLim);
title(sprintf('Pre: %s (%g-%g Hz)',bandInfo.name,bandInfo.range(1),bandInfo.range(2)));

subplot(2,3,2);
plotHeatmap(gca,Gpost(goodElectrodes,goodElectrodes),tickPos,tickLabels,commonCLim);
title(sprintf('Post: %s (%g-%g Hz)',bandInfo.name,bandInfo.range(1),bandInfo.range(2)));

subplot(2,3,3);
plotHeatmap(gca,Gdiff(goodElectrodes,goodElectrodes),tickPos,tickLabels,diffCLim);
title('Post - Pre');

subplot(2,3,4);
plotHeatmap(gca,roiDiff,1:numel(roiLabels),roiLabels,symmetricCLim(roiDiff));
title('ROI Granger Change');

subplot(2,3,5);
plot(x.freqPost,pairSpectra.preForward,'Color',[0.8 0.2 0.2],'LineWidth',1.5); hold on;
plot(x.freqPost,pairSpectra.postForward,'Color',[0.5 0 0],'LineWidth',2);
plot(x.freqPost,pairSpectra.preReverse,'Color',[0.2 0.4 0.8],'LineWidth',1.5);
plot(x.freqPost,pairSpectra.postReverse,'Color',[0 0 0.5],'LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('Granger');
grid on;
legend({sprintf('%s -> %s (pre)',x.label{sourceElectrode},x.label{targetElectrode}), ...
        sprintf('%s -> %s (post)',x.label{sourceElectrode},x.label{targetElectrode}), ...
        sprintf('%s -> %s (pre)',x.label{targetElectrode},x.label{sourceElectrode}), ...
        sprintf('%s -> %s (post)',x.label{targetElectrode},x.label{sourceElectrode})}, ...
        'Interpreter','none','Location','best');
title('Bidirectional Pair Spectrum');

subplot(2,3,6);
bar(netFlow(plotElectrodes),'FaceColor',[0.2 0.6 0.8]); hold on;
line(xlim,[0 0],'Color',[0 0 0],'LineStyle','--');
set(gca,'XTick',1:numel(plotElectrodes),'XTickLabel',x.label(plotElectrodes));
xtickangle(90);
ylabel('Net Flow');
title('Top Net-Flow Electrodes');
grid on;

if exist('sgtitle','file') == 2
    sgtitle(sprintf('%s | %s | %s | Good elecs: %d | Trials: %d',subjectName,protocolName,bandInfo.name,numel(goodElectrodes),x.metaPost.numTrials));
end

viz = struct;
viz.subjectName = subjectName;
viz.protocolName = protocolName;
viz.bandInfo = bandInfo;
viz.goodElectrodes = goodElectrodes;
viz.goodLabels = x.label(goodElectrodes);
viz.sourceElectrode = sourceElectrode;
viz.targetElectrode = targetElectrode;
viz.Gpre = Gpre;
viz.Gpost = Gpost;
viz.Gdiff = Gdiff;
viz.roiPre = roiPre;
viz.roiPost = roiPost;
viz.roiDiff = roiDiff;
viz.netFlow = netFlow;
viz.figureHandle = hFig;

if saveFigureFlag
    figureFolder = fullfile(grangerFolder,'figures','singleSubject');
    ensureFolder(figureFolder);
    saveas(hFig,fullfile(figureFolder,[subjectName '_' protocolName '_' bandInfo.name '_' badEyeCondition '_' badTrialVersion '.png']));
end
end

function selectedElectrodes = getSelectedElectrodes(metaStruct)
selectedElectrodes = [];
if isfield(metaStruct,'selectedElectrodes')
    selectedElectrodes = metaStruct.selectedElectrodes(:)';
end
end

function bandInfo = resolveBandSpec(bandSpec)
defaultBands = getDefaultBands();

if isnumeric(bandSpec) && numel(bandSpec) == 2
    bandInfo.name = sprintf('%g_%gHz',bandSpec(1),bandSpec(2));
    bandInfo.range = bandSpec(:)';
    return;
end

if ~ischar(bandSpec) && ~(isstring(bandSpec) && isscalar(bandSpec))
    error('bandSpec must be a band name or a [fLow fHigh] range.');
end

bandName = char(bandSpec);
for iBand = 1:numel(defaultBands)
    if strcmpi(defaultBands(iBand).name,bandName)
        bandInfo = defaultBands(iBand);
        return;
    end
end

error('Unknown band name: %s',bandName);
end

function defaultBands = getDefaultBands()
defaultBands(1).name = 'alpha';      defaultBands(1).range = [7 10];
defaultBands(2).name = 'beta';       defaultBands(2).range = [20 32];
defaultBands(3).name = 'highgamma';  defaultBands(3).range = [30 80];
end

function G = averageBandMatrix(grangerData,freqMask)
G = mean(grangerData(:,:,freqMask),3,'omitnan');
end

function G = removeSelfConnections(G)
n = size(G,1);
G(1:n+1:end) = NaN;
end

function [roiGroups,roiLabels] = getDefaultROIGroups()
roiGroups{1} = [14:16 32+(12:15)];         roiLabels{1} = 'Occipital_L';
roiGroups{2} = [18:20 32+(17:20)];         roiLabels{2} = 'Occipital_R';
roiGroups{3} = [6:8 11 12 32+[7:9 11]];    roiLabels{3} = 'Central_L';
roiGroups{4} = [22 23 25 28 29 32+[22 24:26]]; roiLabels{4} = 'Central_R';
roiGroups{5} = [1 3 4 32+[1 2 4 5]];       roiLabels{5} = 'Frontal_L';
roiGroups{6} = [30:32 32+(28:31)];         roiLabels{6} = 'Frontal_R';
end

function roiMatrix = computeROIMatrix(channelMatrix,roiGroups)
numROIs = numel(roiGroups);
roiMatrix = nan(numROIs,numROIs);

for iROI = 1:numROIs
    srcIdx = roiGroups{iROI};
    for jROI = 1:numROIs
        dstIdx = roiGroups{jROI};
        block = channelMatrix(srcIdx,dstIdx);
        if iROI == jROI
            block(1:size(block,1)+1:end) = NaN;
        end
        roiMatrix(iROI,jROI) = mean(block(:),'omitnan');
    end
end
end

function [sourceElectrode,targetElectrode] = findStrongestPair(Gdiff,goodElectrodes)
block = abs(Gdiff(goodElectrodes,goodElectrodes));
block(~isfinite(block)) = -Inf;

[~,linearIndex] = max(block(:));
if isinf(block(linearIndex))
    error('Could not find a valid directed pair to plot.');
end

[rowPos,colPos] = ind2sub(size(block),linearIndex);
sourceElectrode = goodElectrodes(rowPos);
targetElectrode = goodElectrodes(colPos);
end

function validatePair(sourceElectrode,targetElectrode,goodElectrodes)
if ~ismember(sourceElectrode,goodElectrodes)
    error('Source electrode %d is not in the selected electrode set.',sourceElectrode);
end
if ~ismember(targetElectrode,goodElectrodes)
    error('Target electrode %d is not in the selected electrode set.',targetElectrode);
end
if sourceElectrode == targetElectrode
    error('Source and target electrodes must be different.');
end
end

function netFlow = computeNetFlow(channelMatrix)
outFlow = mean(channelMatrix,2,'omitnan');
inFlow = mean(channelMatrix,1,'omitnan')';
netFlow = outFlow - inFlow;
end

function plotElectrodes = chooseTopElectrodes(netFlow,goodElectrodes,maxCount)
[~,order] = sort(abs(netFlow(goodElectrodes)),'descend');
numToPlot = min(maxCount,numel(order));
plotElectrodes = goodElectrodes(order(1:numToPlot));
end

function plotHeatmap(ax,matrix,tickPos,tickLabels,cLims)
imagesc(ax,matrix,'AlphaData',isfinite(matrix));
set(ax,'YDir','normal');
axis(ax,'square');
colorbar(ax);
colormap(ax,parula);
set(ax,'XTick',tickPos,'XTickLabel',tickLabels,'YTick',tickPos,'YTickLabel',tickLabels);
xtickangle(ax,90);
if ~isempty(cLims) && all(isfinite(cLims))
    caxis(ax,cLims);
end
end

function tickLabels = makeTickLabels(labelList,electrodeNumbers)
tickLabels = cell(size(labelList));
for iLabel = 1:numel(labelList)
    tickLabels{iLabel} = sprintf('%s(%d)',labelList{iLabel},electrodeNumbers(iLabel));
end
end

function cLim = symmetricCLim(matrix)
maxAbs = max(abs(matrix(isfinite(matrix))));
if isempty(maxAbs) || ~isfinite(maxAbs)
    cLim = [];
else
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

function ensureFolder(folderName)
if ~exist(folderName,'dir')
    mkdir(folderName);
end
end
