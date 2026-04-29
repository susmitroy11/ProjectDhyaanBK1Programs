function hFig = createBandPairGridFigure(gridData,pairInfo,roiLabels,cLims,figureTitle,cMap,sigMask,qGrid)
% createBandPairGridFigure Render one heatmap tile for each phase-amplitude pair.
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
