function plotHeatmap(ax,matrix,roiLabels,cLims,cMap)
% plotHeatmap Draw a labeled ROI heatmap with consistent styling.
imagesc(ax,matrix,'AlphaData',isfinite(matrix));
set(ax,'YDir','normal');
axis(ax,'square');
colormap(ax,cfGCUtils.resolveColormap(cMap));
colorbar(ax);
set(ax,'XTick',1:numel(roiLabels),'XTickLabel',roiLabels, ...
    'YTick',1:numel(roiLabels),'YTickLabel',roiLabels);
xtickangle(ax,45);
if ~isempty(cLims)
    caxis(ax,cLims);
end
end
