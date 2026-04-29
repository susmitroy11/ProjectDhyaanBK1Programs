function collapsed = meanAcrossPairs(gridData)
% meanAcrossPairs Average across all phase/amplitude band-pair slices.
collapsed = nan(size(gridData,1),size(gridData,2));
for iRow = 1:size(gridData,1)
    for iCol = 1:size(gridData,2)
        tmp = squeeze(gridData(iRow,iCol,:,:));
        collapsed(iRow,iCol) = mean(tmp(:),'omitnan');
    end
end
end
