function collapsed = minAcrossPairs(gridData)
% minAcrossPairs Take the minimum across all phase/amplitude band-pair slices.
collapsed = nan(size(gridData,1),size(gridData,2));
for iRow = 1:size(gridData,1)
    for iCol = 1:size(gridData,2)
        tmp = squeeze(gridData(iRow,iCol,:,:));
        tmp = tmp(isfinite(tmp));
        if ~isempty(tmp)
            collapsed(iRow,iCol) = min(tmp);
        end
    end
end
end
