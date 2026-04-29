function normalizedGrid = normalizeSubjectGrid(gridData)
% normalizeSubjectGrid Z-score a subject grid using its finite entries only.
normalizedGrid = nan(size(gridData));
validValues = gridData(isfinite(gridData));

if numel(validValues) < 2
    return;
end

mu = mean(validValues);
sigma = std(validValues,0);
if sigma <= eps
    return;
end

normalizedGrid = (gridData - mu) / sigma;
end
