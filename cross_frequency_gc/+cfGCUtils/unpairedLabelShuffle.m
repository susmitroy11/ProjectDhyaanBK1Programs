function [pValue,effectD,nUsed] = unpairedLabelShuffle(xValues,yValues,numPermutations)
% unpairedLabelShuffle Permutation test for independent-group differences.
xValues = xValues(isfinite(xValues));
yValues = yValues(isfinite(yValues));
xValues = xValues(:);
yValues = yValues(:);

nX = numel(xValues);
nY = numel(yValues);
nUsed = min(nX,nY);

if nX < 2 || nY < 2
    pValue = NaN;
    effectD = NaN;
    return;
end

observed = mean(xValues) - mean(yValues);
combined = [xValues; yValues];
permStats = nan(numPermutations,1);

for iPerm = 1:numPermutations
    order = randperm(numel(combined));
    xPerm = combined(order(1:nX));
    yPerm = combined(order(nX+1:end));
    permStats(iPerm) = mean(xPerm) - mean(yPerm);
end

pValue = (sum(abs(permStats) >= abs(observed)) + 1) / (numPermutations + 1);
pooledStd = sqrt(((nX - 1) * var(xValues,0) + (nY - 1) * var(yValues,0)) / max(nX + nY - 2,1));
effectD = observed / max(pooledStd,eps);
end
