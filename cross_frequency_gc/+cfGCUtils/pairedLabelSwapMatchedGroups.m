function [pValue,effectDz,nUsed] = pairedLabelSwapMatchedGroups(xValues,yValues,numPermutations)
% pairedLabelSwapMatchedGroups Permutation test for matched-group differences.
validMask = isfinite(xValues) & isfinite(yValues);
xValues = xValues(validMask);
yValues = yValues(validMask);
xValues = xValues(:);
yValues = yValues(:);
nUsed = numel(xValues);

if nUsed < 2
    pValue = NaN;
    effectDz = NaN;
    return;
end

pairDiff = xValues - yValues;
observed = mean(pairDiff);
permStats = nan(numPermutations,1);

for iPerm = 1:numPermutations
    swapMask = rand(nUsed,1) > 0.5;
    xPerm = xValues;
    yPerm = yValues;
    xPerm(swapMask) = yValues(swapMask);
    yPerm(swapMask) = xValues(swapMask);
    permStats(iPerm) = mean(xPerm - yPerm);
end

pValue = (sum(abs(permStats) >= abs(observed)) + 1) / (numPermutations + 1);
effectDz = observed / max(std(pairDiff,0),eps);
end
