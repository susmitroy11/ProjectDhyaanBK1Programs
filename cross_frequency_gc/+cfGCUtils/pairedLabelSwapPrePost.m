function [pValue,effectDz,nUsed] = pairedLabelSwapPrePost(preValues,postValues,numPermutations)
% pairedLabelSwapPrePost Permutation test for paired pre/post changes.
validMask = isfinite(preValues) & isfinite(postValues);
preValues = preValues(validMask);
postValues = postValues(validMask);
preValues = preValues(:);
postValues = postValues(:);
nUsed = numel(preValues);

if nUsed < 2
    pValue = NaN;
    effectDz = NaN;
    return;
end

diffValues = postValues - preValues;
observed = mean(diffValues);
permStats = nan(numPermutations,1);

for iPerm = 1:numPermutations
    swapMask = rand(nUsed,1) > 0.5;
    prePerm = preValues;
    postPerm = postValues;
    prePerm(swapMask) = postValues(swapMask);
    postPerm(swapMask) = preValues(swapMask);
    permStats(iPerm) = mean(postPerm - prePerm);
end

pValue = (sum(abs(permStats) >= abs(observed)) + 1) / (numPermutations + 1);
effectDz = observed / max(std(diffValues,0),eps);
end
