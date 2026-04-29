function [qValues,sigMask] = bhFDRByPair(pValues,alpha)
% bhFDRByPair Apply BH-FDR separately within each phase/amplitude band pair.
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.05;
end

qValues = nan(size(pValues));
sigMask = false(size(pValues));

for iPhase = 1:size(pValues,3)
    for iAmp = 1:size(pValues,4)
        thisP = pValues(:,:,iPhase,iAmp);
        qThis = localBH(thisP(:));
        qThis = reshape(qThis,size(thisP));
        qValues(:,:,iPhase,iAmp) = qThis;
        sigMask(:,:,iPhase,iAmp) = qThis < alpha;
    end
end
end

function qValues = localBH(pVector)
qValues = nan(size(pVector));
validMask = isfinite(pVector);
validP = pVector(validMask);
if isempty(validP)
    return;
end

[sortedP,order] = sort(validP(:));
n = numel(sortedP);
qSorted = nan(n,1);
qSorted(n) = sortedP(n);
for i = n-1:-1:1
    qSorted(i) = min((n/i) * sortedP(i), qSorted(i+1));
end
qSorted = min(qSorted,1);

tmp = nan(n,1);
tmp(order) = qSorted;
qValues(validMask) = tmp;
end
