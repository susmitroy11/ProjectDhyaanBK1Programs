function cLims = resolveCommonCLim(matrixList,inputCLim,useSymmetricRange)
% resolveCommonCLim Compute shared color limits across one or more matrices.
if ~exist('useSymmetricRange','var') || isempty(useSymmetricRange)
    useSymmetricRange = 0;
end

if ~isempty(inputCLim)
    cLims = inputCLim;
    return;
end

if ~iscell(matrixList)
    matrixList = {matrixList};
end

if useSymmetricRange
    cMax = NaN;
    for iMatrix = 1:numel(matrixList)
        thisMax = cfGCUtils.maxFiniteValue(abs(matrixList{iMatrix}(:)));
        if ~isfinite(cMax) || (isfinite(thisMax) && thisMax > cMax)
            cMax = thisMax;
        end
    end
    cLims = [-cMax cMax];
else
    cMax = NaN;
    for iMatrix = 1:numel(matrixList)
        thisMax = cfGCUtils.maxFiniteValue(matrixList{iMatrix}(:));
        if ~isfinite(cMax) || (isfinite(thisMax) && thisMax > cMax)
            cMax = thisMax;
        end
    end
    cLims = [0 cMax];
end

if ~all(isfinite(cLims))
    cLims = [];
end
end
