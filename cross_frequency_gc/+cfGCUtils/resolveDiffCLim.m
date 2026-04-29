function cLims = resolveDiffCLim(matrixList,inputCLim)
% resolveDiffCLim Compute symmetric color limits for difference maps.
cLims = cfGCUtils.resolveCommonCLim(matrixList,inputCLim,1);
end
