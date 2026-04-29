function [cachedResult,wasLoaded] = tryLoadCachedAnalysis(summaryFileName,requireOutputFiles)
% tryLoadCachedAnalysis Reuse a saved analysisResult when cached outputs exist.
cachedResult = [];
wasLoaded = false;

if ~exist(summaryFileName,'file')
    return;
end

tmp = load(summaryFileName,'analysisResult');
if ~isfield(tmp,'analysisResult')
    return;
end

cachedResult = tmp.analysisResult;
if ~requireOutputFiles
    wasLoaded = true;
    return;
end

if ~isfield(cachedResult,'outputManifest') || isempty(cachedResult.outputManifest)
    return;
end

manifest = cachedResult.outputManifest;
if ischar(manifest)
    manifest = {manifest};
end

for iFile = 1:numel(manifest)
    if ~exist(manifest{iFile},'file')
        return;
    end
end

wasLoaded = true;
end
