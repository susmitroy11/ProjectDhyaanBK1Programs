function results = loadConditionResults(subjectName,protocolName,conditionType,badEyeCondition,badTrialVersion,dataFolder)
% loadConditionResults Load one saved subject/protocol/condition CF-GC MAT file.
fileName = fullfile(dataFolder,subjectName,[protocolName '_' conditionType '_' badEyeCondition '_' badTrialVersion '_crossFreqPhaseAmpGC.mat']);
if ~exist(fileName,'file')
    error('Could not find CF-GC file: %s',fileName);
end

tmp = load(fileName,'results');
if ~isfield(tmp,'results') || isempty(tmp.results)
    error('No results structure found in file: %s',fileName);
end

results = tmp.results;
results.loadedFrom = fileName;
end
