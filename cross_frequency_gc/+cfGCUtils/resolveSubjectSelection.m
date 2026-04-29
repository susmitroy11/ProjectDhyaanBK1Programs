function [subjectNameList,selectionLabel] = resolveSubjectSelection(subjectSelection,dataFolder,protocolNameList,badEyeCondition,badTrialVersion)
% resolveSubjectSelection Validate subject input against available CF-GC files.
protocolNameList = cfGCUtils.resolveProtocolNameList(protocolNameList);

if nargin < 1 || isempty(subjectSelection)
    subjectSelection = 'all';
end

if (ischar(subjectSelection) || (isstring(subjectSelection) && isscalar(subjectSelection))) && strcmpi(char(subjectSelection),'all')
    dirInfo = dir(dataFolder);
    candidateList = {dirInfo([dirInfo.isdir]).name};
    candidateList = setdiff(candidateList,{'.','..','groupStats','figures'});
    candidateList = sort(candidateList);
    subjectNameList = {};

    for iSubject = 1:numel(candidateList)
        if hasCompleteResults(candidateList{iSubject},protocolNameList,dataFolder,badEyeCondition,badTrialVersion)
            subjectNameList{end+1} = candidateList{iSubject}; %#ok<AGROW>
        end
    end
    selectionLabel = 'all_subjects';
else
    subjectNameList = cfGCUtils.normalizeList(subjectSelection);
    for iSubject = 1:numel(subjectNameList)
        if ~hasCompleteResults(subjectNameList{iSubject},protocolNameList,dataFolder,badEyeCondition,badTrialVersion)
            error('Subject %s does not have complete PRE/POST CF-GC files for the requested protocol list.',subjectNameList{iSubject});
        end
    end
    selectionLabel = 'selected_subjects';
end

if isempty(subjectNameList)
    error('No subjects were available for the requested visualization setup.');
end
end

function tf = hasCompleteResults(subjectName,protocolNameList,dataFolder,badEyeCondition,badTrialVersion)
tf = true;
for iProtocol = 1:numel(protocolNameList)
    preFile = fullfile(dataFolder,subjectName,[protocolNameList{iProtocol} '_pre_' badEyeCondition '_' badTrialVersion '_crossFreqPhaseAmpGC.mat']);
    postFile = fullfile(dataFolder,subjectName,[protocolNameList{iProtocol} '_post_' badEyeCondition '_' badTrialVersion '_crossFreqPhaseAmpGC.mat']);
    if ~exist(preFile,'file') || ~exist(postFile,'file')
        tf = false;
        return;
    end
end
end
