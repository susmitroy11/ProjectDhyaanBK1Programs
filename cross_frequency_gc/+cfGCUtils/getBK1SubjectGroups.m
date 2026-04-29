function [subjectNameList,groupNameList,unknownSubjects] = getBK1SubjectGroups(projectRoot)
% getBK1SubjectGroups Load BK1 subjects and assign meditator/control labels.
if ~exist('projectRoot','var') || isempty(projectRoot)
    projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

infoFileName = fullfile(projectRoot,'ProjectDhyaanBK1Programs','commonAnalysisCodes','informationFiles','BK1AllSubjectList.mat');
tmp = load(infoFileName);

declaredBadSubjects = {'004P','081SN','069MG','092KB','063VK'};
allSubjectList = setdiff(normalizeList(tmp.allSubjectList),declaredBadSubjects,'stable');
meditatorList = setdiff(normalizeList(tmp.meditatorList),declaredBadSubjects,'stable');
controlList = setdiff(normalizeList(tmp.controlList),declaredBadSubjects,'stable');

subjectNameList = allSubjectList(:);
groupNameList = cell(size(subjectNameList));
unknownSubjects = {};

for iSubject = 1:numel(subjectNameList)
    if ismember(subjectNameList{iSubject},meditatorList)
        groupNameList{iSubject} = 'Meditator';
    elseif ismember(subjectNameList{iSubject},controlList)
        groupNameList{iSubject} = 'Control';
    else
        groupNameList{iSubject} = 'Unknown';
        unknownSubjects{end+1,1} = subjectNameList{iSubject}; %#ok<AGROW>
    end
end

if ~isempty(unknownSubjects)
    warning('cfGCUtils:getBK1SubjectGroups:UnknownSubjects', ...
        'Subjects excluded because group was Unknown: %s', strjoin(unknownSubjects',', '));
end
end

function out = normalizeList(in)
if iscell(in)
    out = in(:);
elseif isstring(in)
    out = cellstr(in(:));
else
    out = cellstr(in);
end
end
