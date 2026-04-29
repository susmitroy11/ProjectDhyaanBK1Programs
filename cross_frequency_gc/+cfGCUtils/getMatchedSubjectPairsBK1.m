function matchedPairList = getMatchedSubjectPairsBK1(projectRoot)
% getMatchedSubjectPairsBK1 Load the BK1 matched meditator-control pair list.
if ~exist('projectRoot','var') || isempty(projectRoot)
    projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

infoFolder = fullfile(projectRoot,'ProjectDhyaanBK1Programs','commonAnalysisCodes','informationFiles');
addpath(infoFolder);
cleanupObj = onCleanup(@() rmpath(infoFolder)); %#ok<NASGU>

pairedSubjectNameList = getPairedSubjectsBK1();
numPairs = size(pairedSubjectNameList,1);

matchedPairList = repmat(struct( ...
    'pairIndex',[], ...
    'meditator','', ...
    'control','', ...
    'pairLabel',''), numPairs, 1);

for iPair = 1:numPairs
    matchedPairList(iPair).pairIndex = iPair;
    matchedPairList(iPair).meditator = pairedSubjectNameList{iPair,1};
    matchedPairList(iPair).control = pairedSubjectNameList{iPair,2};
    matchedPairList(iPair).pairLabel = sprintf('%s_vs_%s', ...
        pairedSubjectNameList{iPair,1},pairedSubjectNameList{iPair,2});
end
end
