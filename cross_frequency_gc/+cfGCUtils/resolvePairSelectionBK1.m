function [pairList,selectionLabel] = resolvePairSelectionBK1(pairSelectionIn,projectRoot)
% resolvePairSelectionBK1 Normalize user input into a BK1 matched-pair list.
matchedPairList = cfGCUtils.getMatchedSubjectPairsBK1(projectRoot);

if nargin < 1 || isempty(pairSelectionIn) || ...
        ((ischar(pairSelectionIn) || (isstring(pairSelectionIn) && isscalar(pairSelectionIn))) && strcmpi(char(pairSelectionIn),'all'))
    pairList = matchedPairList;
    selectionLabel = 'all_pairs';
    return;
end

if isnumeric(pairSelectionIn)
    pairIndexList = pairSelectionIn(:)';
    validatePairIndices(pairIndexList,numel(matchedPairList));
    pairList = matchedPairList(pairIndexList);
    selectionLabel = 'selected_pair_indices';
    return;
end

if iscell(pairSelectionIn) && numel(pairSelectionIn) == 2 && all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)),pairSelectionIn))
    pairList = resolveByNestedCells({pairSelectionIn},matchedPairList);
    selectionLabel = 'selected_pairs';
    return;
end

if iscell(pairSelectionIn)
    pairList = resolveByNestedCells(pairSelectionIn,matchedPairList);
    selectionLabel = 'selected_pairs';
    return;
end

error('pair selection must be ''all'', numeric pair indices, a single pair {meditator, control}, or a cell array of such pairs.');
end

function pairList = resolveByNestedCells(inputCell,matchedPairList)
pairList = repmat(matchedPairList(1),0,1);
for iPair = 1:numel(inputCell)
    pairSpec = inputCell{iPair};
    if ~(iscell(pairSpec) && numel(pairSpec) == 2)
        error('Each pair specification must be a 1x2 cell: {meditatorSubject, controlSubject}.');
    end

    meditatorName = char(pairSpec{1});
    controlName = char(pairSpec{2});
    matchIndex = find(strcmp({matchedPairList.meditator},meditatorName) & strcmp({matchedPairList.control},controlName),1);
    if isempty(matchIndex)
        error('Pair %s vs %s was not found in BK1 matched pairs.',meditatorName,controlName);
    end
    pairList(end+1,1) = matchedPairList(matchIndex); %#ok<AGROW>
end
end

function validatePairIndices(pairIndexList,numPairs)
if any(pairIndexList < 1) || any(pairIndexList > numPairs) || any(pairIndexList ~= round(pairIndexList))
    error('Pair indices must be integers between 1 and %d.',numPairs);
end
end
