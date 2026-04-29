function protocolNameList = resolveProtocolNameList(protocolNameListIn)
% resolveProtocolNameList Normalize protocol input and expand 'all'.
defaultProtocolNameList = cfGCUtils.getDefaultProtocolNameList();

if ~exist('protocolNameListIn','var') || isempty(protocolNameListIn)
    protocolNameList = defaultProtocolNameList;
    return;
end

if ischar(protocolNameListIn) || (isstring(protocolNameListIn) && isscalar(protocolNameListIn))
    protocolToken = char(protocolNameListIn);
    if strcmpi(protocolToken,'all')
        protocolNameList = defaultProtocolNameList;
        return;
    end
    protocolNameList = {protocolToken};
elseif isstring(protocolNameListIn)
    protocolNameList = cellstr(protocolNameListIn(:));
elseif iscell(protocolNameListIn)
    protocolNameList = protocolNameListIn(:);
else
    error('protocolNameList must be ''all'', a protocol name, or a cell/string array of protocol names.');
end

for iProtocol = 1:numel(protocolNameList)
    protocolNameList{iProtocol} = char(string(protocolNameList{iProtocol}));
end

[~,uniqueIndices] = unique(protocolNameList,'stable');
protocolNameList = protocolNameList(sort(uniqueIndices));
protocolNameList = protocolNameList(:)';

invalidProtocols = setdiff(protocolNameList,defaultProtocolNameList);
if ~isempty(invalidProtocols)
    error('Unknown protocol(s): %s',strjoin(invalidProtocols,', '));
end
end
