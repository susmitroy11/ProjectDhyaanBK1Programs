function out = normalizeList(in)
% normalizeList Convert common text inputs to a row cell array of char vectors.
if nargin == 0 || isempty(in)
    out = {};
    return;
end

if ischar(in)
    out = {in};
elseif isstring(in)
    out = cellstr(in(:));
elseif iscell(in)
    out = in(:)';
else
    out = cellstr(in);
end

for iItem = 1:numel(out)
    if isstring(out{iItem})
        out{iItem} = char(out{iItem});
    end
end
end
