function token = sanitizePathToken(tokenIn)
% sanitizePathToken Turn free text into a filesystem-safe token.
token = char(tokenIn);
token = regexprep(token,'[^a-zA-Z0-9_\-]+','_');
token = regexprep(token,'_+','_');
token = regexprep(token,'(^_)|(_$)','');
if isempty(token)
    token = 'default';
end
end
