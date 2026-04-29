function hashText = makeOptionsHash(textIn)
% makeOptionsHash Generate a short repeatable hash for cached analyses.
textIn = char(textIn);

try
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(textIn));
    hash = typecast(md.digest,'uint8');
    hashText = lower(reshape(dec2hex(hash)',1,[]));
catch
    hashText = sprintf('fallback_%u',sum(double(textIn)));
end
end
