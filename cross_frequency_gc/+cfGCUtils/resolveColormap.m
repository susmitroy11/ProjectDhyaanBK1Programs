function cMapOut = resolveColormap(cMapIn)
% resolveColormap Resolve a function handle, name, or raw colormap array.
if isa(cMapIn,'function_handle')
    cMapOut = cMapIn(256);
    return;
end

if ischar(cMapIn) || (isstring(cMapIn) && isscalar(cMapIn))
    cMapOut = feval(char(cMapIn),256);
    return;
end

cMapOut = cMapIn;
if isnumeric(cMapOut) && size(cMapOut,2) == 3 && size(cMapOut,1) > 1 && size(cMapOut,1) < 256
    xOld = linspace(0,1,size(cMapOut,1));
    xNew = linspace(0,1,256);
    cMapOut = interp1(xOld,cMapOut,xNew);
end
end
