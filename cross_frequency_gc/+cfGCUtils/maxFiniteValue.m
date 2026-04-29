function value = maxFiniteValue(x)
% maxFiniteValue Return the maximum finite value in an array.
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = max(x);
end
end
