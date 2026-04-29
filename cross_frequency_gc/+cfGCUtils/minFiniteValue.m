function value = minFiniteValue(x)
% minFiniteValue Return the minimum finite value in an array.
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = min(x);
end
end
