function cMap = blueWhiteRed(n)
% blueWhiteRed Create a diverging blue-white-red colormap.
if ~exist('n','var') || isempty(n)
    n = 256;
end

if mod(n,2) == 0
    nLeft = n/2;
    nRight = n/2;
else
    nLeft = floor(n/2);
    nRight = nLeft + 1;
end

blue = [0.1922 0.2118 0.5843];
white = [1 1 1];
red = [0.6471 0 0.1490];

left = [linspace(blue(1),white(1),nLeft)' ...
        linspace(blue(2),white(2),nLeft)' ...
        linspace(blue(3),white(3),nLeft)'];
right = [linspace(white(1),red(1),nRight)' ...
         linspace(white(2),red(2),nRight)' ...
         linspace(white(3),red(3),nRight)'];

if nLeft > 0
    cMap = [left; right(2:end,:)];
else
    cMap = right;
end

if size(cMap,1) ~= n
    xOld = linspace(0,1,size(cMap,1));
    xNew = linspace(0,1,n);
    cMap = interp1(xOld,cMap,xNew);
end
end
