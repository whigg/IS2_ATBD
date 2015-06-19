function  [low, high]=iqr(x, f, skipNaN)

if exist('skipNaN','var') && skipNaN
    x=x(isfinite(x));
end

if nargin==1;
   f=0.68; 
end

if min(size(x))==1 & size(x,2)>size(x,1);
   x=x(:);
end

f=0.5+[-f/2;f/2];

[x,ind]=sort(x);

I=(0.5:size(x,1)-0.5)'/(size(x,1));   % this is now the abscissa for the CDF

low=interp1(I, x, f);

low(f<min(I))=min(x);
low(f>max(I))=max(x);

if nargout==1;
   low=diff(low);
elseif nargout==2;
   high=low(2,:); low=low(1,:);
end
