function count=my_histc(x, bins)
x=sort(x);
% case 1: bins are of equal spacing
 if abs(max(diff(bins))-min(diff(bins)))/mean(diff(bins)) > 10*eps
    ind=round((x-bins(1))/(bins(2)-bins(1)))+1;
    [bin_num, ia]=unique(ind,'first');
    [~, ib]=unique(ind,'last');
    count=zeros(size(bins));
    count(bin_num)=ib-ia+1;
else
    count=zeros(size(bins));
    bin_edges=[bins(1)-(bins(2)-bins(1))/2; 
        (bins(2:end)+bins(1:end-1))/2; 
        bins(end)+(bins(end)-bins(end-1))/2];
    for k=1:length(bin_edges)-1;
        count=sum(x>=bin_edges(k) & x <bin_edges(k+1));
    end
end
    
