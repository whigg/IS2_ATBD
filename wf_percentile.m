function t0=wf_percentile(t, z, p)

C=cumsum(z)/sum(z);
t0=NaN(size(p));
for k=1:length(p);
    samp0=find(C<p(k), 1, 'last');
    samp1=find(C>p(k), 1, 'first');
    
    if isempty(samp0) && ~isempty(samp1)
        t0(k)=t(samp1); continue
    end
    if isempty(samp1) && ~isempty(samp0)
        t0(k)=t(samp0); continue;
    end
    if isempty(samp0) && isempty(samp1);
        t0(k)=NaN; continue
    end
    
    t0(k)=interp1(C([samp0, samp1]), t([samp0,samp1]), p(k));
end