function t0=wf_median(t, z);
p=0.5;
C=cumsum(z)/sum(z);
samp0=find(C<p, 1, 'last');
samp1=find(C>p, 1, 'first');
t0=interp1(C([samp0, samp1]), t([samp0,samp1]), p);