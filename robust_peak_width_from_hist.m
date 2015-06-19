function [sigma_hat, med]=robust_peak_width_from_hist(x, z, N_BG, XR)

% C_tot=C_sig+C_noise
%      =C_sig(t)+(t-t0)BGR
% want t1 such that C_sig(t1)=f
% Ctot(t1)=f+(t1-t0)BGR
% look for first point where Ctot > 0.16 Nsig + (t1-t0)*BGR 
% and last point where Ctot < 0.84 Nsig + (t1-t0)*BGR
% the difference gives the peak width

if ~exist('TR','var')
    XR=[min(x) max(x)];
end
BGR=N_BG/diff(XR);

N_sig=sum(z)-N_BG;

C=cumsum(z); 

i0=find(C<0.16*N_sig + (x-XR(1))*BGR, 1, 'last'); 
if isempty(i0); i0=1;end
i1=find(C>0.84*N_sig + (x-XR(1))*BGR, 1, 'first');
if isempty(i1); i1=length(x);end
sigma_hat=diff(x([i0 i1]))/2;

if nargout==2;
    i0=find(C<0.5*N_sig + (x-XR(1))*BGR, 1, 'last');
    if isempty(i0); i0=1;end
    i1=find(C>0.5*N_sig + (x-XR(1))*BGR, 1, 'first');
    if isempty(i1); i1=length(x);end
    
    med=interp1(C([i0 i1])-(x([i0 i1])-XR(1))*BGR, x([i0 i1]), 0.5*N_sig);
end

