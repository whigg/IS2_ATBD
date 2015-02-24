function [med, centroid, count, t_WF, N0, N_fpb_corr, sigma_med, sigma_centroid]=fpb_corr_ATBD(dh, chan, pulse, N_chan, N_pulses, t_dead, skip_fpb_corr)
c=3e8;

t=-2/c*dh;

dt=2e-12; % approx 1/3  mm
t_WF=min(t):dt:(max(t)+1e-9);
t_WF=t_WF(:);

tbin=floor((t-t_WF(1))/dt)+2; % the +2 is: 1 for 1-indexed array, 1 for the timing of the event vs. the time of the bin  
good=tbin>0 & tbin < length(t_WF);
tbin=tbin(good);
t=t(good);
chan=chan(good);
pulse=pulse(good);
pulse=pulse-min(pulse)+1;

pd=pulse+1i*chan;
upd=unique(pd);
% calc the signal rate, in photons per pulse per chan per time bin as a fn of time
N0=histc(t, t_WF)/N_pulses/N_chan;  

N_dt_bins=floor(t_dead/dt)-1;

% note: could speed up the calculation by calculating:
%SigRatePre=conv(N0, [zeros(1, N_dt_bins), ones(1, N_dt_bins)]/N_dt_bins/dt,'same')
% the gain is only likely to differ from 1 when SigRatePre > 0.05 /t_dead

ind=1:length(t);

N=N0;
els=cell(length(upd),1);
if ~skip_fpb_corr
    for kk=1:3;
        dead_all=zeros(length(t_WF), N_pulses*N_chan);
        for kpd=1:length(upd);  % loop over pulses and channels
            if kk==1
                els{kpd}=ind(pd==upd(kpd));
            end
            these=els{kpd};
            dead=false(size(t_WF));
            P_dead=zeros(size(t_WF));
            dead_bins=cell(length(these), 1);
            for kp=1:length(these);  %assign definitely-dead and definitely-not-still-dead times
                dead_bins{kp}=tbin(these(kp)):min(length(t_WF), tbin(these(kp))+N_dt_bins);
                dead(dead_bins{kp})=true;
            end
            % now loop again and assign probabilities of dead bins from extended deadtime for bins
            % after the end of the dead time of each detected photon
            for kp=1:length(these);
                E_hits=N(dead_bins{kp}); % expected number of hits during dead time
                % the probability that the channel is still dead at time dt after the
                % end of the dead time is the probability of getting at least 1 hit in
                % the previous t_dead
                % or 1 - P(zero hits in prev. dt)
                %
                % take the cumsum of the reverse of E_hits
                %
                % note that if the expected number of hits is E_hits, the probability of no hits in a bin is 0^k exp(-E_hits) / 0! or exp(-E_hits).
                %this gives this_P_dead = 1-exp(cumsum(log(exp(-E_hits(end:-1:1)) or...
                this_P_dead = 1-exp(cumsum(-E_hits(end:-1:1)));
                this_P_dead = this_P_dead(end:-1:1);
                out_bins=((tbin(these(kp))+N_dt_bins+1):min(tbin(these(kp))+1+2*N_dt_bins, length(t_WF)));
                if ~isempty(out_bins)
                    P_dead(out_bins)=max(P_dead(out_bins), this_P_dead(1:length(out_bins)));
                end
            end
            
            P_dead(dead)=1; %we know when the channel was definitely dead.
            dead_all(:, real(upd(kpd))*N_chan + imag(upd(kpd)))=P_dead;
        end
        
        gain = 1-mean(dead_all,2);
        N= N0./gain;
    end
else
    gain=ones(size(N0));
    N=N0;
end

N_fpb_corr=N0./gain*N_pulses*N_chan;
[med, centroid, count, sigma_med, sigma_centroid]=calc_stats(N0*N_pulses*N_chan, gain, t_WF*(-1.5e8) );

%----------------------------------------------------
function [med, centroid, N, sigma_med, sigma_centroid]=calc_stats(WF, gain, t_WF)

WFc=WF(:)./gain(:);
centroid=sum(t_WF.*WFc)./sum(WFc);
N=sum(WF./gain);

t_40_50_60=percentile_of_histogram([0.4 0.5 0.6], t_WF, WF./gain);
med=t_40_50_60(2);

sigma_WF=sqrt(WF)./gain;
sigma_CDF_med=sqrt([0.5*sum(sigma_WF.^2)./N^2]);
sigma_med=abs(diff(t_40_50_60([1 3])))/0.2*sigma_CDF_med;

sigma_centroid=sqrt(sum((WF./gain.*(t_WF-centroid)/N).^2));


%----------------------------------------------------
function X=percentile_of_histogram(P, bins, counts)
bin_width=[diff(bins(:)); bins(end)-bins(end-1)];
edges=[bins(1)-bin_width(1)/2; bins(:)+bin_width/2];
C=[0; cumsum(counts(:))]; C=C/C(end);
for k=1:length(P)
    i_minus=find(C<=P(k), 1, 'last');
    i_plus=find(C>=P(k), 1, 'first');
    if C(i_plus)==C(i_minus)
        X(k)=0.5*(edges(i_plus)+edges(i_minus));
    else
        X(k)=((P(k)-C(i_minus))*edges(i_plus)+(C(i_plus)-P(k))*edges(i_minus))/(C(i_plus)-C(i_minus));
    end
end

