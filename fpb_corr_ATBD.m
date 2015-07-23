function [med, centroid, count, t_WF_full, N0_full, N_fpb_corr, sigma_med, sigma_centroid, minGain]=fpb_corr_ATBD(dh, chan, pulse, N_chan, N_pulses, t_dead, skip_fpb_corr,dt)
c=3e8;

t=-2/c*dh;

if ~exist('dt','var')
    dt=100e-12; %approx 1.5 cm
end

% calculate a full-waveform time vector:
t_WF_full=(min(t)-dt/2):dt:(max(t)+1e-9);
t_WF_full=t_WF_full(:);

% calc the signal rate, in photons per pulse per chan per time bin as a fn of time
N0_full=my_histc(t, t_WF_full)/N_pulses/N_chan;  

% make sure nothing funny happens at the start and end of the WF
N0_full(1)=0; 
N0_full(end)=0;
N_dt_bins=floor(t_dead/dt)-1;

% calculate the number of Ph per dead time
N_per_dt=conv(N0_full,  [zeros(N_dt_bins,1); ones(N_dt_bins,1)],'same');

% only calculate the gain for bins for which the PH/deadtime rate is greater than 0.1 
if ~any(N_per_dt>0.05)
    % generate a default gain estimate for the full waveform
    % If we haven't calculated a gain value, assume it's equal to 1
    gain_full=ones(size(N0_full));
     
    N_fpb_corr=N0_full./gain_full*N_pulses*N_chan;
    [med, centroid, count, sigma_med, sigma_centroid]=calc_stats(N0_full*N_pulses*N_chan, gain_full, t_WF_full*(-1.5e8) );
    minGain=1;
    return
end
TR=range(t_WF_full(N_per_dt>0.05))+[-1 1]*t_dead;
gain_calc_bins=t_WF_full >=TR(1) & t_WF_full <= TR(2);
t_WF=t_WF_full(gain_calc_bins);
N0=N0_full(gain_calc_bins);

tbin=floor((t-t_WF(1))/dt)+2; % the +2 is: 1 for 1-indexed array, 1 for the timing of the event vs. the time of the bin  
good=tbin>0 & tbin < length(t_WF);
tbin=tbin(good);
t=t(good);
chan=chan(good);
pulse=pulse(good);
pulse=pulse-min(pulse)+1;

pd=pulse+1i*chan;
upd=unique(pd);
 
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

% generate a default gain estimate for the full waveform
% If we haven't calculated a gain value, assume it's equal to 1
gain_full=ones(size(N0_full));
gain_full(gain_calc_bins)=gain;
minGain=min(gain_full);

N_fpb_corr=N0_full./gain_full*N_pulses*N_chan;
[med, centroid, count, sigma_med, sigma_centroid]=calc_stats(N0_full*N_pulses*N_chan, gain_full, t_WF_full*(-1.5e8) );

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

