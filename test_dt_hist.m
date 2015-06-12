function [R, dt_hist, bar_out, med_out, bar_0, med_0, D3_corr_1, D3_uncorr_1, D3_no_fpb_1]=test_dt_hist(N_chan, WF)
% run a first-photon-bias simulation for an N-channel detector, for a
% range of signal strengths (N_per_pulse), and surface roughnesses (R).
% output values give the median and mean height offset.
%
% inputs:
%    N_chan: the number of channels in the detector
%    WF:   The transmit pulse.  struct: t: time, p: power.  t should be offset so that the waveform centroid is zero
%
% outputs:
%   R: the roughness values for which the simulation was calculated, m
%   N_per_pulse: the signal strengths (photons per pulse) for the simulation
%   bar: means of output values, struct with fields:
%      med: corrected median photon height
%      centroid: corrected centroid heights
%      med_uncorr: uncorrected median photon heights
%      centroid_uncorr: uncorrected centroid photon heights
%      N: corrected photon counts
%      N_uncorr: uncorrected photon counts
%   med: medians of output values, struct with same fields as bar

% to display the output:
%
%    subplot(1,2,1);  pcolor(R(1,:), N_per_pulse(:,1),  med.med_uncorr); xlabel('Roughness'); ylabel('PE per pulse'); h_cb=colorbar; ylabel(h_cb,'uncorrected height error, m');
%    subplot(1,2,2);  pcolor(R(1,:), N_per_pulse(:,1),  med.med); xlabel('Roughness'); ylabel('PE per pulse'); h_cb=colorbar; ylabel(h_cb,'corrected height error, m');
%

N_per_pulse=2*N_chan*.75;
R=[0 .25 1]; R=0.25;
dt_hist=logspace(-1, 1, 30)*100e-12;

[dt_hist, R]=ndgrid(dt_hist, R);

dt_WF=(WF.t(2)-WF.t(1));

WF.t=WF.t-sum(WF.t.*WF.p)./sum(WF.p);

% time period within which the WF is truncated- pick 4 m so that not too much of the broadest WF gets cut
H_trunc=4;

% N.B.  This can be run as a parfor if you have the parallel tool box
for kR=1:size(R,2)
    if R(kR) > 0
        % generate a synthetic waveform that matches the spread transmit pulse
        t_r=R(kR)/1.5e8;
        n_G=ceil(4*t_r/dt_WF);
        t_G=(-n_G:n_G)*dt_WF;
        G=gaussian(t_G, 0, t_r);
        WF1.p=conv(WF.p, G);
        WF1.t=[(WF.t(1)+(-n_G:-1)*dt_WF)'; WF.t; (WF.t(end)+(1:n_G)*dt_WF)'];
    else
        WF1=WF;
    end
    
    fprintf(1, 'R=%d, N=%d\n', R(kR), N_per_pulse);    
    % truncate the synthetic WF (speeds up the calculation)
    els=abs(WF1.t*1.5e8)<H_trunc;
    WF1.t=WF1.t(els); WF1.p=WF1.p(els);
    % run this for enough iterations to beat down the noise in the means and medians caused by the spread of the RX pulse
    N_pulses=floor(3e4*(R(kR).^2+.24^2)/(.24^2));
    % calculate the mean and median of the RX pulse
    med_0(kR)=-1.5e8*wf_median(WF1.t, WF1.p);
    bar_0(kR)=-1.5e8*sum(WF1.t.*WF1.p)./sum(WF1.p);
    [D2, params]=test_fpb_corr('make_data', N_pulses, N_chan, R(kR), N_per_pulse, WF);
    D2=index_struct(D2, abs(D2.h)<H_trunc);
            
    for kT=1:size(dt_hist,1);
        params1=params;
        params1.dt_hist_bin=dt_hist(kT, kR);
        % truncate the fake data by H_trunc (speeds up the calculation)
        
        % run correction
        [D3_corr, D3_uncorr, D3_no_fpb]=test_fpb_corr('correct_data', D2, params1);
        [D3_corr_1(kT, kR), D3_uncorr_1(kT, kR), D3_no_fpb_1(kT, kR)]=deal(D3_corr, D3_uncorr, D3_no_fpb);

        % save output : calculate the median of each output parameter
        temp_med=median(D3_corr.med-D3_no_fpb_1(kT, 1).med);
        temp_centroid=median(D3_corr.centroid-D3_no_fpb_1(kT, 1).centroid);
        temp_med_uncorr=median(D3_uncorr.med-D3_no_fpb_1(kT, 1).med);
        temp_centroid_uncorr=median(D3_uncorr.centroid-D3_no_fpb_1(kT, 1).centroid);
        temp_N=median(D3_corr.count);
        temp_N_uncorr=median(D3_uncorr.count);
        med(kT, kR)=struct('med', temp_med,'centroid',temp_centroid,'centroid_uncorr', temp_centroid_uncorr, 'med_uncorr', temp_med_uncorr,'N', temp_N,'N_uncorr', temp_N_uncorr);
        
        % save output : calculate the mean of each output parameter
        temp_med=mean(D3_corr.med-D3_no_fpb_1(kT, 1).med);
        temp_centroid=mean(D3_corr.centroid-D3_no_fpb_1(kT, 1).centroid);
        temp_med_uncorr=mean(D3_uncorr.med-D3_no_fpb_1(kT, 1).med);
        temp_centroid_uncorr=mean(D3_uncorr.centroid-D3_no_fpb_1(kT, 1).centroid);
        temp_N=mean(D3_corr.count);
        temp_N_uncorr=mean(D3_uncorr.count);
        bar(kT,kR)=struct('med', temp_med,'centroid',temp_centroid,'centroid_uncorr', temp_centroid_uncorr, 'med_uncorr', temp_med_uncorr,'N', temp_N,'N_uncorr', temp_N_uncorr);
    end
end
% collect the output:

f=fieldnames(bar);
for kf=1:length(f);
    med_out.(f{kf})=reshape([med.(f{kf})], size(R));
    bar_out.(f{kf})=reshape([bar.(f{kf})], size(R));
end

for k=1:length(dt_hist); 
    m(k)=mean(D3_corr_1(k).med)-med_0; 
    c(k)=mean(D3_corr_1(k).centroid)-bar_0; 
end

figure; plot(dt_hist/1e-12, c, dt_hist/1e-12, m); legend('centroid','median'); xlabel('histogram bin size, ps'); ylabel('error bias, m');

