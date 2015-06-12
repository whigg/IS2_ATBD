function varargout=test_fpb_corr(varargin)

if nargout>0
    [varargout{1:nargout}]=feval(varargin{:});
else
    feval(varargin{:});
end

%------------------------------------------------------
function [R, N_per_pulse, bar_out, med_out]=run_sim(N_chan, WF)
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

N_per_pulse=2.^linspace(0, 1, 10)*N_chan*.75;
R=0:.125:1;
[N_per_pulse, R]=ndgrid(N_per_pulse, R);

dt_WF=(WF.t(2)-WF.t(1));

WF.t=WF.t-sum(WF.t.*WF.p)./sum(WF.p);

% time period within which the WF is truncated- pick 4 m so that not too much of the broadest WF gets cut
H_trunc=4;

% N.B.  This can be run as a parfor if you have the parallel tool box
parfor kk=1:numel(R)
        if R(kk) > 0
            % generate a synthetic waveform that matches the spread transmit pulse
            t_r=R(kk)/1.5e8;
            n_G=ceil(4*t_r/dt_WF);
            t_G=(-n_G:n_G)*dt_WF;
            G=gaussian(t_G, 0, t_r);
            WF1.p=conv(WF.p, G);
            WF1.t=[(WF.t(1)+(-n_G:-1)*dt_WF)'; WF.t; (WF.t(end)+(1:n_G)*dt_WF)'];           
        else
            WF1=WF;
        end
        
        fprintf(1, 'R=%d, N=%d\n', R(kk), N_per_pulse(kk));
        
        % truncate the synthetic WF (speeds up the calculation)
        els=abs(WF1.t*1.5e8)<H_trunc;
        WF1.t=WF1.t(els); WF1.p=WF1.p(els);
        % run this for enough iterations to beat down the noise in the means and medians caused by the spread of the RX pulse
        N_pulses=floor(5e3*(R(kk).^2+.24^2)/(.24^2));
        % calculate the mean and median of the RX pulse
        med_0=-1.5e8*wf_median(WF1.t, WF1.p);
        bar_0=-1.5e8*sum(WF1.t.*WF1.p)./sum(WF1.p);
        [D2, params]=make_data(N_pulses, N_chan, R(kk), N_per_pulse(kk), WF);        
        params.dt_hist_bin=100e-12;
        % truncate the fake data by H_trunc (speeds up the calculation)
        D2=index_struct(D2, abs(D2.h)<H_trunc);
        
        % run correction
        [D3_corr, D3_uncorr, D3_no_fpb]=correct_data(D2, params);
        % save output : calculate the median of each output parameter
        temp_med=median(D3_corr.med-D3_no_fpb.med);
        temp_centroid=median(D3_corr.centroid-D3_no_fpb.centroid);
        temp_med_uncorr=median(D3_uncorr.med-D3_no_fpb.med);
        temp_centroid_uncorr=median(D3_uncorr.centroid-D3_no_fpb.centroid);
        temp_N=median(D3_corr.count);
        temp_N_uncorr=median(D3_uncorr.count);      
        med(kk)=struct('med', temp_med,'centroid',temp_centroid,'centroid_uncorr', temp_centroid_uncorr, 'med_uncorr', temp_med_uncorr,'N', temp_N,'N_uncorr', temp_N_uncorr);
       
        % save output : calculate the mean of each output parameter
        temp_med=mean(D3_corr.med-D3_no_fpb.med);
        temp_centroid=mean(D3_corr.centroid-D3_no_fpb.centroid);
        temp_med_uncorr=mean(D3_uncorr.med-D3_no_fpb.med);
        temp_centroid_uncorr=mean(D3_uncorr.centroid-D3_no_fpb.centroid);
        temp_N=mean(D3_corr.count);
        temp_N_uncorr=mean(D3_uncorr.count);          
        bar(kk)=struct('med', temp_med,'centroid',temp_centroid,'centroid_uncorr', temp_centroid_uncorr, 'med_uncorr', temp_med_uncorr,'N', temp_N,'N_uncorr', temp_N_uncorr);
%    end
%end
end
% collect the output:

f=fieldnames(bar);
for kf=1:length(f); 
    med_out.(f{kf})=reshape([med.(f{kf})], size(R));
    bar_out.(f{kf})=reshape([bar.(f{kf})], size(R));
end


%-------------------------------------------------------
function [D2, params]=make_data(N_pulses, N_chan, roughness, N_per_pulse, WF);

DEM=0;
x0=zeros(N_pulses,1);
ATM_xmit=ones(size(x0));
params=struct('roughness', roughness,   'sigma_x', 7.5, 'NoiseRate', 0, 'H_window', 0, 'WF', WF, 'c', 3e8, 'N_channels', N_chan, 'N_per_pulse', N_per_pulse,'t_dead', 3.2e-9); 
D2=det_sim(DEM, x0, params, ATM_xmit);
D2.pulse=D2.pulse_num;
D2=rmfield(D2, 'pulse_num');

%-------------------------------------------------------
function [D3_corr, D3_uncorr, D3_no_fpb]=correct_data(D2, params)

segs=1:58:max(D2.pulse);
[D3_no_fpb.med, D3_no_fpb.centroid, D3_no_fpb.count]=deal(NaN(length(segs)-1,1));
D3_corr=D3_no_fpb;
[D3_corr.sigma_med, D3_corr.sigma_centroid]=deal(NaN(length(segs)-1,1));
D3_uncorr=D3_corr;

for k=1:length(segs)-1;
    in_bin=D2.pulse>=segs(k) & D2.pulse <=segs(k+1);
    D3_no_fpb.med(k)=median(D2.h(in_bin));
    D3_no_fpb.centroid(k)=mean(D2.h(in_bin));
    D3_no_fpb.count(k)=sum(in_bin);
    
    these=D2.detected & in_bin;
    D2sub=index_struct(D2, these);  
    [D3_corr.med(k), D3_corr.centroid(k), D3_corr.count(k), t_WF, N0, N_fpb_corr, D3_corr.sigma_med(k), D3_corr.sigma_centroid(k)]=fpb_corr_ATBD(D2sub.h, D2sub.channel, D2sub.pulse, params.N_channels,57, params.t_dead, false, params.dt_hist_bin);
    [D3_uncorr.med(k), D3_uncorr.centroid(k), D3_uncorr.count(k), t_WF, N0, N_fpb_corr, D3_uncorr.sigma_med(k), D3_uncorr.sigma_centroid(k)]=fpb_corr_ATBD(D2sub.h, D2sub.channel, D2sub.pulse, params.N_channels, 57, params.t_dead, true, params.dt_hist_bin);

end
