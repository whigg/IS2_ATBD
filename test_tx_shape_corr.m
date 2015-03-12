function varargout=test_fpb_corr(varargin)

if nargout>0
    [varargout{1:nargout}]=feval(varargin{:});
else
    feval(varargin{:});
end

%------------------------------------------------------
function [R, HW, bar_out, med_out]=run_sim( WF)
% run a tx-pulse-shape correction simulation for a range of surface roughnesses (R).
% output values give the median and mean height offset.
%
% inputs: 
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
%   med: medians of output values, struct with same fields as bar

% use a moderately-unsaturated strong WF.     
N_chan=12;
N_per_pulse=N_chan*1.25;
R=0:.125:1;
HW=[2:.5:4];
[R, HW]=meshgrid(R, HW);
% aligh the transmit waveform with time zero at its centroid.
dt_WF=(WF.t(2)-WF.t(1));
WF.t=WF.t-sum(WF.t.*WF.p)./sum(WF.p);
 
% N.B.  This can be run as a parfor if you have the parallel tool box
parfor kk=1:numel(R)
     
        fprintf(1, 'R=%d, HW=%3.1d\n', R(kk), HW(kk));
        
        % run this for enough iterations to beat down the noise in the means and medians caused by the spread of the RX pulse
        N_pulses=floor(5e3*(R(kk).^2+.24^2)/(.24^2));
        % Generate the fake data
        [D2, params]=make_data(N_pulses, N_chan, R(kk), N_per_pulse, WF);
        
        % test the correction
        [D3, D3u, D3a]=calc_h(D2, HW(kk), params);
        % save output : calculate the median of each output parameter
        temp_med=median(D3.med);
        temp_centroid=median(D3.centroid);
        temp_med_uncorr=median(D3u.med);
        temp_centroid_uncorr=median(D3u.centroid);
        med(kk)=struct('med', temp_med,'centroid',temp_centroid,'centroid_uncorr', temp_centroid_uncorr, 'med_uncorr', temp_med_uncorr);
       
        % save output : calculate the mean of each output parameter
        temp_med=mean(D3.med-D3a.med);
        temp_centroid=mean(D3.centroid);
        temp_med_uncorr=mean(D3u.med);
        temp_centroid_uncorr=mean(D3u.centroid);
        bar(kk)=struct('med', temp_med,'centroid',temp_centroid,'centroid_uncorr', temp_centroid_uncorr, 'med_uncorr', temp_med_uncorr);
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
function [D3, D3u, D3a]=calc_h(D2, H_HW, params)

% calculate h from received-pulse data.  
% inputs:
% D2: structure array giving PE heights and pulse numbers
% H_HW: surface-window half-height, in meters, used to truncate the data
% params: a structured array that includes (at minimum) the transmit-waveform shape and the background rate
%
% outputs:
% D3: corrected heights,based on the truncated PE distribution,  including the TX-shape correction
% D3u: Uncorrected heights based on the truncated PE distribution not including the TX-shape correction
% D3a: heights based on all of the PE, with no correction and no truncation.

segs=1:58:max(D2.pulse);
[D3a.med, D3a.centroid, D3a.count]=deal(NaN(length(segs)-1,1));
D3=D3a;
[D3.sigma_med, D3.sigma_centroid]=deal(NaN(length(segs)-1,1));
D3u=D3;

for k=1:length(segs)-1;
    in_bin=D2.pulse>=segs(k) & D2.pulse <=segs(k+1);
    D3a.med(k)=median(D2.h(in_bin));
    D3a.centroid(k)=mean(D2.h(in_bin));
    D3a.count(k)=sum(in_bin);
    
    % take an initial subset of pulses in the segment
    D2sub0=index_struct(D2, in_bin);
    
    % iterate 5x to find the truncated median:
    Hmed=median(D2sub0.h);
    in_HW=abs(D2sub0.h-Hmed)<H_HW;
    for k_it=1:5;
        Hmed=median(D2sub0.h(in_HW));
        in_HW=abs(D2sub0.h-Hmed)<H_HW;
    end
    D2sub=index_struct(D2sub0, in_HW);
    
    % calc the uncorrected statistics;
    D3u.med(k)=median(D2sub.h);
    D3u.centroid(k)=mean(D2sub.h);
    D3u.count(k)=length(D2sub.h);
    
    % correct the data
    [t_bin, bin_count]=quick_hist(D2sub.h/(-1.5e8), 1e-12);
    [dM, dCtr]=correct_for_TX_shape(t_bin, bin_count, params.WF.t, params.WF.p, 2*H_HW/1.5e8, params.NoiseRate, params.N_per_pulse*58);
    D3.med(k)=D3u.med(k)+dM;
    D3.centroid(k)=D3u.centroid(k)+dCtr;
end


%-----------------------------------------------------------
function [bin_ctr, count]=quick_hist(x, dx, x0)

% fast algorithm for generating histograms

if ~exist('x0','var');
    x0=0;
end
xb=sort(round((x-x0)/dx));
bin_num=xb-xb(1)+1;

[uB, ind1]=unique(bin_num,'first');
[uB, ind2]=unique(bin_num,'last');
bin_ctr=(xb(1):xb(end))*dx+x0;
count=zeros(size(bin_ctr));
count(uB)=ind2-ind1+1;




