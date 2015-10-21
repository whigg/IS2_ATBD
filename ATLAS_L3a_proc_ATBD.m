function [D3]=ATLAS_L3a_proc_ATBD(D2, params, which_L)

% ATLAS_L3a_proc_ATBD(D2, params, which_L)
% Inputs: 
% D2: a 2-element structure array containing two beams worth of L2b
%    ICESat-2 data
%    must contain, at minimum, fields:
%        x_RGT: along-track coordinate
%        h: surface height
%        detected: whether each photon was detected (set to all true if no
%        detector model used
%         ...probably others...
% params: a 2-element (one per beam structure containing metadata 
%    about the D2 data and the L2->L3 processing
%      must contain, at minimum, fields:
%        WF: A structure array describing the transmit WF, with fields:
%            t: photon-bin time
%            p: mean power in the bin
% which_L : An optional parameter specifying the segment to be processed.
%    If specified, only the segment centered on which_L will be processed.
%    This is useful for plotting and debugging, can be ignored for batch 
%    processing.



for k=1:length(params); 
    params(k).sigma_pulse=diff(wf_percentile(params(k).WF.t, params(k).WF.p, [0.16 .84]))/2;
end

%allow an effectively unlimited number of iterations
max_ground_bin_iterations=100;

if isfield(D2,'elev');
    for k=1:length(D2);
        D2(k).h=D2(k).elev;
    end
end

% NEW 10/2015: Throw away PE around gaps in the DEM
for kB=1:2;
    D2(kB)=index_struct(D2(kB), isfinite(D2(kB).h) & isfinite(D2(kB).z0)); 
end


% use Anita's ground-finding algo?
if isfield(params(1),'ATL03_sig_find') && params(1).ATL03_sig_find 
    % NEW 10/2015: preserve a copy of the unfiltered data to use with the
    % fallback ground-finder
    D2_unfilt=D2;
    for k=1:2
        D2(k)=index_struct(D2(k),D2(k).ph_class>0);        
    end
end

%NEW 10/2015:removed the veg file option

% NEW 10/2015: exit if x_RGT is of zero length
x_AT_range=range(real(cat(1, D2.x_RGT)));
if isempty(x_AT_range); 
    D3=[];
    return
end

L0=round_to(x_AT_range(1), 20):20:round_to(x_AT_range(2), 20);

if exist('which_L','var');
    L0=L0(find(abs(L0-which_L)==min(abs(L0-which_L)), 1, 'first'));
end

if isfield(params(1),'skip_fpb_corr') && params(1).skip_fpb_corr
    skip_fpb_corr=true;
else
    skip_fpb_corr=false;
end

% NEW 10/2015:  Added signal_selection_status and ATL06_status, 
% added x_PS_ctr, y_PS_ctr, x_RGT, y_RGT, SNR, exit_iteration
D3_fields={ 'signal_selection_status', 'ATL06_status', 'h_expected_rms', 'sigma_geo_AT', 'sigma_geo_XT', ...
    'RGT','GT','PT','cycle','orbit_number','seg_count', ...
    'track','beam', 'BGR','h_initial','W_surface_window_final', ...
    'dh_fit_dx', 'dh_fit_dy', 'h_robust_spread', 'h_rms','h_mean','sigma_h_fit','h_med', ...
    'sigma_dh_fit_dx', 'fpb_error' ,'fpb_med_corr','fpb_mean_corr', 'med_r_fit', ...
    'N_initial', 'N_noise', 'n_fit_photons','fpb_N_corr', 'sigma_photon_est', ...
    'TX_med_corr','TX_mean_corr', 'h_LI', 'z0','m0', 'x_RGT','y_RGT', ...
    'x_PS_ctr', 'y_PS_ctr', 'lat_ctr','lon_ctr', ...
    'fpb_med_corr_sigma','fpb_mean_corr_sigma', ...
    'first_seg_pulse','N_seg_pulses','time','SNR','exit_iteration'};

for kf=1:length(D3_fields);
    D3_empty.(D3_fields{kf})=NaN;
end
% NEW 10/2015: initialize ATL06 status and signal_selection_status
D3_empty.signal_selection_status=uint8(0);
D3_empty.ATL06_status=uint8(0);
D3=repmat(D3_empty, [length(L0),2]);
tic
for k0=1:length(L0);
    % initial processing: find the intial vertical bin centers
    for kB=1:2;
        D3(k0, kB).x_RGT= L0(k0);
        els=(abs(real(D2(kB).x_RGT)-L0(k0))<20) & isfinite(D2(kB).zground) & isfinite(D2(kB).z0);
       
        if isfield(params(kB),'skip') &&  params(kB).skip;
            els=[];
        end
        
        % NEW 10/2015: Eetenive modifications 
        initial_LS_fit_options= struct('initial',true,'Nsigma',3, 'Hwin_min', 10);
        
        if any(els)
            D2sub_all(kB)=index_struct(D2(kB), els);
            [initial_fit_els, D3(k0, kB).signal_selection_status]=choose_ground_strategy(D2sub_all(kB));
            if D3(k0, kB).signal_selection_status <1  % enough confident PE found to determine slope and height of the bin
                initial_LS_fit_options.Hwin_min=4;
            end
        else
            D3(k0, kB).signal_selection_status=bitor(D3(k0, kB).signal_selection_status, 31);
        end
        if D3(k0, kB).signal_selection_status >=4  % not enough confident or padded PE have been found, fall back to alternate strategies
            if D3(k0, kB).signal_selection_status < 16 % at least some detected PE found, select in a window around them
                h0=mean(D2sub_all(kB).h);
                W0=max(diff(range(D2sub_all(kB).h)),  initial_LS_fit_options.Hwin_min);
                els=(abs(real(D2_unfilt(kB).x_RGT)-L0(k0))<20) & isfinite(D2_unfilt(kB).zground) & isfinite(D2_unfilt(kB).z0) & abs(D2_unfilt(kB).h-h0) < W0/2;
                D2sub_all(kB)=index_struct(D2_unfilt(kB), els);
                initial_fit_els=true(size(D2sub_all(kB).h));
                if bitand(D3(k0, kB).signal_selection_status, 8)==8  % not good along-track spread in flagged PE
                    initial_LS_fit_options.NoSlope=true;
                end
            else
                D2sub_all(kB)=find_best_ground_from_hist(D2_unfilt(kB), L0(k0), 20, 10);
                initial_LS_fit_options.NoSlope=true;
                initial_fit_els=true(size(D2sub_all(kB).h));
            end
        end
        
        if ~skip_fpb_corr
            % NEW 10/2015: identify initial fit elements
            initial_fit_els=initial_fit_els(D2sub_all(kB).detected==1);
            D2sub(kB)=index_struct(D2sub_all(kB), D2sub_all(kB).detected==1);
        else
            D2sub(kB)=D2sub_all(kB);
        end
        
        % NEW 10/2015: change in initial fit els calculation (deletion
        % here)
        
        % NEW 10/2015: change in way LS_fit options are passed
        D3(k0, kB)=ATLAS_LS_fit(index_struct(D2sub(kB), initial_fit_els), L0(k0), [0, 0], 1, params(kB), D3(k0, kB), initial_LS_fit_options); % get initial h_med
        D3(k0, kB).h_initial=D3(k0, kB).h_med;
        if isfield(D2sub,'y_RGT');
            ybar(kB)=median(D2sub(kB).y_RGT);
        else
            ybar(kB)=median(imag(D2sub(kB).x_RGT));
        end
    end
    
    dh_fit_dy=diff([D3(k0, :).h_initial])./diff(ybar);
    if ~isfinite(dh_fit_dy)
        dhdy = 0;
    end
     
    for kB=1:2
        % second-round fit: iterate LS fit to convergence
        if kB==1
            h_and_y_other_beam=[D3(k0, 2).h_med, ybar(2)];
        else
            h_and_y_other_beam=[D3(k0, 1).h_med, ybar(1)];
        end
        % NEW 10/2015: set options for final LS fit
        LS_fit_options=struct('initial',false,'Nsigma', 3, 'Hwin_min', 4);
        % try to improve the fit unless the signal_selection_status is
        % equal to zero
        if D3(k0, kB).signal_selection_status ==0      
            LS_fit_options.m_initial=[D3(k0, kB).h_mean; D3(k0, kB).dh_fit_dx];
            max_ground_bin_iterations=1;
        end
                
        [D3(k0, kB), r, els]=ATLAS_LS_fit(D2sub(kB), L0(k0), h_and_y_other_beam, max_ground_bin_iterations, params(kB), D3(k0, kB), LS_fit_options);
   
        D3(k0, kB).med_r_fit=median(r);
        % NEW 10/2015:  Changes in the way the signal selection flag is
        % implemented
        % check PE distribution
        if   isempty(D2sub(kB).x_RGT) || diff(range(D2sub(kB).x_RGT))<20
            D3(k0, kB).h_LI=median(D2sub(kB).h);
            D3(k0, kB).x_RGT=L0(k0);
            D3(k0, kB).signal_selection_status=bitor(D3(k0, kB).signal_selection_status, 64);
         end
        
        if sum(els)<10; 
            D3(k0, kB).signal_selection_status=bitor(D3(k0, kB).signal_selection_status, 32);    
            continue;
        end
        D3(k0, kB).n_fit_photons=length(r);
        [D3(k0, kB).fpb_med_corr, D3(k0, kB).fpb_mean_corr, D3(k0, kB).fpb_N_corr, t_WF, N_WF, N_WF_corr, D3(k0, kB).fpb_med_corr_sigma, D3(k0, kB).fpb_mean_corr_sigma, minGain]=...
            fpb_corr_ATBD(r, D2sub(kB).channel(els), D2sub(kB).pulse_num(els), params(kB).N_det, 57, params(kB).t_dead, skip_fpb_corr, 100e-12);
        % NEW 10/2015: update the ATL06_status flag
        % check if gain correction is valid
        if minGain < 1/(2*params(kB).N_det); 
            D3(k0, kB).ATL06_status=bitor(D3(k0, kB).ATL06_status, 64); 
        end
        D3(k0, kB).N_noise=median(D2sub(kB).BGR)*D3(k0, kB).W_surface_window_final/1.5e8*57;
        sigma_hat_robust=robust_peak_width_from_hist(t_WF, N_WF_corr, D3(k0, kB).N_noise, D3(k0, kB).W_surface_window_final*[-0.5 0.5]/1.5e8);
        
        [D3(k0, kB).TX_med_corr, D3(k0, kB).TX_mean_corr]=correct_for_TX_shape(sigma_hat_robust,[],  params(kB).WF.t, params(kB).WF.p, D3(k0, kB).W_surface_window_final/(1.5e8));
        D3(k0, kB).h_LI=D3(k0, kB).h_mean + D3(k0, kB).fpb_med_corr + D3(k0, kB).TX_med_corr;
        
        % extra code to estimate the true fpb error
        GG=[ones(size(D2sub_all(kB).x0)), real(D2sub_all(kB).x_RGT)-L0(k0)];
        els_all=abs(D2sub_all(kB).h-GG*[D3(k0, kB).h_mean; D3(k0, kB).dh_fit_dx]) < D3(k0, kB).W_surface_window_final/2;
        m_all=GG(els_all,:)\D2sub_all(kB).h(els_all);
        %fprintf(1, 'widths are: [sigma_hat_robust, sigma_hat_robust_all]=\t[%4.2d %4.2d], slope spreading is %3.2d\n', sigma_hat_robust*1.5e8, sigma_hat_robust_all, std(D2sub(kB).zground(D2sub(kB).SigNoise==1)-D2sub(kB).z0(D2sub(kB).SigNoise==1)));
        D3(k0, kB).fpb_error=D3(k0, kB).h_mean-m_all(1);
        
        % now report the across-track coordinates for the segment
         % NEW 10/2015:  Changed order of operations (x_RGT is calculated
         % earlier
        if isfield(D2sub,'y_RGT');
            D3(k0, kB).y_RGT=median(D2sub(kB).y_RGT);
        else
            D3(k0, kB).y_RGT=median(imag(D2sub(kB).x_RGT));
        end
        temp=unique(D2sub(kB).beam(isfinite(D2sub(kB).beam)));
        D3(k0, kB).beam=temp(1);
        temp=unique(D2sub(kB).track(isfinite(D2sub(kB).track)));
        D3(k0, kB).track=temp(1);
        D3(k0, kB).first_seg_pulse=min(D2sub(kB).pulse_num);
        D3(k0, kB).N_seg_pulses=diff(range(D2sub(kB).pulse_num))+1;
        D3(k0, kB).time=median(D2sub(kB).time);
        
        % NEW 10/2015:  Use regression to calculate reference center
        % parameters
        % regress WRT along-track distance to get parameters for segment
        % center
        %'x_PS_ctr', 'y_PS_ctr', 'lat_ctr','lon_ctr'
        Ginv_AT=[ones(size(D2sub(kB).x_RGT(els))) D2sub(kB).x_RGT(els)-L0(k0)]\eye(sum(els));
        Ginv_AT=Ginv_AT(1,:);
        D3(k0, kB).x_PS_ctr=Ginv_AT*D2sub(kB).x0(els);
        D3(k0, kB).y_PS_ctr=Ginv_AT*D2sub(kB).y0(els);
        D3(k0, kB).lat_ctr=Ginv_AT*D2sub(kB).lat(els);
        D3(k0, kB).lon_ctr=Ginv_AT*D2sub(kB).lon(els);
        
    end
    
    dh_fit_dy=(D3(k0,2).h_LI-D3(k0, 1).h_LI)/(diff(ybar));
    if ~isfinite(dh_fit_dy); dh_fit_dy=0; end
    for kB=1:2;
        D3(k0, kB).dh_fit_dy=dh_fit_dy;         
        N_sig=max(0,D3(k0, kB).n_fit_photons-D3(k0, kB).N_noise);
        sigma_signal=sqrt((D3(k0, kB).dh_fit_dy.^2+D3(k0, kB).dh_fit_dx.^2)*params(kB).sigma_x.^2 + (params(kB).sigma_pulse*1.5e8).^2);
        D3(k0, kB).sigma_photon_est=sqrt((D3(k0, kB).N_noise*(D3(k0, kB).W_surface_window_final*.287).^2+N_sig*sigma_signal.^2)/D3(k0, kB).n_fit_photons);
        sigma_per_photon=max(D3(k0, kB).sigma_photon_est, D3(k0, kB).h_robust_spread);
        D3(k0, kB).sigma_h_fit=D3(k0, kB).sigma_h_fit*sigma_per_photon;
        D3(k0, kB).sigma_dh_fit_dx=D3(k0, kB).sigma_dh_fit_dx*sigma_per_photon;     
        D3(k0, kB).sigma_h_LI=max(D3(k0, kB).sigma_h_fit, D3(k0, kB).fpb_med_corr_sigma);
        
        % add in dummy parameter values
        if isfield(params,'RGT');
            D3(k0, kB).RGT=params(kB).RGT;
            D3(k0, kB).GT=params(kB).GT;
            D3(k0, kB).PT=params(kB).PT;
            D3(k0, kB).cycle=params(kB).cycle;
            D3(k0, kB).orbit_number=params(kB).orbit_number;
        end
        
        
    end
    if mod(k0, 250)==0; 
        disp([num2str(k0) ' out of ' num2str(length(L0))]);
        clear D3a;
        f=fieldnames(D3); for kF=1:length(f); for kB=1:2; D3a.(f{kF})(:,kB)=cat(1, D3(:, kB).(f{kF})); end; end
        r=[D3a.h_LI-D3a.z0];fprintf(1, 'mean r =%3.2d, std r =%3.2d\n', [mean(r(isfinite(r))), std(r(isfinite(r)))])
%         if exist('save_file','var');
%             save(save_file, 'D3');
%         end
    end
end
     
%---------------------------------------------------------------------------------------
function [D3, r0, els]=ATLAS_LS_fit(D2, L0, H_and_y_other_beam, N_it, params, D3, options)


if ~exist('options','var');
    options=[];
end

c2=3e8/2;
if isfield(D2,'elev');
    D2.h=D2.elev;
end

if ~options.initial
    XT_slope_est=abs(mean(D2.h)-H_and_y_other_beam(1))/(mean(imag(D2.x_RGT))-H_and_y_other_beam(2));
else
    XT_slope_est=0;
end

m=[]; sub1=[]; r0=[];

% fit an along-track polynomial
%s_ctr=mean(range(real(D2.x_LC)));
ds=real(D2.x_RGT)-L0;
els=true(size(D2.x_RGT));
if isfield(options,'NoSlope') && options.NoSlope;
    G=ones(size(ds(:)));
else
    G=[ones(size(ds(:))), ds(:)];
end

% NEW 10/2015: More initial seg model if it is defined
% if the initial model is defined, use it to initialize the fit
if isfield(options,'m_initial')
    r=D2.h-G*options.m_initial;
    sigma_expected=sqrt((c2*params.sigma_pulse).^2+params.sigma_x.^2*(options.m_initial(2).^2+XT_slope_est.^2)); 
    H_win=max([2*sigma_expected*options.Nsigma, options.Hwin_min]);
    els=abs(r) < H_win/2;
end

D3.N_initial=sum(els);

% iterate to reduce residuals
Noise_Ph_per_m=length(unique(D2.pulse_num))*median(D2.BGR)/c2;
H_win=diff(range(D2.h));
%peak_stats_bin_size=5/cum_noise_rate;
for k=1:N_it;
    if sum(els) < 10 || max(real(D2.x_RGT))-min(real(D2.x_RGT))< 20;
        els=false(size(els)); r=NaN(size(els));
        return
    end
    m=G(els,:)\D2.h(els);
    if ~options.initial
        XT_slope_est=abs(m(1)-H_and_y_other_beam(1))/(mean(imag(D2.x_RGT))-H_and_y_other_beam(2));
        if ~isfinite(XT_slope_est);
            XT_slope_est=0;
        end
    else
        XT_slope_est=0;
    end
    
    r_all=D2.h-G*m;
    r0=r_all(els);
    sigma_r=robust_peak_width_CDF(r0, Noise_Ph_per_m*H_win, [0 H_win]); 
    % NEW 10/2015:  Add the no-slope option
    if isfield(options,'NoSlope') && options.NoSlope; m(2)=0;end
    sigma_expected=sqrt((c2*params.sigma_pulse).^2+params.sigma_x.^2*(m(2).^2+XT_slope_est.^2));
 
    els_last=els;
    H_win=max([2*[sigma_expected, sigma_r]*options.Nsigma, options.Hwin_min]);

    els=abs(r_all) < H_win/2;
    if sum(abs(els_last-els))==0 || sum(els) < 10
        if H_win > options.Hwin_min && options.initial
            H_win=options.Hwin_min;
            els=abs(r_all) < H_win/2;
        else
            break
        end
    end
end

% plotting command:
%clf; plot(real(D2.x0), D2.h,'.'); hold on; plot(real(D2.x0(els)), D2.h(els),'ro'), plot(real(D2.x0(els)), G(els,:)*m,'g')

D3.h_mean=m(1);
D3.dh_fit_dx=m(2);
D3.BGR=median(D2.BGR);
r0=r_all(els);
D3.W_surface_window_final=H_win;

D3.h_robust_spread=iqr(r0)/2;
D3.h_rms=std(r0);
D3.h_med=m(1)+median(r0);

% NEW 10/2015: Calculate the noise estimate and the SNR
Noise_est=H_win*Noise_Ph_per_m-length(unique(D2.pulse_num))*median(D2.BGR)/c2;
D3.SNR=(sum(els)-Noise_est)./Noise_est;

% NEW 10/2015: Report the exit iteration
D3.exit_iteration=k;
if N_it==1;
    D3.exit_iteration=0;
end

if options.initial
    [D3.sigma_h_fit, D3.sigma_dh_fit_dx, D3.z0]=deal(NaN);
    return
end

G1=G(els,:);
Ginv=(G1'*G1)\G1';
CovMat=Ginv*Ginv';

D3.sigma_h_fit=sqrt(CovMat(1,1));
D3.sigma_dh_fit_dx=sqrt(CovMat(2,2));
 
% true height: sample surface every meter, fit with slope 
ds0=round_to(ds, 1); [~, ds0_ind]=unique(ds0);
G=[ones(length(ds0_ind),1), ds(ds0_ind)];
m=G\D2.z0(ds0_ind); 
D3.z0=m(1); D3.m0=m(2);

% NEW 10/2015: choose_ground_strategy scheme
%---------------------
function [els, signal_selection_status]=choose_ground_strategy(D2)

signal_selection_status=uint8(0);
% check if there are enough flagged (not padded) photons
els=D2.ph_class > 1 &  D2.detected & isfinite(D2.h) & isfinite(D2.z0);   
if sum(els) <10
    % not enough flagged PE
    signal_selection_status=bitor(signal_selection_status, 1);
end
if ~any(els) || diff(range(D2.x_RGT(els))) < 20
    % not enough space between first and last PE
    signal_selection_status=bitor(signal_selection_status, 2);
end
% if all good, return
if signal_selection_status==0
     return
end

% now check if there are enough PE including the pad
els==D2.ph_class >= 1 &  D2.detected & isfinite(D2.h) & isfinite(D2.z0);
if sum(els) <10
    % number of detected PE too small
    signal_selection_status=bitor(signal_selection_status, 4);
end
if  ~any(els) || diff(range(D2.x_RGT(els))) < 20  
    % along-track spread of PE too small
    signal_selection_status=bitor(signal_selection_status, 8);
end
if signal_selection_status < 4
    % good spread and good count -> return
    return
end
if sum(els)==0
    % if no flagged PE present, set the 4th bit.
    signal_selection_status=bitor(signal_selection_status, 16);
end


%-------------------------------------------
function D2=find_best_ground_from_hist(D2, L0, W, dz)

these=abs(D2.x_RGT-L0) < 2*W;
if sum(these)<10
    D2=index_struct(D2, []);
    return
end
D2=index_struct(D2, these);
bins=(floor(min(D2.h))+0.5):ceil(max(D2.h));
count=my_histc(D2.h, bins);
C1=conv_corrected(count(:), ones(dz, 1)');
z0=mean(bins(C1==max(C1)));
D2=index_struct(D2, abs(D2.x_RGT-L0) < W & abs(D2.h-z0)<dz & isfinite(D2.h) & isfinite(D2.z0));



