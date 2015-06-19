function [D3]=ATLAS_L3a_proc_ATBD(D2, params, which_L)

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
  
% use Anita's ground-finding algo?
if isfield(params(1),'ATL03_sig_find') && params(1).ATL03_sig_find 
    for k=1:2
        D2(k)=index_struct(D2(k),D2(k).ph_class>0);
    end
end


for k=1:2
    if isfield(params(k),'veg_ground_file') && exist(params(k).veg_ground_file,'file');
        veg_gf=load(params(k).veg_ground_file);
        D2(k)=index_struct(D2(k), veg_gf.filter_idx);
    end
end
   
x_AT_range=range(real(cat(1, D2.x_RPT)));

L0=round_to(x_AT_range(1), 20):20:round_to(x_AT_range(2), 20);

if exist('which_L','var');
    L0=L0(find(abs(L0-which_L)==min(abs(L0-which_L)), 1, 'first'));
end

if isfield(params(1),'skip_fpb_corr') && params(1).skip_fpb_corr
    skip_fpb_corr=true;
else
    skip_fpb_corr=false;
end

D3_fields={'track','beam', 'BGR','h_initial','W_surface_window_final', ...
    'dh_fit_dx', 'dh_fit_dy', 'h_robust_spread', 'h_rms','h_mean','sigma_h_fit','h_med', ...
    'sigma_dh_fit_dx', 'fpb_error' ,'fpb_med_corr','fpb_mean_corr', ...
    'N_initial', 'N_noise', 'N_window','fpb_N_corr', 'sigma_photon_est', ...
    'TX_med_corr','TX_mean_corr', 'h_LI', 'z0','m0', 'x_RPT','y_RPT', ...
    'fpb_med_corr_sigma','fpb_mean_corr_sigma', ...
    'first_seg_pulse','N_seg_pulses','time'};

for kf=1:length(D3_fields);
    D3_empty.(D3_fields{kf})=NaN;
end

D3=repmat(D3_empty, [length(L0),2]);
tic
for k0=1:length(L0);
    % initial processing: find the intial bin centers
    for kB=1:2;      
        els=(abs(real(D2(kB).x_RPT)-L0(k0))<20) & isfinite(D2(kB).zground) & isfinite(D2(kB).z0);
        if isfield(params(kB),'skip') &&  params(kB).skip;
           els=[];
        end
        
        D3(k0, kB).N_initial=sum(els);
        if D3(k0, kB).N_initial==0
            ybar(kB)=NaN;
            continue
        end
        
        
        D2sub_all(kB)=index_struct(D2(kB), els);
        if ~skip_fpb_corr
            D2sub(kB)=index_struct(D2sub_all(kB), D2sub_all(kB).detected==1);
        else
            D2sub(kB)=D2sub_all(kB);
        end
                
        D3(k0, kB)=ATLAS_LS_fit(D2sub(kB), L0(k0), [0, 0], 1, params(kB), D3(k0, kB), struct('initial',true,'Nsigma',3, 'Hwin_min', 10)); % get initial h_med
        D3(k0, kB).h_initial=D3(k0, kB).h_med;
        if isfield(D2sub,'y_RPT');
            ybar(kB)=median(D2sub(kB).y_RPT);
        else
            ybar(kB)=median(imag(D2sub(kB).x_RPT));
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
         
        [D3(k0, kB), r, els ]=ATLAS_LS_fit(D2sub(kB), L0(k0), h_and_y_other_beam, max_ground_bin_iterations, params(kB), D3(k0, kB), struct('initial',false,'Nsigma',3,'Hwin_min', 4));
  
        if sum(els)<10; continue; end
        D3(k0, kB).N_window=length(r);
        [D3(k0, kB).fpb_med_corr, D3(k0, kB).fpb_mean_corr, D3(k0, kB).fpb_N_corr, t_WF, N_WF, N_WF_corr, D3(k0, kB).fpb_med_corr_sigma, D3(k0, kB).fpb_mean_corr_sigma]=fpb_corr_ATBD(r, D2sub(kB).channel(els), D2sub(kB).pulse_num(els), params(kB).N_det, 57, params(kB).t_dead, skip_fpb_corr, 100e-12);
        D3(k0, kB).N_noise=median(D2sub(kB).BGR)*D3(k0, kB).W_surface_window_final/1.5e8*57;
        sigma_hat_robust=robust_peak_width_from_hist(t_WF, N_WF_corr, D3(k0, kB).N_noise, D3(k0, kB).W_surface_window_final*[-0.5 0.5]/1.5e8);
        
        [D3(k0, kB).TX_med_corr, D3(k0, kB).TX_mean_corr]=correct_for_TX_shape(sigma_hat_robust,[],  params(kB).WF.t, params(kB).WF.p, D3(k0, kB).W_surface_window_final/(1.5e8));
        D3(k0, kB).h_LI=D3(k0, kB).h_mean + D3(k0, kB).fpb_med_corr + D3(k0, kB).TX_med_corr;
        
        % extra code to estimate the true fpb error
        GG=[ones(size(D2sub_all(kB).x0)), real(D2sub_all(kB).x_RPT)-L0(k0)];
        els_all=abs(D2sub_all(kB).h-GG*[D3(k0, kB).h_mean; D3(k0, kB).dh_fit_dx]) < D3(k0, kB).W_surface_window_final/2;
        m_all=GG(els_all,:)\D2sub_all(kB).h(els_all);
        %fprintf(1, 'widths are: [sigma_hat_robust, sigma_hat_robust_all]=\t[%4.2d %4.2d], slope spreading is %3.2d\n', sigma_hat_robust*1.5e8, sigma_hat_robust_all, std(D2sub(kB).zground(D2sub(kB).SigNoise==1)-D2sub(kB).z0(D2sub(kB).SigNoise==1)));
        D3(k0, kB).fpb_error=D3(k0, kB).h_mean-m_all(1);
        
        % now report the along-track coordinates for the segment
        D3(k0, kB).x_RPT= L0(k0);
        if isfield(D2sub,'y_RPT');
            D3(k0, kB).y_RPT=median(D2sub(kB).y_RPT);
        else
            D3(k0, kB).y_RPT=median(imag(D2sub(kB).x_RPT));
        end
        temp=unique(D2sub(kB).beam(isfinite(D2sub(kB).beam)));
        D3(k0, kB).beam=temp(1);
        temp=unique(D2sub(kB).track(isfinite(D2sub(kB).track)));
        D3(k0, kB).track=temp(1);
        D3(k0, kB).first_seg_pulse=min(D2sub(kB).pulse_num);
        D3(k0, kB).N_seg_pulses=diff(range(D2sub(kB).pulse_num))+1;
        D3(k0, kB).time=median(D2sub(kB).time);
    end
    
    dh_fit_dy=(D3(k0,2).h_LI-D3(k0, 1).h_LI)/(diff(ybar));
    if ~isfinite(dh_fit_dy); dh_fit_dy=0; end
    for kB=1:2;
        D3(k0, kB).dh_fit_dy=dh_fit_dy;         
        N_sig=max(0,D3(k0, kB).N_window-D3(k0, kB).N_noise);
        sigma_signal=sqrt((D3(k0, kB).dh_fit_dy.^2+D3(k0, kB).dh_fit_dx.^2)*params(kB).sigma_x.^2 + (params(kB).sigma_pulse*1.5e8).^2);
        D3(k0, kB).sigma_photon_est=sqrt((D3(k0, kB).N_noise*(D3(k0, kB).W_surface_window_final*.287).^2+N_sig*sigma_signal.^2)/D3(k0, kB).N_window);
        sigma_per_photon=max(D3(k0, kB).sigma_photon_est, D3(k0, kB).h_robust_spread);
        D3(k0, kB).sigma_h_fit=D3(k0, kB).sigma_h_fit*sigma_per_photon;
        D3(k0, kB).sigma_dh_fit_dx=D3(k0, kB).sigma_dh_fit_dx*sigma_per_photon;     
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

c2=3e8/2;
if isfield(D2,'elev');
    D2.h=D2.elev;
end

if ~options.initial
    XT_slope_est=abs(mean(D2.h)-H_and_y_other_beam(1))/(mean(imag(D2.x_RPT))-H_and_y_other_beam(2));
else
    XT_slope_est=0;
end

h_expected_rms0=sqrt((c2*params.sigma_pulse).^2+params.sigma_x.^2*XT_slope_est.^2);
m=[]; sub1=[]; r0=[];

% fit an along-track polynomial
%s_ctr=mean(range(real(D2.x_LC)));
ds=real(D2.x_RPT)-L0;
els=true(size(D2.x_RPT));
G=[ones(size(ds(:))), ds(:)];  
% iterate to reduce residuals

Noise_Ph_per_m=length(unique(D2.pulse_num))*median(D2.BGR)*2/3e8;
H_win=diff(range(D2.h));
%peak_stats_bin_size=5/cum_noise_rate;
for k=1:N_it;
    if sum(els) < 10 || max(real(D2.x_RPT))-min(real(D2.x_RPT))< 20;
        els=false(size(els)); r=NaN(size(els));
        return
    end
    m=G(els,:)\D2.h(els);
    if ~options.initial
        XT_slope_est=abs(m(1)-H_and_y_other_beam(1))/(mean(imag(D2.x_RPT))-H_and_y_other_beam(2));
        if ~isfinite(XT_slope_est);
            XT_slope_est=0;
        end
    else
        XT_slope_est=0;
    end
    
    r_all=D2.h-G*m;
    r0=r_all(els);
    sigma_r=robust_peak_width_CDF(r0, Noise_Ph_per_m*H_win, [0 H_win]);  
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
