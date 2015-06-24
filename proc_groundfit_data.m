function D3=proc_groundfit_data(in_file, DOPLOT)

%proc_groundfit_data
% given an input file (saved by generate_groundfit_data) process the
% ATL03-type data to ATL06-type data and evaluate the results.
% To run the examples from the 'parameters' section
%  D3=proc_groundfit_data('groundfit_test_D2_slope_0.1_tau_-2.mat', true)
%  D3=proc_groundfit_data('groundfit_test_D2_slope_0.1_tau_0.mat', true)

% fields to look at in the output. all are N_segs x 2.
%  --the first column gives weak-beam values, the second gives strong-beam values
% z0:  segment-mean elevation for signal photons for the center of the footprints.  this is 'truth'
% h_LI: land-ice elevation with all corrections applied.
% h_mean: least-squares fit elevation with no corrections applied.  this is
%         the output of the ground-fit section.
% dh_fit_dx : along-track slope
% sigma_h_fit : the propagated error in the least-squares fit
% W_surface_window_final : the total height of the ground-fit window.
% N_window : the number of photons in the refined window

load(in_file,'D2','params','DEM');

[D3a]=ATLAS_L3a_proc_ATBD(D2, params);
f=fieldnames(D3a);
for kF=1:length(f);
    for kB=1:2;
        D3.(f{kF})(:, kB)=cat(1, D3a(:, kB).(f{kF}));
    end
end

if exist('DOPLOT','var') && DOPLOT
    figure
    clf
    
    for k=1:2;
        subplot(2,1,k)
        plot(real(D2(k).x0), D2(k).h,'.'); hold on;
        plot([1; 1]*D3.x_RPT(:,k)'+[-20; 20]*ones(size(D3.x_RPT(:,k)')), [1; 1]*D3.h_LI(:,k)'+[-20; 20]*D3.dh_fit_dx(:,k)','r--','linewidth', 1)
        plot([1; 1]*D3.x_RPT(:,1)'+[-20; 20]*ones(size(D3.x_RPT(:,k)')), [1; 1]*D3.h_LI(:,k)'+[-20; 20]*D3.dh_fit_dx(:,k)'-[1; 1]*D3.W_surface_window_final(:,k)'/2,'g--','linewidth', 1)
        plot([1; 1]*D3.x_RPT(:,1)'+[-20; 20]*ones(size(D3.x_RPT(:,k)')), [1; 1]*D3.h_LI(:,k)'+[-20; 20]*D3.dh_fit_dx(:,k)'+[1; 1]*D3.W_surface_window_final(:,k)'/2,'g--','linewidth', 1)
        good=poisson_p_table(D3.N_noise(:,k), D3.N_window(:,k) ) > 0.95;
        plot([1; 1]*D3.x_RPT(good,k)'+[-20; 20]*ones(size(D3.x_RPT(good,k)')), [1; 1]*D3.h_LI(good,k)'+[-20; 20]*D3.dh_fit_dx(good,k)','r-','linewidth', 2)
        plot([1; 1]*D3.x_RPT(good,1)'+[-20; 20]*ones(size(D3.x_RPT(good,k)')), [1; 1]*D3.h_LI(good,k)'+[-20; 20]*D3.dh_fit_dx(good,k)'-[1; 1]*D3.W_surface_window_final(good,k)'/2,'g-','linewidth', 2)
        plot([1; 1]*D3.x_RPT(good,1)'+[-20; 20]*ones(size(D3.x_RPT(good,k)')), [1; 1]*D3.h_LI(good,k)'+[-20; 20]*D3.dh_fit_dx(good,k)'+[1; 1]*D3.W_surface_window_final(good,k)'/2,'g-','linewidth', 2)       
    end
    
    
    figure; set(gcf,'defaultaxesfontsize', 12);
    sigma_h_LI=max(D3.sigma_h_fit, D3.fpb_med_corr_sigma);
    r=D3.h_LI-D3.z0;
    error_bins=0:.0125:1;
    for k=1:length(error_bins);
        els=abs(sigma_h_LI-error_bins(k))<(error_bins(2)-error_bins(1))/2 & isfinite(r);
        if sum(els) > 6
            rmsR(k)=sqrt(mean((r(els)).^2));
            rmed(k)=median(abs(r(els)))/.667;
            Riqr(k)=iqr(r(els))/2;
        else
            [rmsR(k), rmed(k), Riqr(k)]=deal(NaN);
        end
    end
    plot(error_bins, rmsR, error_bins, rmed, error_bins, Riqr);
    set(findobj(gca,'type','line'),'marker','.','linewidth', 2)
    axis equal tight; set(gca,'ylim', [0 max([rmsR(:); rmed(:); Riqr(:)])*1.5]);
    XR=get(gca,'xlim'); set(gca,'xlim', [0 XR(2)*1.1]);
    hold on;
    plot([0 XR(2)], [0 XR(2)]);
    legend('RMS','MAD','IQR', '1:1');
    xlabel('error estimate, m'); ylabel('residual magnitude, m');
    
    % report biases
    bias_check_fields={'all', 'fpb_med_corr','TX_med_corr'};
    [gx, gy]=gradient(DEM.z, DEM.x, DEM.y);
    slope_vals=abs(interp2(DEM.x, DEM.y, gx+1i*gy, real(D3.x_RPT), imag(D3.x_RPT)));
    slope_ints=[0:0.01:.1];
    for kF=1:length(bias_check_fields);
        if ~strcmp(bias_check_fields{kF},'all');
            fm=D3.h_LI-D3.(bias_check_fields{kF})-D3.z0;
        else
            fm=D3.h_LI-D3.z0;
        end
        for kB=1:2
            for kS=1:length(slope_ints)-1;
                these=slope_vals(:,kB) > slope_ints(kS) & slope_vals(:, kB) < slope_ints(kS+1) & isfinite(fm(:, kB));
                Rminus(kF, kS, kB)=mean(fm(these, kB));
                Sminus(kF, kS, kB)=std(fm(these, kB))/sqrt(sum(these));
            end
        end
    end
    
    fprintf(1, 'biases in mm\n')
    for kB=1:2;
        fprintf(1, '\nBeam %d\n', kB);
        fprintf(1,'pulse width\t%s\t-%s\t-%s\n', bias_check_fields{:});
        for kS=1:length(slope_ints)-1;
            fprintf(1, '%3.0f',  1000*sqrt((7.2*(slope_ints(kS)+slope_ints(kS+1))/2).^2 + (1.6*.15).^2));
            for kF=1:length(bias_check_fields);
                fprintf(1, '\t%3.1f (%3.1f)', Rminus(kF, kS, kB)*1000, Sminus(kF, kS, kB)*1000);
            end
            fprintf(1, '\n');
        end
    end
end
