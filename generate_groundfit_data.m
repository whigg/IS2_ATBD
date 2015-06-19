%-----------------------
% Generate groundfit data
% make a dataset for testing the ground fitting routines.  Given a set of
% parameters

%parameter section: comment or uncomment the next two lines to generate
%output for different parameter sets
%out_file='groundfit_test_D2_slope_0.1_tau_0.mat'; this_tau=0; max_slope_mag=0.1;
out_file='groundfit_test_D2_slope_0.1_tau_-2.mat'; this_tau=-2; max_slope_mag=0.1;

% DEM geometry.  Generate a sine wave in x with a gradually increasing
% cross-track slope in y
Lx=20*3000;
Ly=300;
dx=5;
DEM.x=0:dx:Lx;
DEM.y=-Ly/2:dx:Ly/2;
[xg,yg]=meshgrid(DEM.x, DEM.y);

y_slope=2*xg/Lx*max_slope_mag/sqrt(2);
x_slope_mag=0.2/sqrt(2);
lambda=2000;

DEM.z=y_slope.*yg + x_slope_mag*lambda/2/pi*cos(2*pi*xg/lambda);

load WF_est
params_L=struct('N_per_pulse', 12, 't_dead', 3.2e-9, 'sigma_x', 7.2,'sigma_pulse', 1.6e-9,'c', 3e8, 'N_det', 16, 'NoiseRate', 1e7,'H_window', 39, 'WF', WF,'DEBUG', false);
params_R=params_L; params_R.N_per_pulse=3; params_R.N_det=4;
params_L.refine_ground_bin_threshold=50;
params_R.refine_ground_bin_threshold=50;

x0=min(DEM.x)+34 : 0.7 : max(DEM.x)-34;
% set the atmospheric transmittance
ATM_xmit=zeros(size(x0))+exp(this_tau);

% generate the photon-elevation data.  Note that geographic locations are
% expressed in complex coordinates
clear D2;
D2(1)=det_sim(DEM, x0-1i*45, params_R, ATM_xmit);
D2(2)=det_sim(DEM, x0+1i*45, params_L, ATM_xmit);

D2(1).x_RPT=real(D2(1).x0)-1i*45;
D2(2).x_RPT=real(D2(2).x0)+1i*45;
for k=1:2
    D2(k).beam=ones(size(D2(k).h));
    D2(k).track=ones(size(D2(k).h));
    D2(k).time=D2(k).pulse_num*1e-4;
end

params=[params_R params_L];

save(out_file,'D2','DEM','params')




