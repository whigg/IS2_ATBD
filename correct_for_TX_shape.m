function [dM, dCtr]=correct_for_TX_shape(t, z, t_TX, TX, HW, Nrate, N_sig)

if ~isempty(z);
    if exist('Nrate','var')
        dt=t(2)-t(1);
        z1=max(0, z-Nrate*dt);
        sigma_rx=(wf_percentile(t, z1, 0.84)-wf_percentile(t, z1, 0.16))/2;
    else
        sigma_rx=(wf_percentile(t, z, 0.84)-wf_percentile(t, z, 0.16))/2;
    end
else
    sigma_rx=t;
end

sigma_tx=(wf_percentile(t_TX, TX, 0.84)-wf_percentile(t_TX, TX, 0.16))/2;
dt=t_TX(2)-t_TX(1);

% generate a WF that looks like the tx pulse convolved with the pulse-spreading function
% broaden the TX waveform if it's narrower than the RX waveform
if sigma_rx > sigma_tx
    sigma_G=sqrt(sigma_rx.^2-sigma_tx.^2);
    Gsamps=ceil(sigma_G/dt);
    G=gaussian(dt*(-4*Gsamps:4*Gsamps), 0, sigma_G);
    t1=[(-4*Gsamps:-1)*dt+t_TX(1) t_TX(:)' (1:4*Gsamps)*dt+t_TX(end)]; t1=t1(:);
else
    G=1;
    t1=t_TX(:);
end
TX_WF_broadened=conv(TX(:), G(:));
 
% find the 16th and 84th percentiles of the broadened WF, use these to truncate it
LCH=wf_percentile(t1, TX_WF_broadened, [0.16 0.5 0.84]);

if exist('Nrate', 'var')
    N_add_early=ceil((1.25*HW/2-(LCH(2)-t1(1)))/dt);
    if N_add_early > 0;
        % add samples to the front of the WF
        t1=[t1(1)+(-N_add_early:-1)'*dt ; t1];
        TX_WF_broadened=[zeros(N_add_early,1); TX_WF_broadened];
    end
    N_add_late=ceil((1.25*HW/2-(t1(end)-LCH(2)))/dt);
    if N_add_late > 0;
        % add samples to the front of the WF
        t1=[t1; t1(end)+(1:N_add_late)'*dt ];
        TX_WF_broadened=[TX_WF_broadened; zeros(N_add_late,1)  ];
    end
       
    TX_WF_broadened=N_sig*TX_WF_broadened/sum(TX_WF_broadened)+Nrate*dt;

    % iterate to find the centroid
    last_t_ctr=t1(1);
    t_ctr=LCH(2);
    while abs(t_ctr-last_t_ctr) > 1e-3/3e8;
        els=abs(t1-t_ctr)<HW/2;
        last_t_ctr=t_ctr;
        t_ctr=sum(TX_WF_broadened(els).*t1(els))/sum(TX_WF_broadened(els));
    end
    TXBm0=t_ctr;
else
    % truncate around the median
    TXB_els=abs(t1-LCH(2)) < HW/2;%1.5*(LCH(3)-LCH(1));
    
    %find the mean of the RX waveform (edited by twice the window)
    %TXBm0=wf_percentile(t1(TXB_els), TX_WF_broadened(TXB_els), 0.5);
    TXBm0=sum(t1(TXB_els).*TX_WF_broadened(TXB_els))./sum(TX_WF_broadened(TXB_els));
end

%truncate the synthetic WF around its edited mean by HW
TXB_els_HW=abs(t1-TXBm0) < HW/2;
TXBm=wf_percentile(t1(TXB_els_HW), TX_WF_broadened(TXB_els_HW), 0.5);
 
%dM=TXm-RXm;
TimeToH=1.5e8;
dM=TXBm*TimeToH;
dCtr=sum(t1(TXB_els_HW).*TX_WF_broadened(TXB_els_HW))./sum(TX_WF_broadened(TXB_els_HW))*TimeToH;
 