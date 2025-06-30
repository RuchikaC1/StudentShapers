%%%%%%%%%%%%%% Matlab code for Example 7.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
L0         = 1;     % original length, [m]

%%% part (a), constant strain rate:
epsdot     = 0.01;   % strain rate, [1/s]
dt         = 0.01;
La(1)      = L0;
ta(1)       = 0;
k          = 1;
while La(end) < 2*L0
k          = k+1;
ta(k)      = ta(k-1) + dt;
La(k)      = L0*exp(epsdot*ta(k));  % length, [m]
EPSDOTa(k) = epsdot;                % strain rate, [1/s]
end

%%% part (b), constant elongation speed:
V          = 0.01;   % elongation speed, [m/s]
dt         = 0.01;
Lb(1)      = L0;
tb(1)      = 0;
k          = 1;
while Lb(end) < 2*L0
k          = k+1;
tb(k)      = tb(k-1) + dt;
Lb(k)      = L0 + V*tb(k);     % length, [m]
EPSDOTb(k) = 1/((L0/V)+tb(k));          % strain rate, [1/s]
end

figure(1)
subplot(2,1,1), plot(ta,La,'b-',tb,Lb,'k--','linewidth',2)
xlabel('Time, t [seconds]'), ylabel('Length, L [m]'), grid on 
legend('d\epsilon / dt = constant = 0.1 1/s','dL / dt = V = constant = 0.1 m/s',0)

subplot(2,1,2), plot(ta,EPSDOTa,'b-',tb,EPSDOTb,'k--','linewidth',2)
xlabel('Time, t [seconds]'), ylabel('Strain rate, d\epsilon / dt, [1/s]')
grid on
legend('d\epsilon / dt = constant = 0.1 1/s','dL / dt = V = constant = 0.1 m/s',0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%