%%%%%%%% Matlab code for Example 5.9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
mu    = 0.24;
 h    = 0.005;
Vf    = 0.30;
Q     = 4e-6;
R1    = 0.006;
R2    = 0.012;
rmold = 0.40;
t     = [0:5:30]';
P     = [0 125 147 160 169 176 182]' * 1000;
C     = [log(1+(Q*t)./(pi*h*R1^2*(1-Vf)))];
X     = (C'*C)\(C'*P);
K     = ((mu*Q)/(4*pi*h))/X

EmptyVolume1 = (pi*(rmold^2 - R1^2))*h*(1-Vf)
t_fill       = EmptyVolume1/Q
T1           = linspace(0,     t(end),1000);
T2           = linspace(t(end),t_fill,1000);
T3           = linspace(0     ,t_fill,1000);
PP1          = (mu*Q)/(4*pi*h*K) * log(1+(Q*T1)./(pi*h*R1^2*(1-Vf)));
PP2          = (mu*Q)/(4*pi*h*K) * log(1+(Q*T2)./(pi*h*R1^2*(1-Vf)));
PP3          = (mu*Q)/(4*pi*h*K) * log(1+(Q*T3)./(pi*h*R2^2*(1-Vf)));
plot(t,P,'ro', T1,PP1,'r-.', T2,PP2,'b-', T3,PP3,'k-.')
grid, xlabel('Time, t [s]'), ylabel('Injection Pressure, P_{inj}(t) [Pa]')
legend('Exp. data when R = R_1', ...
       'Curve fit when R = R_1', ...
       'Extension of curve fit when R = R_1', ...
       'When R = R_2',0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%