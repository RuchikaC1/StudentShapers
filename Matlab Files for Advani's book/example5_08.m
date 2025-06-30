%%%%%%%%%% Matlab code for Example 5.8 %%%%%%%%
clear all, close all, clc
t    = 0:5:30;
P    = [0 11000 22400 33300 42600 53450 60850];
C    = polyfit(t,P,1);
dPdt = C(1)
T    = linspace(t(1),t(end),1000);
Pfit = polyval(C,T);
plot(t,P,'bo', T,Pfit,'k-'), grid on
xlabel('Time, t [s]')
ylabel('Injection Pressure, P_{inj}(t) [Pa]')
legend('Exp. data','1st order curve fit',0)

mu   = 0.24;
Cw   = 0.200;
Ch   = 0.005;
Vf   = 0.474;
Q    = 4e-6;
Kxx  = (mu/(1-Vf))*(Q/(Cw*Ch))^2*(1/dPdt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%