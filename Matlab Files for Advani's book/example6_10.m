%%%%%%%%%% Matlab code for Example 6.10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc

T1    = 280 + 273.2;        % inlet temperature, [oK]
mu    = 50000;              % viscosity of the resin, [Pa.s]
h     = 0.005;              % thickness of the channel, [m]
W     = 0.100;              % width of the channel, [m]
L     = 0.5;                % length of the channel, [m]
k     = 0.8;                % thermal conductivity of composite, [W/(m oK)]
Pi    = [1:2:50] * 1.0e6;   % inlet gauge pressure, [Pa]
Pe    = 0.0e6;              % exit  pressure, [Pa]
DP    = Pe - Pi;

Q     = -DP/(2*mu*L) * ((4*W*h^3)/3);                 % flow rate, [m^3/s]

Tmid  = T1 * ( 1 + ( (3*mu)/(16*k*T1*W^2*h^2)*Q.^2) );

figure(1)
subplot(2,1,1), plot(Pi*1.0e-6,Q*1.0e6,'b-o','linewidth',2)
grid on
xlabel('Inlet gauge pressure, P_i [MPa]')
ylabel('Flow rate, Q [cc/s]')  

subplot(2,1,2), plot(Pi*1.0e-6,Tmid,'k-o','linewidth',2)
grid on
xlabel('Inlet gauge pressure, P_i [MPa]')
ylabel('Mid-plane temperature, T_{mid-plane} [^oK]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
