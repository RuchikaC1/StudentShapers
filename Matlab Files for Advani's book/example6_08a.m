%%%%%%%%%% Matlab code for Example 6.8 (a) %%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc

r0    = 0.010;              % radius of central gate, [m]
R     = 0.100;              % outer radius of disk, [m]
h     = 0.001;              % half of the thickness, [m]
MU    = [100 500 1000];     % viscosity of the resin, [Pa.s]
P0    = [5:0.2:10] * 1.0e5; % gage inlet pressure, [Pa]

for i = 1:length(MU)
mu    = MU(i);
rs    = r0/R;
Ps    = P0/(3*mu*R^2/h^2);

tfill = 1./(2*Ps) * (-log(rs)-(1-rs^2)/2);  % time to fill the mold

subplot(3,1,i), plot(P0/100000,tfill,'k-o','linewidth',2)
str = ['viscosity = ',num2str(mu),' Pa.s'];
grid on, legend(str,0)
end
xlabel('Gauge inlet pressure, P_0 [bar]')
ylabel('Time to fill the mold, t_{fill}, [s]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
