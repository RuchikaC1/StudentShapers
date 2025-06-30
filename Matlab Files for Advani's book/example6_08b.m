%%%%%%%%%% Matlab code for Example 6.8 (b) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc

r0      = 0.010;              % radius of central gate, [m]
R       = 0.100;              % outer radius of disk, [m]
h       = 0.001;              % half of the thickness, [m]
MU      = [100 500 1000];     % viscosity of the resin, [Pa.s]
QQ      = [1 2 5]  * 1.0e-6;  % gage inlet pressure, [Pa]

for i = 1:length(MU)
for j = 1:length(QQ)
mu      = MU(i);
Q       = QQ(j);

clear t Rf P0
t(1)    = 0;
Rf(1)   = r0;
P0(1)   = 0;
k       = 1;
dt      = 0.01;
str_cur = ['k- ';'b-.';'r--'];

while Rf < R
t(k+1)  = t(k) + dt;
Rf(k+1) = sqrt(Q*t(k+1)/(2*pi*h) + r0^2); 
P0(k+1) = (3*Q*mu)/(8*pi*h^3)*log(1+(Q*t(k+1))/(2*pi*h*r0^2));    
k       = k+1;
end

figure(1)
subplot(3,1,i), plot(t,P0/100000,str_cur(j,:),'linewidth',2), hold on
str_title  = ['Viscosity = ',num2str(mu),' Pa.s'];
grid on, title(str_title)
xlabel('Time, t [s]')
ylabel('Inlet pressure, P_{0}(t), [bar]')

figure(2)
subplot(3,1,i), plot(t,Rf*1000,str_cur(j,:),'linewidth',2), hold on
str_title  = ['Viscosity = ',num2str(mu),' Pa.s'];
grid on, title(str_title)
xlabel('Time, t [s]')
ylabel('Radial flow front position, R_{f}(t), [mm]')

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
