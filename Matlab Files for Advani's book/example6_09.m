%%%%%%%%%% Matlab code for Example 6.9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc

r0    = 0.010;              % radius of central gate, [m]
R     = 0.100;              % outer radius of disk, [m]
h     = 0.001;              % half of the thickness, [m]
mu    = 100;                % viscosity of the resin, [Pa.s]
P0    = 5 * 1.0e5;          % gage inlet pressure, [Pa]

epsil = 0.001;

Rf(1) = r0 * (1 + epsil);
t(1)  = 0;
Q(1)  = 0;
i     = 1;

while Rf < R
        if     Rf(end) < 1.1*r0, dt = 0.00001;
        elseif Rf(end) < 1.5*r0, dt = 0.0001;
        elseif Rf(end) < 2.0*r0, dt = 0.001;
        elseif Rf(end) < 5.0*r0, dt = 0.01;
        end
f       = (P0*h^2)/(3*mu) * (1/(Rf(i)*log(Rf(i)/r0)));
Q(i)    = 4*pi*h*Rf(i)*f; 
Rf(i+1) = Rf(i) + f*dt;
t(i+1)  = t(i) + dt;
i       = i+1;    
end
Q(i)    = 4*pi*h*Rf(i)*f; 

figure(1)
subplot(2,1,1), plot(t,Rf*1000,'k-','linewidth',2)
grid on
xlabel('Time, t [s]')
ylabel('Radial flow front position, R_f [mm]')

subplot(2,1,2), plot(t(2:end),Q(2:end)*1e6,'k-','linewidth',2)
grid on
xlabel('Time, t [s]')
ylabel('Flow rate, Q [cc/s]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
