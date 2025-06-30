%%%%%%%%%%%%%%% Matlab code for Example 8.5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
Vf        = 0.40;        % average fiber volume fraction
kf        = 1.20;        % thermal conductivity of fiber,  [W/(m K)]
km        = 0.20;        % thermal conductivity of matrix, [W/(m K)]
rhof      = 2500;        % density of fiber,  [kg/m3]
rhom      = 1200;        % density of matrix, [kg/m3]
cf        =  800;        % specific heat capacity of fiber,  [J/(kg K)]
cm        =  700;        % specific heat capacity of matrix, [J/(kg K)]
Tin       =  20 + 273.2; % initial temperature, [K]
Tw        = 120 + 273.2; % mold wall temperature,  [K]
h         = 0.020;       % thickness of the part, [m]
R         = 8.3145;      % constant to be used in cure kinetics
CA0       = 0.5;         % constant to be used in cure kinetics
k0        = 2;           % constant to be used in cure kinetics
E         = 1000 * R;    % constant to be used in cure kinetics
LAMBDA    = 1;           % constant to be used in cure kinetics
Hr        = 280000;      % constant to be used in cure kinetics, [J/kg]
k         =  1/((1-Vf)/km + Vf/kf);           % thermal conductivity of composite
rho_c     = Vf*(rhof*cf) + (1-Vf)*(rhom*cm);  % (density * specific heat) of composite
Nz        = 21;                               % # of nodes in z
z         = linspace(0,h,Nz);                 % z' array, [m]
dz        = z(2)-z(1);                        % Delta z', [m]

T         = Tin * ones(size(z));
T(1)      = Tw;
T(end)    = Tw;
C         = zeros(size(z));
To(1)     = T(1                   );
Tq(1)     = T((round(Nz-1)*0.25)+1);
Tm(1)     = T((round(Nz-1)*0.50)+1);
Co(1)     = C(1                   );
Cq(1)     = C((round(Nz-1)*0.25)+1);
Cm(1)     = C((round(Nz-1)*0.50)+1);
t(1)      = 0;
dt        = 0.01;
j         = 1;
Cnext     = C;
Tnext     = T;

while min(C) < 0.85
for i = 2:Nz-1    
dCdt      = k0 * exp(-E/(R*T(i))) * CA0 * (1-C(i)) * (LAMBDA-C(i));
dTdt      = (1/rho_c)*(k*((T(i+1)-2*T(i)+T(i-1))/(dz^2))+(1-Vf)*rhom*Hr*dCdt);
Cnext(i)  = C(i) + dCdt*dt;
Tnext(i)  = T(i) + dTdt*dt;
end
Cs        = C(end-1) + (C(end-1)-C(end-2));
Cnext(1)  = Cs;
Cnext(end)= Cs;
C         = Cnext;
T         = Tnext;
t(j+1)    = t(j) + dt;
j         = j + 1;
To(j)     = T(1                   );
Tq(j)     = T((round(Nz-1)*0.25)+1);
Tm(j)     = T((round(Nz-1)*0.50)+1);
Co(j)     = C(1                   );
Cq(j)     = C((round(Nz-1)*0.25)+1);
Cm(j)     = C((round(Nz-1)*0.50)+1);
    if abs(mod(t(end),10) - 10) < 1.0e-8 |  mod(t(end),10) < 1.0e-8 
    figure(1)
    subplot(121), plot(z*1000,T-273.2,'k-','linewidth',1), grid on, hold on
    xlabel('z [mm]'), ylabel('Temperature, T [^oC]')
    subplot(122), plot(z*1000,C,'k-','linewidth',1), grid on, hold on
    xlabel('z [mm]'), ylabel('Extent of reaction, C^*')
    pause(0.1)
    else
    end
end

figure(2)
subplot(121), plot(t,Tq-273.2,'k-',t,Tm-273.2,'k--','linewidth',1)
grid on, legend('z = h/4','z = h/2',0)
xlabel('Time, t [s]'), ylabel('Temperature, T [^oC]')
subplot(122), plot(t,Cq,'k-',t,Cm,'k--','linewidth',1), grid on, hold on
xlabel('Time, t [s]'), ylabel('Extent of reaction, C^*')
grid on, legend('z = h/4','z = h/2',0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%