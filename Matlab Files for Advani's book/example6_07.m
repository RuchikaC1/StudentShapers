%%%%%%%%%% Matlab code for Example 6.7  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
L     = 1.920;              % screw length, [m]
D     = 0.080;              % screw diameter, [m]
Nscrew= 0.20:0.01:0.50;     % screw speed, [rev/s] = [rpm/60]
Ho    = 0.020;              % inlet channel depth, [m]
Hl    = 0.004;              % exit  channel depth, [m]
theta = (20/180)*pi;        % helix angle, [rad]
KG    = 1.2e-8;             % related to die geometry and characteristics, [m^3]
e     = 0.008;              % flight width, [m]

Vf    = 0.20;        % Fiber volume fraction
rhof  = 2500;        % density of fiber,  [kg/m3]
rhom  = 1200;        % density of matrix, [kg/m3]
kf    = 1.20;        % thermal conductivity of fiber,  [W/(m K)]
km    = 0.20;        % thermal conductivity of matrix, [W/(m K)]
cf    =  800;        % specific heat capacity of fiber,  [J/(kg K)]
cm    =  700;        % specific heat capacity of matrix, [J/(kg K)]

a     = 6.5e5;       % where mu = a exp(-b T), [Pa.s]
b     = 0.009;       %                       , [1/oK]

To    = 50 + 273.2;  % inlet temperature, [oK]
Vo    = a*exp(-b*To);% viscosity at the inlet, [Pa.s]

for j = 1:length(Nscrew)
N     = Nscrew(j);
Ls    = tan(theta)*(pi*D);  % screw lead, [m]
W     = Ls*cos(theta) - e;  % channel width, [m]
l     = L/sin(theta);       % spiral channel length, [m]
Vz    = pi*D*N*cos(theta);  % relative velocity in z direction, [m/s]
Vx    = pi*D*N*sin(theta);  % relative velocity in z direction, [m/s]

k     =  1/((1-Vf)/km  + Vf/kf);      % thermal conductivity of composite
rho_c = (1-Vf)*rhom*cm + Vf*rhof*cf;  % (density * specific heat) of composite

disp('Consider varying channel thickness:')
H1    = (Ho*Hl)  /(0.5*(Ho+Hl));
H     = H1;
H3    = (Ho*Hl)^2/(0.5*(Ho+Hl));
A     = 0.5*pi*D*W*H1*cos(theta);
C     = (W*H3)/(12*l);
E     = (pi^3 * D^3 * sin(theta))/H * (1+(sin(theta))^2); 
N1    = C/KG
N2    = (A*rho_c)/(N*l*Vo*b*E)

X     = linspace(0.01,40,4000);
F     = X.*log(X) .* (N2 - 1./ (X-1)) - N1;
figure(1)
plot(X,F), grid on, hold on
xlabel('x = e^{b \Delta T}')
ylabel('f(x) = x ln(x) (N_2 - 1/(x-1)) - N_1') 

for i = 2:length(F)
if sign(F(i))*sign(F(1)) < 0, break, end
end
x     = X(i)
str   = ['The root is at x = ',num2str(x),' when N = ',num2str(N*60),' rpm'];
plot(x,F(i),'ro','linewidth',2), text(7,0.85*max(F),str), hold off
pause(0.2)

Q(j)  = A*N - (E*C*l*N^2*b*Vo)/(rho_c*KG*x*log(x));  % flow rate, [m^3/s]
DP(j) = (Q(j)*Vo)/(KG*x);                            % pressure build up, [Pa]
DT(j) = log(x)/b;                                    % temperature rise, [oK]
Te(j) = To + DT(j);                                  % exit temperature, [oK]
Ve(j) = a*exp(-b*Te(j));                             % exit viscosity, [Pa.s]
end

figure(2), plot(Nscrew*60,Q*1.0e6,'b-o','linewidth',2), grid on
xlabel('Screw speed, N [rpm]'), ylabel('Flow rate, Q [cc/s]')

figure(3), plot(Nscrew*60,DP*1.0e-6,'b-o','linewidth',2), grid on
xlabel('Screw speed, N [rpm]'), ylabel('Pressure build-up, \Delta P [MPa]')

figure(4), plot(Nscrew*60,DT,'b-o','linewidth',2), grid on
xlabel('Screw speed, N [rpm]'), ylabel('Temperature rise, \Delta T [^oK]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%