%%%%%%%%%% Matlab code for Example 6.5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
mu    = 120000;             % viscosity of polymer, [Pa.s]
L     = 1.920;              % screw length, [m]
D     = 0.080;              % screw diameter, [m]
N     = 1;                  % screw speed, [rev/s]
Ho    = 0.015;              % inlet channel depth, [m]
Hl    = 0.002;              % exit  channel depth, [m]
theta = (20/180)*pi;        % helix angle, [rad]
KG    = 5.0e-7;             % related to die geometry and characteristics, [m^3]
e     = 0.008;              % flight width, [m]

Ls    = tan(theta)*(pi*D);  % screw lead, [m]
W     = Ls*cos(theta) - e;  % channel width, [m]
l     = L/sin(theta);       % spiral channel length, [m]
Vz    = pi*D*N*cos(theta);  % relative velocity in z direction, [m/s]
Vx    = pi*D*N*sin(theta);  % relative velocity in z direction, [m/s]

disp('Assume constant channel thickness:')
H     = 0.5*(Ho+Hl);
Fp    = 1;
Fd    = 1;
A     = 0.5*W*H*Fd*cos(theta)*pi*D;
C     = (W*H^3)/(12*l)*Fp;
Q     = (A*KG)/(C+KG)*N
DP    = (mu*A)/(C+KG)*N

disp('Consider varying channel thickness:')
slope = atan((Ho-Hl)/l) * (180/pi)     % slope of the channel, [degree]
H1    = (Ho*Hl)  /(0.5*(Ho+Hl));
H3    = (Ho*Hl)^2/(0.5*(Ho+Hl));
A     = 0.5*pi*D*W*H1*cos(theta);
C     = (W*H3)/(12*l);
DP    = DP                             % Pressure built up, [Pa]
Q     = A*N-C*DP/mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


