%%%%%%%%%% Matlab code for Example 7.5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
v0      = 0.40;        % inlet fiber volume fraction
v1      = 0.50;        % exit  fiber volume fraction
Rf      = 0.008;       % radius of the fiber bundle, [m]
Vf      = (v0+v1)/2;   % average fiber volume fraction
mu      = 400000;      % resin viscosity, [Pa.s]
kf      = 1.20;        % thermal conductivity of fiber,  [W/(m K)]
km      = 0.20;        % thermal conductivity of matrix, [W/(m K)]
rhof    = 2500;        % density of fiber,  [kg/m3]
rhom    = 1200;        % density of matrix, [kg/m3]
cf      =  800;        % specific heat capacity of fiber,  [J/(kg K)]
cm      =  700;        % specific heat capacity of matrix, [J/(kg K)]
L       = 0.500;       % length of die, [m]
Tin     =  20 + 273.2; % preheater temperature, [K]
Tw      = 300 + 273.2; % die wall temperature,  [K]
U       = 0.005;       % pulling speed, [m/s]

Ri      = sqrt(Rf^2/v0);            % inlet radius of the die, [m]
Re      = sqrt(Rf^2/v1);            % exit  radius of the die, [m]
z       = linspace(0,L,1001);       % z array, [m]
dz      = z(2)-z(1);                % Delta z, [m]
R       = (Re-Ri)*(z/L)+Ri;         % radius of the die, [m]
Rd      = mean(R);                  % average radius of the die, [m]

s       = 3 - 4*(Rf/Re)^2 + (Rf/Re)^4 + 4*log(Rf/Re);
dPdz    = (8*mu*U*pi*(Re^2-R.^2))./(pi*R.^4*s);
lamd    = Rf ./ R;
QDOT    = 1/(2*mu) * dPdz.^2 .* R.^2 .* (-0.75 - log(lamd) - 0.25*lamd.^4 + lamd.^2);

P(1) = 0;
for i = 2:length(z)
P(i) = P(i-1) + dPdz(i-1)*dz;
end

figure(11), plot(z*1000,dPdz*1.0e-6,'b-','linewidth',2), grid on
xlabel('z [mm]'), ylabel('Pressure gradient, dP/dz [MPa/m]')
legend('when U = 5 mm/s and L = 0.5 m')

figure(12), plot(z*1000,P*1.0e-6,'b-','linewidth',2), grid on
xlabel('z [mm]'), ylabel('Pressure, P [MPa]')
legend('when U = 5 mm/s and L = 0.5 m')

figure(13), plot(z*1000,QDOT,'b-','linewidth',2), grid on
xlabel('z [mm]'), ylabel('Energy generation term, dq/dt [W/m^3]')
legend('when U = 5 mm/s and L = 0.5 m')

krr     =  1/((1-Vf)/km + Vf/kf);       % thermal conductivity of composite
rho_c   = (1-Vf)*rhom*cm + Vf*rhof*cf;  % (density * specific heat) of composite

Nr      =  41;                % number of nodes in r 
Nz      = 1001;               % number of nodes in z       
hr      = Rd/(Nr-1);          % Delta_r 
hz      = L /(Nz-1);          % Delta_z 
r       = linspace(0,Rd,Nr);  % r array
z       = linspace(0,L ,Nz);  % z array
T       = Tin * ones(Nz,Nr);  % initial assignment for T matrix
T(:,Nr) = Tw;                 % boundary condition along r = Rd

c1      = krr/(rho_c*U);
c2      = hz/(2*hr);
c3      = hz/(hr^2);

for i = 2:Nz
qdot  = QDOT(i);
c4    = (qdot*hz)/(rho_c * U);
   for j = 2:Nr-1    
   T(i,j)=  T(i-1,j)                                                 ...
          + c1* (   (1/r(j))*c2*(T(i-1,j+1)           -T(i-1,j-1)  ) ...
                  +          c3*(T(i-1,j+1)-2*T(i-1,j)+T(i-1,j-1)) ) ...
          + c4;
   end
   T(i,1) = (1/3)*(4*T(i,2)-T(i,3));
end

TC          = T - 273.2;
T_corner    = TC(end,1)

T_physical = TC';
T_physical = T_physical(end:-1:1,:);
[Z,R]      = meshgrid(0:hz:L, Rd:-hr:0);

figure(1)
[pp,qq] = contour(Z*1000,R*1000,T_physical,[20:20:300]); clabel(pp,qq)
title('Temperature distribution [^oC]'), xlabel('z [mm]'), ylabel('r [mm]'), grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%