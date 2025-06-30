%%%%%%%%%% Matlab code for Example 8.11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
Vf       = 0.40;        % average fiber volume fraction
mu       = 0.25;        % resin viscosity, [Pa.s]
kf       = 1.20;        % thermal conductivity of fiber,  [W/(m K)]
km       = 0.20;        % thermal conductivity of matrix, [W/(m K)]
rhof     = 2500;        % density of fiber,  [kg/m3]
rhom     = 1200;        % density of matrix, [kg/m3]
cf       =  800;        % specific heat capacity of fiber,  [J/(kg K)]
cm       =  700;        % specific heat capacity of matrix, [J/(kg K)]
Tin      =  20 + 273.2; % inlet temperature, [K]
To       = 120 + 273.2; % mold wall temperature,  [K]
u        = 0.0714;      % average 1D speed, [m/s]
L        = 0.500;       % length of die, [m]
h        = 0.005;       % half-thickness of the mold cavity, [m]

U        = [u 0.10*u 0.01*u];
for icase=1:1
u        = U(icase);
k        = 1/((1-Vf)/km + Vf/kf);        % thermal conductivity of composite
rho_c    = rhom*cm;                      % (density * specific heat) of resin
alpha    = k/rho_c;                      % thermal diffusivity, m^2/s
Gz       = (u*h^2)/(alpha*L);            % Graetz number
Nz       =  81;                          % # of nodes in z'
z        = linspace(-1,1,Nz);            % z' array, [m]
dz       = z(2)-z(1);                    % Delta z', [m]
r        = 0.20;                         % the ratio, dx/(Gz dz^2); use < 0.5 
dx       = r*(Gz*dz^2);
Nx       = round(1/dx) + 1;              % # of nodes in x'
x        = linspace( 0,1,Nx);            % x' array, [m]

TH       = zeros(Nx,Nz);                 % theta array
TH(1,: ) = (Tin - Tin)/(To - Tin);       % inlet  B.C.
TH(:,1 ) = (To  - Tin)/(To - Tin);       % bottom B.C.
TH(:,Nz) = (To  - Tin)/(To - Tin);       % top    B.C.

for i = 1:Nx-1
for j = 2:Nz-1    
 TH(i+1,j)=  (1-2*r)*TH(i,j) +  r * (TH(i,j+1) + TH(i,j-1));
end
end

T           = TH   * (To - Tin) + Tin;
TC          = T - 273.2;
T_physical  = TC';
T_physical  = T_physical(end:-1:1,:);
[X,Z]       = meshgrid(linspace(0,1,Nx)*L,linspace(1,-1,Nz)*h);
Tiso        = Tin + 0.001* (To  - Tin) - 273.2;
figure(icase)
[pp,qq] = contour(X*1000,Z*1000,T_physical,[Tiso,21,40:20:120]); clabel(pp,qq)
xlabel('x [mm]'), ylabel('z [mm]'), grid
pause(0.1)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%