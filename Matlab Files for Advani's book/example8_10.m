%%%%%%%%%% Matlab code for Example 8.10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

k        =  1/((1-Vf)/km + Vf/kf);       % thermal conductivity of composite
rho_c    = rhom*cm;                      % (density * specific heat) of resin
alpha    = k/rho_c;                      % thermal diffusivity, m^2/s
Gz       = (u*h^2)/(alpha*L);            % Graetz number
Nz       =  81;                          % # of nodes in z'
z        = linspace(-1,1,Nz);            % z' array, [m]
dz       = z(2)-z(1);                    % Delta z', [m]
r        = 0.2;
dx       = r*(Gz*dz^2);
Nx       = round(1/dx) + 1;              % # of nodes in x'
x        = linspace( 0,1,Nx);            % x' array, [m]

Tasy        = zeros(Nx,Nz);              % asymptotic temp. array
Tasy(1,: )  = Tin;                       % inlet  B.C.
Tasy(:,1 )  = To;                        % bottom B.C.
Tasy(:,Nz)  = To;                        % top    B.C.

for i = 2:Nx
for j = 2:Nz-1
    if z(j) < 0, zhat = (1+z(j))*h;
    else
                 zhat = (1-z(j))*h;
    end
    ksi    = zhat/(2*sqrt(alpha*(x(i)*L)/u));
 Tasy(i,j) = (1-erf(ksi)) * (To-Tin) + Tin;
end
end

TTC         = Tasy - 273.2;
TT_physical = TTC';
TT_physical = TT_physical(end:-1:1,:);
[X,Z]       = meshgrid(linspace(0,1,Nx)*L,linspace(1,-1,Nz)*h);
Tiso        = Tin + 0.001* (To  - Tin) - 273.2;
figure(1)
[pp,qq] = contour(X*1000,Z*1000,TT_physical,[Tiso,21,40:20:120]); clabel(pp,qq)
xlabel('x [mm]'), ylabel('z [mm]'), grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%