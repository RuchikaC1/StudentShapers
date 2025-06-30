%%%%%%%%% Matlab code for Example 5.6 %%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all, close all
%%% INPUT: 
k_r     =    0.1;   % in [W/(m C)]
k_f     =   10.0;  
Rho_r   =  900;     % in [kg / m^3]
Rho_f   = 2500;
C_p_r   = 2000;     % in [J / (kg C)]
C_p_f   =  410;
Vf      =    0.50;  % fiber volume fraction
h_thick =    0.010; % half thickness in [m]
T_w     =   25;     % temperature in [C]
T_i     =  200;

RhoCp_r = Rho_r * C_p_r;
RhoCp_f = Rho_f * C_p_f;
k       = (1-Vf)*k_r   + Vf*k_f;
RhoCp   = (1-Vf)*RhoCp_r + Vf*RhoCp_f;
alpha   = k/(RhoCp);             % thermal diffusivity
t_c     = h_thick^2 / alpha;     % reference time
DeltaT  = (T_w - T_i);

r          = 0.1; % = Delta t^ / (Delta z^) ^2
N          = 101; % Number of increments in z direction
%%% I.C.:
theta      = zeros(1,N);
theta(1)   = 1;
theta(end) = 1;
theta_next = theta;
dz_non     = (1-(-1))/(N-1);
dt_non     = r * (dz_non)^2;
dt         = dt_non * t_c;
z_non      = linspace(-1,1,N);

t_non  = 0.0;
t      = 0.0;
k_time = 1;
mid    = (N+1)/2;
TIME_non(k_time)    = 0;
TIME(k_time)        = 0;
THETA(k_time)       = 0;
TEMPERATURE(k_time) = T_i;

figure(1)
while theta(mid) < (70 - T_i)/(T_w - T_i)
        
for i = 2:N-1
theta_next(i) = r*theta(i-1) + (1-2*r)*theta(i) + r*theta(i+1);
end

theta               = theta_next;
temperature         = theta * DeltaT + T_i;
t_non               = t_non + dt_non;
t                   = t     + dt;
k_time              = k_time + 1;
TIME_non(k_time)    = t_non;
TIME(k_time)        = t;
THETA(k_time)       = theta(mid);
TEMPERATURE(k_time) = temperature(mid);

    if (mod(k_time,round(1/dt)) == 1)
    text_plot = ['time [seconds] = ',num2str(t,3)];
    plot(temperature,z_non*h_thick,'linewidth',2); grid on
    axis([0 1.1*T_i -1.2*h_thick 1.2*h_thick])
    text(0.55*T_i,0.95*h_thick,text_plot);
    xlabel('Temperature, T(z)  [C]')
    ylabel('Height, h [m]')
    pause(0.10)
    else
    end
end

t_non_critical = t_non(end); % non-dime. time when T_mid = 70 C

figure(2)
plot(TIME_non,THETA,'linewidth',2)
grid
xlabel('Nondimensional Time')
ylabel('Nondimensional Temperature of Mid Point')

%%%% Investigate the effect of k_f
%%%% INPUT: 
k_r     =    0.1;   % in [W/(m C)]
k_f     =   linspace(1,20,191);  
Rho_r   =  900;     % in [kg / m^3]
Rho_f   = 2500;
C_p_r   = 2000;     % in [J / (kg C)]
C_p_f   =  410;
Vf      =    0.50;  % fiber volume fraction
h_thick =    0.010; % half thickness in [m]
T_w     =   25;     % temperature in [C]
T_i     =  200;

RhoCp_r = Rho_r * C_p_r;
RhoCp_f = Rho_f * C_p_f;
k       = (1-Vf)*k_r   + Vf*k_f;
RhoCp   = (1-Vf)*RhoCp_r + Vf*RhoCp_f;
alpha   = k./(RhoCp);            % thermal diffusivity
t_c     = h_thick^2 ./ alpha;    % reference time
DeltaT  = (T_w - T_i);

t_critical     = t_non_critical*t_c;
figure(3)
plot(k_f,t_critical,'b-','linewidth',2)
grid
xlabel('Thermal Conductivity of Fiber, k_f [W/(mC)]')
ylabel('Critical Time, t_{T70} [s] when T_{mid} = 70 degrees')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%