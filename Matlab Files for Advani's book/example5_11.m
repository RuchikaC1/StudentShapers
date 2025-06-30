%%%%%%%%%% Matlab code for Example 5.11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all, close all
mu    = 0.1;
Ri    = 0.006;
L     = 2.0;
h     = 0.004;
rmold = 0.50;
Vf    = 0.35;
k     = [6.0e-10 3.0e-10];
Patm  = 101000;
str   = 'bk';

for icase = 1:2
clear R
clear t
K    = k(icase);
R(1) = Ri;
t(1) = 0;
dt   = 0.001;
i=1;

while  R(end) < rmold
    if R(end) > 0.01, dt = 0.01; end 
    if R(end) > 0.10, dt = 0.10; end 
    a1   = (Ri^4*Patm)/(16*mu*L*(1-Vf)*h)*(1/R(i));
    a2   = 1 - 1/(1+(16*K*L*h)/(Ri^4*log(R(i)/Ri)) );
    dRdt = a1*a2;
    R(i+1) = R(i) + dRdt*dt;
    t(i+1) = t(i) + dt;
    i      = i+1;
end

t_fill(icase) = t(end)
figure(1)
plot(t,R,str(icase)), grid on, hold on
xlabel('Time, t [seconds]'), ylabel('Radius of flow front, R(t) [m]')
end
legend('K = K_1 = 6.0e-10','K = K_2 = 3.0e-10',0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%