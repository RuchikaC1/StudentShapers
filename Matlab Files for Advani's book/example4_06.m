%%%%%%%%%%%%%%% Matlab code for Example 4.6 %%%%%%%%%%%%%%%%
clear all, close all, clc
R = 8.3145;
CA0 = 0.5;
k0 = 2;
E = 1000 * R;
LAMBDA = 1;
T = 400;
C(1) = 0;
t(1) = 0;
dt = 0.01;
i =1;
while C(end) < 0.8
dCdt = k0 * exp(-E/(R*T)) * CA0 * (1-C(i)) * (LAMBDA-C(i));
C(i+1) = C(i) + dCdt*dt;
t(i+1) = t(i) + dt;
i = i + 1;
end
figure(1), plot(t,C,'b-','linewidth',2)
xlabel('Time, t [seconds]'), ylabel('Extent of reaction, C^*')
grid, legend('T = 400^oK',0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%