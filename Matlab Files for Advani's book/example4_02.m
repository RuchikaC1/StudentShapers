%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Matlab code for Example 4.2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all
R          = 8.3145; 
MU0        = [ 7.93e-14  1.06e-06];
U          = [ 9.08e+04  3.76e+04];
K          = [14.1      18.8     ];
ALPHA(1,:) = linspace(0,0.2,1001);      
ALPHA(2,:) = linspace(0,0.2,1001);
T          = [350 375 400];

for i = 1:2
    mu0   = MU0(i); u = U(i); k = K(i);
    alpha = ALPHA(i,:);
for j = 1:3
    t = T(j)
    mu(j,:) = mu0*exp(u/(R*t)+k*alpha);
end
figure(i)
subplot(1,3,1), plot(alpha,mu(1,:),'b-', 'linewidth',2), grid on%, legend('T = 350^oK',0) 
subplot(1,3,2), plot(alpha,mu(2,:),'b-', 'linewidth',2), grid on%, legend('T = 375^oK',0) 
subplot(1,3,3), plot(alpha,mu(3,:),'b-', 'linewidth',2), grid on%, legend('T = 400^oK',0) 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
