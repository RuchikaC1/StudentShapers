%%%%%%%%%% Matlab code for Example 6.6 %%%%%%%%%%%%%%
clear all, close all, clc

X = [1.2 1.3 1.5 1.7 2 2.4 3 3.5 4 5 6 8 10];

for i = 1:length(X)
x  = X(i);
N1 = logspace(log10(0.01),log10(10),1001);
N2 = N1./(x*log(x)) + 1/(x-1);
loglog(N1,N2,'k-','linewidth',2), grid on, hold on
axis([0.01 20 0.09 40])
xlabel('N_1'), ylabel('N_2')
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
