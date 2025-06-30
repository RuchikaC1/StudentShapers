%%%%%%%%%% Matlab code for Example 8.4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
Vfo      = 0.30;        % initial fiber volume fraction
mu       = 0.40;        % resin viscosity, [Pa.s]
C        = 5.0e-10;     % coefficient [m^2] in Kozeny-Carman permeability: 
                        %                      Kzz = C (1-Vf)^3/Vf^2
Va       = 0.70;        % maximum attainable fiber volume fraction
As       = 280;         % coefficient [Pa] in Fiber in a Box model:
                        %           Sigma_zz = As (sqrt(Vf/Vo)-1)/(sqrt(Va/Vf)-1)^4
PTo      = 0.0e5;       % initial autoclave pressure, [Pa]
dPTdt    = 10000;       % slope of the applied pressure ramp, [Pa/s]
h        = 0.02;        % initial part thickness, [m]
Nz       = 21;                            % # of nodes in z'
z        = linspace(0,1,Nz);              % z' array, [m]
dz       = z(2)-z(1);                     % Delta z', [m]
DT       = [0.00001 0.00002 0.00005];     % Delta t,  [s]
Vf       = Vfo*ones(Nz,1);                % fiber volume fraction
PT(1)    = PTo;                           % autoclave pressure, [Pa]
Pr       = PTo*ones(Nz,1);                % resin pressure, [Pa]
e        = (1-Vf)./Vf;
e0       = e(1);
Kzz      = C*(1-Vf).^3./Vf.^2;            % permeability, [m^2]
Sigmazz  = As*(sqrt(Vf./Vfo)-1)./(sqrt(Va./Vf)-1).^4;
t(1)     = 0;                             % time, [s]
k        = 1;                             % time index
Pr_0(1)  = Pr(1);
Pr_mid(1)= Pr(round((Nz+1)/2));
Pr_h(1)  = Pr(end);
Sigmazz_0(1)  = Sigmazz(1);
Sigmazz_mid(1)= Sigmazz(round((Nz+1)/2));
Sigmazz_h(1)  = Sigmazz(end);
Vf_0(1)  = Vf(1);
Vf_mid(1)= Vf(round((Nz+1)/2));
Vf_h(1)  = Vf(end);
H(1)     = h;
e_next   = e;
Vf_next  = Vf;
VF       = Vfo:0.0001:Va-0.0001;
SIGMAZZ  = As*(sqrt(VF./Vfo)-1)./(sqrt(Va./VF)-1).^4;
VFstop   = 0.47;

while min(min(Vf)) < VFstop
    if t(end) < 0.5, dt = DT(1); end
    if t(end) > 0.5, dt = DT(2); end
    if t(end) > 1.0, dt = DT(3); end        
    if abs(mod(t(end),0.5) - 0.5) < 0.000000001 |  mod(t(end),0.5) < 0.000000001 
    figure(1)
    subplot(311), plot(z,Pr,'k-','linewidth',1), hold on
    xlabel('Dimensionless z/h(t)'), ylabel('Resin pressure, P_r [Pa]'), grid on
    subplot(312), plot(z,Sigmazz,'k-','linewidth',1), hold on
    xlabel('Dimensionless z/h(t)'), ylabel('Fiber stress, \sigma_{zz} [Pa]'), grid on
    subplot(313), plot(z,Vf,'k-','linewidth',1), hold on
    xlabel('Dimensionless z/h(t)'), ylabel('Fiber volume fraction, V_f'), grid on
    pause(0.1)
    else
    end
h             = H(end);
t(    k+1)    = t(k) + dt;
PT(k+1)       = PT(k)+ dPTdt*dt;
e_next(1)     = (1/3)*(-e(3)+4*e(2));       % B.C. at z = 0  
c             = VF(find(abs(SIGMAZZ-PT(end)) == min(abs(SIGMAZZ-PT(end)))));          
e_next( end)  = (1-c)/c;                    % B.C. at z = h
for i = 2:Nz-1
c1            = - ( (Kzz(i+1)/(1+e(i+1))) - (Kzz(i-1)/(1+e(i-1))))/(2*dz  );
c2            =   ( Sigmazz(i+1)          - Sigmazz(i-1)         )/(2*dz  );
c3            = -    Kzz(i)/(1+e(i));
c4            =   ( Sigmazz(i+1) - 2*Sigmazz(i) + Sigmazz(i-1)   )/(  dz^2);
e_next(i)     =   e(i) + (dt*(1+e0^2)/mu) * (1/h^2) * (c1*c2+c3*c4);
end
Vf_next       = 1./(1+e_next);
fprintf('Time = %10.5f s   min.Vf = %7.4f  max.Vf = %7.4f \n',t(end),min(Vf),max(Vf))
Kzz_next      = C*(1-Vf_next.^3)./Vf_next.^2;
Sigmazz_next  = As*(sqrt(Vf_next./Vfo)-1)./(sqrt(Va./Vf_next)-1).^4;
Pr_next       = PT(k+1) - Sigmazz_next;
Pr_next(end)  = 0;
Sigmazz_next(end) = PT(k+1);
k             = k+1;
Pr_0(k)       = Pr(1);
Pr_mid(k)     = Pr(round((Nz+1)/2));
Pr_h(k)       = Pr(end);
Vf_0(k)       = Vf(1);
Vf_mid(k)     = Vf(round((Nz+1)/2));
Vf_h(k)       = Vf(end);
Sigmazz_0(k)  = Sigmazz(1);
Sigmazz_mid(k)= Sigmazz(round((Nz+1)/2));
Sigmazz_h(k)  = Sigmazz(end);
e             = e_next;
Vf            = Vf_next;
Pr            = Pr_next;
Kzz           = Kzz_next;
Sigmazz       = Sigmazz_next;
H(k)          = (H(1)*Vfo) * (1/mean(Vf));
end    

figure(2)
subplot(411), plot(t,Pr_0,'k-', t,Pr_mid,'k-.', t,Pr_h,'k--','linewidth',1) 
grid on, xlabel('Time, t [s]'), ylabel('Resin pressure, P_r [Pa]')
legend('z = 0', 'z = h/2', 'z = h', 0)
subplot(412), plot(t,Sigmazz_0,'k-', t,Sigmazz_mid,'k-.', ...
                   t,Sigmazz_h,'k--','linewidth',1) 
grid on, xlabel('Time, t [s]'), ylabel('Fiber stress, \sigma_{zz} [Pa]')
subplot(413), plot(t,Vf_0,'k-', t,Vf_mid,'k-.', t,Vf_h,'k--','linewidth',1) 
grid on, xlabel('Time, t [s]'), ylabel('Fiber volume fraction, V_f')
subplot(414), plot(t,H*1000,'k-','linewidth',1) 
grid on, xlabel('Time, t [s]'), ylabel('Part thickness, h [mm]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%