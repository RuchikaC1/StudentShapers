%%%%%%%%%% Matlab code for Example 7.2 %%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc
v0   = 0.40;
v1   = 0.50;
Rf   = 0.020;
LL   = linspace(0.1,0.5,41);
mu   = 100000;
P0   = 0;
Pmax = 1.0e6;
Ri   = sqrt(Rf^2/v0);
Re   = sqrt(Rf^2/v1);

for i = 1:length(LL)
L    = LL(i);
z    = linspace(0,L,1001);
dz   = z(2)-z(1);
R    = (Re-Ri)*(z/L)+Ri; 
f    = (Re^2-R.^2)./(R.^4);
I    = dz*(1/3)*(  f(1) + 4*sum(f(2:2:end-1)) ...
                        + 2*sum(f(2:2:end-1)) ...
                 + f(end)                    );
s    = 3 - 4*(Rf/Re)^2 + (Rf/Re)^4 + 4*log(Rf/Re);
U(i) = (Pmax*s)/(8*mu*I);

if i == 11
    k = 1;
    for j = 1:10:length(z)
    zf = z(j);
    zz = linspace(0,zf,1001);
    RR = (Re-Ri)*(zz/L)+Ri; 
    dz = zz(2)-zz(1);
    ff = (Re^2-RR.^2)./(RR.^4);
    II = dz*(1/3)*( ff(1) + 4*sum(ff(2:2:end-1)) ...
                          + 2*sum(ff(2:2:end-1)) ...
                  + ff(end)                     );
    PP(k) = P0 + (8*mu*U(i))/s*II;
    ZZ(k) = z(j);
    k     = k+1;
    end
    figure(1)
    plot(ZZ,PP,'b-','linewidth',2), grid on
    xlabel('z [m]')
    ylabel('Pressure, P [Pa]')
    legend('when P_{max} = 1 MPa and L = 0.2 m')
end

end
    figure(2)
    plot(LL,U*1000,'bo-',LL(11),U(11)*1000,'rs','linewidth',2)
    grid on
    xlabel('Die Length, L [m]')
    ylabel('Pulling Speed, U [mm/s]')
    legend('when P_{max} = 1 MPa')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%