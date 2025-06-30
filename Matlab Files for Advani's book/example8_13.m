%%%%%%%%%%%%%% Matlab code for Example 8.13 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Solution to resin pressure within a L-shaped domain
clear all, close all, clc
Kxx         =  2.0e-9;     % permeability in x direction
Kyy         =  2.0e-9;     % permeability in y direction
Lx          =  1.2;        % length of mold
Ly          =  1.0;        % width  of mold
Lz          =  0.01;       % depth  of mold = inlet depth
Nx          =  25;         % # of nodes in x
Ny          =  21;         % # of nodes in y
hx          =  Lx /(Nx-1); % dx
hy          =  Ly /(Ny-1); % dy
mu          =  0.2;        % viscosity of resin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abs(hx-hy) > eps, stop, else h = hx; end
p           = zeros(Nx,Ny);
p_prev      = p;
c1          = Kxx/(2*(Kxx+Kyy));
c2          = Kyy/(2*(Kxx+Kyy));
Area        = Lz*h;
err         = inf;

while(err > 1e-5)

%%% Left half of inner points:
for i=2:12
for j=2:Ny-1
p(i,j) = c1*(p(i-1,j)+p(i+1,j)) + c2*(p(i,j-1)+p(i,j+1));
end
end

%%% Right side of inner points: 
for i=13:Nx-1
for j=12:Ny-1
p(i,j) = c1*(p(i-1,j)+p(i+1,j)) + c2*(p(i,j-1)+p(i,j+1));
end
end

p(1    ,:    ) = (4*p(2    ,:    )-p(3   ,:    ))/3;    % B.C. at left         edge
p(13   , 1:11) = (4*p(12   , 1:11)-p(11  , 1:11))/3;    % B.C. at right-bottom edge
p(Nx   ,11:21) = (4*p(Nx-1 ,11:21)-p(Nx-2,11:21))/3;    % B.C. at right-upper  edge
p(:    ,Ny   ) = (4*p(:    ,Ny-1 )-p(:   ,Ny-2 ))/3;    % B.C. at top          edge
p(1:13 ,1    ) = (4*p(1:13 ,2   )-p(1:13 ,3    ))/3;    % B.C. at bottom       edge
p(13:Nx,11   ) = (4*p(13:Nx,12  )-p(13:Nx,13   ))/3;    % B.C. at mid-horiz.   edge
p(7    ,1    ) = 500000;                                % B.C. at injection    node
p(Nx   ,17   ) = 100000;                                % B.C. at ventilation  node
p(13   ,11   ) = (4*p(12   ,12  )-p(11   ,13   ))/3;    % B.C. at kink         node
err            = norm(p-p_prev);
p_prev         = p;
end

p(14:Nx,1:10 ) = -inf; 
p_physical     = p';
p_physical     = p_physical(end:-1:1,:);

figure(1)
[X,Y]          = meshgrid(0:h:Lx,Ly:-h:0);
[pp,q]         = contour(X,Y,p_physical*0.001);clabel(pp,q)
title('Pressure Distribution [kPa]'), xlabel('x(m)'), ylabel('y(m)'), grid on
figure(2)
[X,Y]          = meshgrid(0:h:Lx,Ly:-h:0);
[pp,q]         = contourf(X,Y,p_physical*0.001);clabel(pp,q),colorbar
title('Pressure Distribution [kPa]'), xlabel('x(m)'), ylabel('y(m)'), grid on

delPdelx = zeros(Nx,Ny);
delPdely = zeros(Nx,Ny);
u        = zeros(Nx,Ny);
v        = zeros(Nx,Ny);

    for i=2:12
    for j=2:Ny-1
          delPdelx(i,j)  = (p(i+1,j  )-p(i-1,j  ))/(2*h);
          delPdely(i,j)  = (p(i  ,j+1)-p(i  ,j-1))/(2*h);    
    end
    end
 
    for i=13:Nx-1
    for j=12:Ny-1
          delPdelx(i,j)  = (p(i+1,j  )-p(i-1,j  ))/(2*h);
          delPdely(i,j)  = (p(i  ,j+1)-p(i  ,j-1))/(2*h);  
    end
    end

    for j = 2:20
        delPdelx(1,j) = (-p(3,j)+4*p(2,j)-3*p(1,j))/(2*h);    
        delPdely(1,j) = ( p(1,j+1)-p(1,j-1)       )/(2*h);   
    end

    for i = 2:24
        delPdelx(i,21) = ( p(i+1,21)-p(i-1,21)        )/(2*h);    
        delPdely(i,21) =-(-p(i,19)+4*p(i,20)-3*p(i,21))/(2*h);   
    end
   
    for j = 12:20
        delPdelx(25,j) =-(-p(23,j)+4*p(24,j)-3*p(25,j))/(2*h);    
        delPdely(25,j) = ( p(25,j+1)-p(25,j-1)        )/(2*h);   
    end
   
    for i = 14:24
        delPdelx(i,11) = ( p(i+1,11)-p(i-1,11)        )/(2*h);    
        delPdely(i,11) = (-p(i,13)+4*p(i,12)-3*p(i,11))/(2*h);   
    end      
 
    for j = 2:10
        delPdelx(13,j) =-(-p(11,j)+4*p(12,j)-3*p(13,j))/(2*h);    
        delPdely(13,j) = ( p(13,j+1)-p(13,j-1)        )/(2*h);   
    end
   
    for i = 2:12
        delPdelx(i,1) = ( p(i+1,1)-p(i-1,1)        )/(2*h);    
        delPdely(i,1) = (-p(i ,3)+4*p(i,2)-3*p(i,1))/(2*h);   
    end
   
u          = -(Kxx/mu)*delPdelx;   
v          = -(Kyy/mu)*delPdely;    
u_physical = u';  u_physical = u_physical(end:-1:1,:);
v_physical = v';  v_physical = v_physical(end:-1:1,:);

figure(3)
quiver(X,Y,u_physical,v_physical,3), hold off, axis image
title('Velocity Distribution'), xlabel('x [m]'), ylabel('y [m]')

Qinlet  = v( 7, 1)*Area    % inlet flow rate, [m^3/s]
Qexit   = u(25,17)*Area    % exit  flow rate, [m^3/s]
P_03_05 = p(7,11)          % pressure at (x,y) = (0.3,0.5), [Pa]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%