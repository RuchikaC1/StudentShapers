%%%%%%%%%%%%%% Matlab code for Example 8.12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Solution to resin pressure within a rectangular domain
clear all, close all, clc
Kxx         =  1.0e-9;     % permeability in x direction
Kyy         =  0.2e-9;     % permeability in y direction
Lx          =  0.8;        % length of mold
Ly          =  0.2;        % width  of mold
Lz          =  0.02;       % depth  of mold = inlet depth
Nx          =  41;         % # of nodes in x
Ny          =  11;         % # of nodes in y
hx          =  Lx /(Nx-1); % dx
hy          =  Ly /(Ny-1); % dy
mu          =  0.1;        % viscosity of resin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abs(hx-hy) > eps, stop, else h = hx; end
p           = zeros(Nx,Ny);
p_prev      = p;

c1          = Kxx/(2*(Kxx+Kyy));
c2          = Kyy/(2*(Kxx+Kyy));
err         = 1000;

while (err > 0.01)

for  i = 2:Nx-1
for  j = 2:Ny-1
p(i,j) = c1*(p(i-1,j)+p(i+1,j)) + c2*(p(i,j-1)+p(i,j+1));
end
end

p(1 ,: ) = (4*p(2   ,:)-p(3   ,:))/3;    % B.C. at left   edge
p(Nx,: ) = (4*p(Nx-1,:)-p(Nx-2,:))/3;    % B.C. at right  edge
p(: ,1 ) = (4*p(:,2   )-p(:,3   ))/3;    % B.C. at bottom edge
p(: ,Ny) = (4*p(:,Ny-1)-p(:,Ny-2))/3;    % B.C. at top    edge
p(1 , 6) = 300000;                       % B.C. at injection   node
p(41, 6) =      0;                       % B.C. at ventilation node

err      = norm(p-p_prev);
p_prev   = p; 
end

p_physical = p';
p_physical = p_physical(end:-1:1,:);
[X,Y]      = meshgrid(0:h:Lx, Ly:-h:0);

figure(1)
[pp,qq] = contour(X,Y,p_physical*0.001,[0:30:300]); clabel(pp,qq)
title('Pressure distribution [kPa]'),xlabel('x [m]'),ylabel('y [m]'),grid

delPdelx = zeros(Nx,Ny);
delPdely = zeros(Nx,Ny);
u        = zeros(Nx,Ny);
v        = zeros(Nx,Ny);
      
   for i = 2:Nx-1
   for j = 2:Ny-1                  
          delPdelx(i,j)  = (p(i+1,j  )-p(i-1,j  ))/(2*h);
          delPdely(i,j)  = (p(i  ,j+1)-p(i  ,j-1))/(2*h);    
   end     
   end
   
   for i = 2:Nx-1
        delPdelx(i,1 ) = (p(i+1,1   )-p(i-1,1   ))/(2*h);    
        delPdelx(i,Ny) = (p(i+1,Ny  )-p(i-1,Ny  ))/(2*h);    
   end

   for j = 2:Ny-1
        delPdely(1 ,j) = (p(1   ,j+1)-p(1   ,j-1))/(2*h); 
        delPdely(Nx,j) = (p(Nx  ,j+1)-p(Nx  ,j-1))/(2*h);    
   end
  
   for i = 2:Nx-1
        delPdely(i,1 ) = (-p(i, 1+2)+4*p(i, 1+1)-3*p(i,   1))/(2*h);    
        delPdely(i,Ny) =-(-p(i,Ny-2)+4*p(i,Ny-1)-3*p(i,Ny  ))/(2*h);       
   end

   for j = 2:Ny-1
        delPdelx(1 ,j) = (-p( 1+2,j)+4*p( 1+1,j)-3*p(   1,j))/(2*h);    
        delPdelx(Nx,j) =-(-p(Nx-2,j)+4*p(Nx-1,j)-3*p(Nx  ,j))/(2*h);    
   end
   
u          = -(Kxx/mu)*delPdelx;   
v          = -(Kyy/mu)*delPdely;    
u_physical = u';
u_physical = u_physical(end:-1:1,:);
v_physical = v';
v_physical = v_physical(end:-1:1,:);

figure(3)
quiver(X,Y,u_physical,v_physical,3), hold off, axis image
title('Velocity Distribution'), xlabel('x [m]'), ylabel('y [m]')

Qinlet  = u( 1,6)*(h*Lz)
Qexit   = u(41,6)*(h*Lz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%