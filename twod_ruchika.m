%% 2D Darcy flow around circular preform particles
% Domain: 10 mm × 4 mm rectangle  
% Particles: 6 circles of radius 0.5 mm
% Inlet (left): P = Pin; Outlet (right): P = 0; No‐flow on top/bottom and particle walls.

close all; clear;

Pin   = 1e5;      % inlet pressure [Pa]
K     = 1e-12;    % permeability [m^2]
mu    = 0.1;      % viscosity [Pa·s]

%% 1) Create PDE model
model = createpde(1);    % scalar PDE (Darcy potential P)

%% 2) Construct geometry: rectangle minus circles
W = 10e-3;   H = 4e-3;   % domain [m]
R1 = [3;4; 0;W;W;0; 0;0;H;H];

gd = R1;
ns = 'R1';
sf = 'R1';
nCirc = 6;
rng(0);    % seed

% preallocate for centers
centers = zeros(nCirc,2);
r       = 0.5e-3;

for k = 1:nCirc
    x0 = W*(k/(nCirc+1));
    y0 = H/2 + (rand-0.5)*H*0.4;
    centers(k,:) = [x0,y0];     % store for plotting
    Ci = [1; x0; y0; r; zeros(6,1)];
    gd = [gd, Ci];
    ns = char(ns, sprintf('C%d',k));
    sf = [sf, sprintf('-C%d',k)];
end
ns = ns';
[g,bt] = decsg(gd,sf,ns);
geometryFromEdges(model,g);

%% 3) Mesh once to classify edges
mesh = generateMesh(model,'Hmax',0.05e-3);

numE = model.Geometry.NumEdges;
edgesLeft  = [];
edgesRight = [];
for eid = 1:numE
    nodeIDs = findNodes(mesh,'region','Edge',eid);
    pts     = mesh.Nodes(:,nodeIDs);
    xmean   = mean(pts(1,:));
    if abs(xmean) < 1e-6
        edgesLeft(end+1) = eid;
    elseif abs(xmean - W) < 1e-6
        edgesRight(end+1) = eid;
    end
end
edgesAll       = 1:numE;
edgesTopBottom = setdiff(edgesAll, [edgesLeft, edgesRight]);


%% 4) Boundary conditions
applyBoundaryCondition(model,'dirichlet','Edge',edgesLeft,  'u',Pin);
applyBoundaryCondition(model,'dirichlet','Edge',edgesRight, 'u',0);
applyBoundaryCondition(model,'neumann',  'Edge',edgesTopBottom,'g',0,'q',0);

%% 5) PDE coefficients
specifyCoefficients(model,'m',0,'d',0,'c',K,'a',0,'f',0);

%% 6) Solve
R = solvepde(model);
P = R.NodalSolution;

%% 7) Darcy velocity
[gradPx,gradPy] = evaluateGradient(R);
vx = -K/mu * gradPx;
vy = -K/mu * gradPy;

%% 8) Prepare interpolation grid
ngrid = 50; mgrid = 20;
[xg,yg] = meshgrid(linspace(0,W,ngrid), linspace(0,H,mgrid));
Pg  = griddata(R.Mesh.Nodes(1,:),R.Mesh.Nodes(2,:),P,   xg,yg);
vgx = griddata(R.Mesh.Nodes(1,:),R.Mesh.Nodes(2,:),vx', xg,yg);
vgy = griddata(R.Mesh.Nodes(1,:),R.Mesh.Nodes(2,:),vy', xg,yg);


%% 9) 1×2 Subplots
%% --- 1×2 Subplots: Pressure & Streamlines (no background contour on the right) ---

figure('Color','w','Position',[200 200 900 350])

% Left: pressure field + circles
subplot(1,2,1)
contourf(xg,yg,Pg,20,'LineColor','none')
colormap(jet), colorbar
title('Pressure (Pa)')
xlabel('x (m)'), ylabel('y (m)')
axis equal tight, hold on
for k = 1:nCirc
    viscircles(centers(k,:),r,'EdgeColor','k','LineWidth',1);
end
hold off

% Right: streamlines + circles only
subplot(1,2,2)
hold on
% --- after you have xg, yg, vgx, vgy, centers, r, nCirc ---

% 1) compute speed and interpolate
speed = sqrt(vx.^2 + vy.^2);    % nodal speeds on the PDE mesh
vg  = griddata(R.Mesh.Nodes(1,:), R.Mesh.Nodes(2,:), speed, xg, yg);

% 2) mask out the circle interiors so no “ghost” data inside holes
mask = false(size(xg));
for k = 1:nCirc
    mask = mask | ((xg - centers(k,1)).^2 + (yg - centers(k,2)).^2 < r^2);
end
vg(mask) = NaN;


% 
% % draw streamlines
% starty = linspace(0.1*H,0.9*H,12);
% startx = zeros(size(starty));
% hsl = streamline(xg,yg,vgx,vgy,startx,starty);
% set(hsl,'Color','k','LineWidth',1.2)
% % overlay preform circles
% for k = 1:nCirc
%     viscircles(centers(k,:),r,'EdgeColor','k','LineWidth',1);
% end
% title('Streamlines & Preform')
% xlabel('x (m)'), ylabel('y (m)')
% axis equal tight
% hold off
nodes    = mesh.Nodes';     % N×2 array of (x,y)
elements = mesh.Elements';  % M×3 array of node indices

speed = sqrt( vx.^2 + vy.^2 );  % still one value per node

figure('Color','w','Position',[200 200 600 400])
h = trisurf(...
    elements, ...
    nodes(:,1), nodes(:,2), ...  % XY coordinates
    speed, ...                    % color = speed
    'EdgeColor','none' );           % no mesh line);
view(2)       % look straight down on the XY‐plane
axis tight equal
colormap(jet), colorbar
title('Darcy Velocity Magnitude (m/s)')
xlabel('x (m)'), ylabel('y (m)')
hold on
for k = 1:nCirc
    viscircles( centers(k,:), r, 'EdgeColor','k','LineWidth',1 );
end
hold off
