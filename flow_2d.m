% MATLAB Script for 2D Fluid Flow Around a Cylinder using PDE Toolbox
%
% This script solves the 2D incompressible Navier-Stokes equations using the
% stream-function/vorticity formulation.
%
% The two coupled PDEs are:
% 1. Vorticity Transport Equation (parabolic):
%    dw/dt + u*(dw/dx) + v*(dw/dy) = (1/Re) * (d^2w/dx^2 + d^2w/dy^2)
% 2. Stream-Function Poisson Equation (elliptic):
%    d^2psi/dx^2 + d^2psi/dy^2 = -w
%
% where:
% w   = vorticity
% psi = stream function
% u,v = velocity components (u = d(psi)/dy, v = -d(psi)/dx)
% Re  = Reynolds number

% Clear workspace and close figures
clear all;
close all;
clc;

disp('Setting up the fluid dynamics model...');

% 1. --- Model Parameters ---
Re = 1;              % Reynolds number (moderate value for vortex shedding)
inletVelocity = 0.1;   % Inlet velocity in the x-direction
channelHeight = 20;   % Height of the computational domain
channelLength = 20;    % Length of the computational domain
cylinderRadius = 2;  % Radius of the cylinder
cylinderCenterX = 10; % x-coordinate of the cylinder center
cylinderCenterY = 10; % y-coordinate of the cylinder center
simulationTime = 20;   % Total simulation time
timeStep = 0.1;        % Time step size

% 2. --- PDE Model Setup ---
% We are solving a system of two PDEs for psi (u(1)) and w (u(2)).
numberOfPDE = 2;
model = createpde(numberOfPDE);

% 3. --- Geometry Definition ---
% The geometry is a rectangle with a circle subtracted.
R1 = [3, 4, 0, channelLength, channelLength, 0, 0, 0, channelHeight, channelHeight]';
C1 = [1, cylinderCenterX, cylinderCenterY, cylinderRadius]';
C1 = [C1; zeros(length(R1) - length(C1), 1)]; % Pad with zeros

% Define the geometry description matrix
gdm = [R1, C1];
% Name the shapes
ns = (['R1'; 'C1'])';
% Set formula: Rectangle minus Circle
sf = 'R1-C1';

% Create the geometry and add it to the model
geometryFromEdges(model, decsg(gdm, sf, ns));

% Plot the geometry to verify
figure;
pdegplot(model, 'EdgeLabels', 'on', 'FaceLabels', 'on');
title('Computational Domain Geometry');
axis equal;

% 4. --- PDE Coefficients ---
% Define coefficients
c = [1; 1/Re];
a = [0; 0];
d = [0; 1];
f_handle = @fCoeffFunc; 

specifyCoefficients(model, 'm', 0, 'd', d, 'c', c, 'a', a, 'f', f_handle);

% 5. --- Boundary Conditions ---
% Inlet (Edge 4) - Dirichlet for both psi and w
applyBoundaryCondition(model, 'dirichlet', 'Edge', 4, 'u', @(region,state) [inletVelocity*region.y; zeros(size(region.y))], 'Vectorized', 'on');

% Outlet (Edge 2) - Neumann for both psi and w
applyBoundaryCondition(model, 'neumann', 'Edge', 2, 'g', [0; 0]);

% --- Solid Walls (No-slip condition approximated) ---
% We use a 'mixed' boundary condition to set a Dirichlet condition on psi 
% (u(1)) and a Neumann condition on w (u(2)) in a single call for each boundary.
% General form: h*u = r (Dirichlet part) and n·(c∇u) + q*u = g (Neumann part)

% For Dirichlet on u(1): h(1,1)=1, r(1)=value.
% For Neumann on u(2):  g(2)=0, q(2,:)=0.
h_mat = [1 0; 0 0];
q_mat = [0 0; 0 0];
g_vec = [0; 0];

% Channel Walls (Edges 1, 3)
r_bottom = [0; 0]; % psi = 0 on bottom wall
applyBoundaryCondition(model, 'mixed', 'Edge', 1, 'h', h_mat, 'r', r_bottom, 'g', g_vec, 'q', q_mat);

r_top = [inletVelocity*channelHeight; 0]; % psi = const on top wall
applyBoundaryCondition(model, 'mixed', 'Edge', 3, 'h', h_mat, 'r', r_top, 'g', g_vec, 'q', q_mat);

% Cylinder Surface (Edges 5, 6, 7, 8)
r_cyl = [inletVelocity*cylinderCenterY; 0]; % psi = const on cylinder
applyBoundaryCondition(model, 'mixed', 'Edge', [5,6,7,8], 'h', h_mat, 'r', r_cyl, 'g', g_vec, 'q', q_mat);

% 6. --- Initial Conditions ---
% Start with the stream function corresponding to the inlet flow profile
% and zero initial vorticity. The function handle returns a 2x1 vector
% for the two PDE components [psi; w].
setInitialConditions(model, @(location) [inletVelocity*location.y; zeros(size(location.y))]);

% 7. --- Mesh Generation ---
disp('Generating mesh...');
hmax = 0.2; % Maximum element size
generateMesh(model, 'Hmax', hmax, 'GeometricOrder', 'linear');

% Plot the mesh
figure;
pdemesh(model);
title('Finite Element Mesh');
axis equal;

% 8. --- Solve the PDE System ---
disp('Solving the time-dependent PDE...');
tlist = 0:timeStep:simulationTime;
model.SolverOptions.ReportStatistics = 'on';
results = solvepde(model, tlist);
disp('Solution complete.');

% 9. --- Post-processing and Visualization ---

% Get the solution at the final time step
finalSolution = results.NodalSolution(:, :, end)';
psi = finalSolution(1, :); % This is a 1xNp row vector
w = finalSolution(2, :);   % This is a 1xNp row vector

% Create a grid for interpolation and plotting
x_grid = linspace(0, channelLength, 200);
y_grid = linspace(0, channelHeight, 100);
[X, Y] = meshgrid(x_grid, y_grid);

% --- Interpolation and Gradient Calculation Method ---
% Interpolate the stream function (psi) onto the structured grid first.
x_nodes = model.Mesh.Nodes(1,:)';
y_nodes = model.Mesh.Nodes(2,:)';

% Use griddata to interpolate psi. psi is 1xNp, needs to be Npx1.
psi_grid = griddata(x_nodes, y_nodes, psi', X, Y, 'natural');

% Calculate the velocity components by taking the gradient of the interpolated stream function
hx = x_grid(2) - x_grid(1);
hy = y_grid(2) - y_grid(1);
[psi_gradx_grid, psi_grady_grid] = gradient(psi_grid, hx, hy);

% u = d(psi)/dy, v = -d(psi)/dx
Ux_grid = psi_grady_grid;
Uy_grid = -psi_gradx_grid;

% Replace any NaNs that griddata might produce for points outside the convex hull
Ux_grid(isnan(Ux_grid)) = 0;
Uy_grid(isnan(Uy_grid)) = 0;

% % Mask out the cylinder region for plotting
% in_cylinder = (X - cylinderCenterX).^2 + (Y - cylinderCenterY).^2 < cylinderRadius^2;
% Ux_grid(in_cylinder) = NaN;
% Uy_grid(in_cylinder) = NaN;

% Plot Velocity Magnitude (Contourf)
% disp('Visualizing results...');
% figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);
% sgtitle('Fluid Flow Around a Cylinder Simulation');
% velocity_magnitude = sqrt(Ux_grid.^2 + Uy_grid.^2);
% contourf(X, Y, velocity_magnitude, 20, 'LineStyle', 'none');
% hold on;
% rectangle('Position', [cylinderCenterX-cylinderRadius, cylinderCenterY-cylinderRadius, 2*cylinderRadius, 2*cylinderRadius], 'Curvature', [1 1], 'FaceColor', 'w');
% colorbar;
% title('Velocity Magnitude');
% xlabel('x');
% ylabel('y');
% axis equal;
% axis([0 channelLength 0 channelHeight]);

% 10. --- Flow Front Tracing ---
disp('Tracing flow front...');
figure; % Create a new figure for the animation
hold on;
title('Flow Front Animation');
xlabel('x'); ylabel('y');
rectangle('Position', [cylinderCenterX-cylinderRadius, cylinderCenterY-cylinderRadius, 2*cylinderRadius, 2*cylinderRadius], 'Curvature', [1 1], 'FaceColor', 'w');
axis equal; axis([0 channelLength 0 channelHeight]);

% Time integration parameters
dt = 0.5;
max_time = 200;
% Define the initial line of particles (the "front")
t = 0;
num_points = 200;
initialX = zeros(1, num_points);
initialY = linspace(0, channelHeight, num_points);
points = [initialX; initialY]; % 2xN matrix of particle positions

% Define max distance for connecting lines
initial_spacing = channelHeight / (num_points - 1);
max_line_distance = 2.5 * initial_spacing;

% Plot the initial front with lines and markers
h_front = plot(points(1,:), points(2,:), 'b.-', 'MarkerSize', 8);
maxX = max(points(1,:));
while maxX < channelLength && t < max_time
    % Get the x and y coordinates of the current points
    queryX = points(1,:);
    queryY = points(2,:);
    
    % Interpolate the velocity field at the particle locations
    velX = interp2(X, Y, Ux_grid, queryX, queryY, 'linear', 0);
    velY = interp2(X, Y, Uy_grid, queryX, queryY, 'linear', 0);
    
    % Update particle positions using Euler integration
    points = points + [velX; velY] * dt;
    
    % --- Logic to break lines that are too long ---
    % Calculate the distance between adjacent points
    point_diffs = diff(points, 1, 2);
    distances = sqrt(sum(point_diffs.^2, 1));
    
    % Find indices where the line should break
    break_indices = find(distances > max_line_distance);
    
    % Create new data with NaNs inserted at the break points
    plotX = points(1,:);
    plotY = points(2,:);
    
    % Insert NaNs to break the line. We shift the indices by 1 because
    % diff reduces the array size by 1.
    for idx = fliplr(break_indices) % Go backwards to not mess up indices
        plotX = [plotX(1:idx), NaN, plotX(idx+1:end)];
        plotY = [plotY(1:idx), NaN, plotY(idx+1:end)];
    end
    
    % Update the plot with the new front position
    set(h_front, 'XData', plotX, 'YData', plotY);
    drawnow;
    
    % Update loop conditions
    maxX = max(points(1,:));
    t = t + dt;
end
hold off;
disp('Flow front tracing complete.');


disp('End of script.');

% --- Local Function for PDE Coefficients ---
function fMatrix = fCoeffFunc(region, state)
  % This function computes the 'f' vector for the PDE system.
  % It must return a 2-by-Nr matrix, where Nr is the number of evaluation points.

  % f1 for the first equation (stream function) -> -w
  f_row1 = -state.u(2,:);

  % f2 for the second equation (vorticity) -> convection term
  % u_vel = d(psi)/dy = state.uy(1,:)  (psi is component 1)
  % v_vel = -d(psi)/dx = -state.ux(1,:) (psi is component 1)
  % dw_dx = state.ux(2,:)              (w is component 2)
  % dw_dy = state.uy(2,:)              (w is component 2)
  % Convection term f2 = -(u_vel * dw_dx + v_vel * dw_dy)
  u_vel = state.uy(1, :);
  v_vel = -state.ux(1, :);
  dw_dx = state.ux(2, :);
  dw_dy = state.uy(2, :);
  
  f_row2 = -(u_vel .* dw_dx + v_vel .* dw_dy);

  % Combine into a 2xNr matrix
  fMatrix = [f_row1; f_row2];
end
