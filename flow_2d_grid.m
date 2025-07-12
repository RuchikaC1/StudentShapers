% MATLAB Script for 2D Fluid Flow Around a Cylinder using PDE Toolbox
%
% This script solves the 2D incompressible Navier-Stokes equations using the
% stream-function/vorticity formulation with NO-SLIP top and bottom walls.
% It includes a progress indicator for the solver.

% Clear workspace and close figures
clear all;
close all;
clc;

disp('Setting up the fluid dynamics model...');

%% 1. --- Model Parameters ---
Re = 1;              % Reynolds number (e.g., 100 for vortex shedding)
avgInletVelocity = 0.1; % Average inlet velocity for the parabolic profile
channelHeight = 20;    % Height of the computational domain
channelLength = 20;    % Length of the computational domain
cylinderRadius = 5;    % Radius of the cylinder
cylinderCenterX = 10;  % x-coordinate of the cylinder center
cylinderCenterY = 10;  % y-coordinate of the cylinder center
simulationTime = 50;  % Total simulation time (e.g., 200 for flow development)
timeStep = 0.5;        % Time step size for saving results

%% 2. --- PDE Model Setup ---
numberOfPDE = 2;
model = createpde(numberOfPDE);

%% 3. --- Geometry Definition ---
R1 = [3, 4, 0, channelLength, channelLength, 0, 0, 0, channelHeight, channelHeight]';
C1 = [1, cylinderCenterX, cylinderCenterY, cylinderRadius]';
C1 = [C1; zeros(length(R1) - length(C1), 1)]; % Pad with zeros
gdm = [R1, C1];
ns = (['R1'; 'C1'])';
sf = 'R1-C1';
geometryFromEdges(model, decsg(gdm, sf, ns));

% Plot the geometry to verify
figure;
pdegplot(model, 'EdgeLabels', 'on', 'FaceLabels', 'on');
title('Computational Domain Geometry');
axis equal;

%% 4. --- PDE Coefficients ---
c = [1; 1/Re];
a = [0; 0];
d = [0; 1];
f_handle = @fCoeffFunc; 
specifyCoefficients(model, 'm', 0, 'd', d, 'c', c, 'a', a, 'f', f_handle);

%% 5. --- Boundary Conditions ---
% Define a function handle for the parabolic inlet profile
inlet_bc_handle = @parabolicInlet;

% Inlet (Edge 4) - Dirichlet for both psi and w using the parabolic profile
applyBoundaryCondition(model, 'dirichlet', 'Edge', 4, 'u', @(region,state) inlet_bc_handle(region, avgInletVelocity, channelHeight), 'Vectorized', 'on');

% Outlet (Edge 2) - Neumann for both psi and w (allows flow to exit freely)
applyBoundaryCondition(model, 'neumann', 'Edge', 2, 'g', [0; 0]);

% --- NO-SLIP Walls and Cylinder ---
h_mat = [1 0; 0 0]; q_mat = [0 0; 0 0]; g_vec = [0; 0];

% Bottom Wall (Edge 1): psi = 0
applyBoundaryCondition(model, 'mixed', 'Edge', 1, 'h', h_mat, 'r', [0;0], 'g', g_vec, 'q', q_mat);

% Top Wall (Edge 3): psi = constant, value corresponds to total flow rate
psi_top = (2/3) * avgInletVelocity * channelHeight;
r_top = [psi_top; 0]; 
applyBoundaryCondition(model, 'mixed', 'Edge', 3, 'h', h_mat, 'r', r_top, 'g', g_vec, 'q', q_mat);

% cylinder
psi_cylinder_value = parabolicStreamFunction(cylinderCenterY, avgInletVelocity, channelHeight);
r_cyl = [psi_cylinder_value; 0]; % Create a constant r-vector
applyBoundaryCondition(model, 'mixed', 'Edge', [5,6,7,8], 'h', h_mat, 'r', r_cyl, 'g', g_vec, 'q', q_mat);

%% 6. --- Initial Conditions ---
% Initialize the entire domain with the parabolic flow profile
initial_cond_handle = @(location) inlet_bc_handle(location, avgInletVelocity, channelHeight);
setInitialConditions(model, initial_cond_handle);

%% 7. --- Mesh Generation ---
disp('Generating mesh...');
generateMesh(model, 'Hmax', 0.5, 'GeometricOrder', 'linear');

%% 8. --- Solve the PDE System ---
disp('Solving the time-dependent PDE...');

% Set up a manual time-stepping loop to show progress
tlist = 0:timeStep:simulationTime;
numSteps = length(tlist) - 1;

% Pre-allocate a matrix to store the solution at each requested time step
numNodes = size(model.Mesh.Nodes, 2);
fullSolution = zeros(numNodes, numberOfPDE, length(tlist));

% Evaluate the initial conditions function at the mesh nodes to get the starting values.
location.x = model.Mesh.Nodes(1,:);
location.y = model.Mesh.Nodes(2,:);
fullSolution(:,:,1) = initial_cond_handle(location)';


% Start a timer to track progress
tic; 

for i = 1:numSteps
    % Define the time span for this single step
    tspan = [tlist(i), tlist(i+1)];

    % Display Progress
    percentComplete = (i / numSteps) * 100;
    elapsedTime = toc;
    estimatedTotalTime = (elapsedTime / i) * numSteps;
    remainingTime = estimatedTotalTime - elapsedTime;
    
    fprintf('Solver Progress: %6.2f%% | Elapsed: %s | Remaining (est.): %s\n', ...
            percentComplete, ...
            datestr(seconds(elapsedTime), 'HH:MM:SS'), ...
            datestr(seconds(remainingTime), 'HH:MM:SS'));
    
    % Solve the PDE for this small time span
    results_step = solvepde(model, tspan);
    
    % Store the final solution of this step
    fullSolution(:, :, i+1) = results_step.NodalSolution(:, :, end);
    
    % IMPORTANT: Update the initial conditions for the NEXT loop iteration
    setInitialConditions(model, results_step);
end

disp('Solution complete.');

% Recreate the results object in the format expected by the post-processing code
results.NodalSolution = fullSolution;
results.Tlist = tlist;


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

% Plot Velocity Magnitude Contour
disp('Visualizing velocity field...');
figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);
velocity_magnitude = sqrt(Ux_grid.^2 + Uy_grid.^2);
contourf(X, Y, velocity_magnitude, 20, 'LineStyle', 'none');
hold on;
rectangle('Position', [cylinderCenterX-cylinderRadius, cylinderCenterY-cylinderRadius, 2*cylinderRadius, 2*cylinderRadius], 'Curvature', [1 1], 'FaceColor', 'w');
colorbar;
title(sprintf('Velocity Magnitude at t = %.1f s (Re = %d)', simulationTime, Re));
xlabel('x'); ylabel('y');
axis equal; axis([0 channelLength 0 channelHeight]);
hold off;

%% --- ADDED: Flow Front Animation ---
disp('Starting flow front animation...');
figure('Units','normalized','OuterPosition',[0 0 1 1]); % Maximize figure window
hold on;
title('Flow Front Animation');
xlabel('x'); ylabel('y');
rectangle('Position', [cylinderCenterX-cylinderRadius, cylinderCenterY-cylinderRadius, 2*cylinderRadius, 2*cylinderRadius], ...
          'Curvature', [1 1], 'FaceColor', 'w');
axis equal; axis([0 channelLength 0 channelHeight]);

% Time integration parameters for the animation
dt = 0.25;
max_time = 300; % Set animation time
% Define the initial line of particles (the "front")
num_points = 200;
points = [zeros(1, num_points); linspace(0, channelHeight, num_points)];
% Define max distance for connecting lines
initial_spacing = channelHeight / (num_points - 1);
max_line_distance = 10 * initial_spacing;

% Plot the initial front
h_front = plot(points(1,:), points(2,:), 'k.-', 'MarkerSize', 2);
maxX = max(points(1,:));
t = 0;

% Animation loop
while maxX < 2*channelLength && t < max_time
    % Interpolate the velocity field at the particle locations
    velX = interp2(X, Y, Ux_grid, points(1,:), points(2,:), 'linear', 0);
    velY = interp2(X, Y, Uy_grid, points(1,:), points(2,:), 'linear', 0);
    
    % Update particle positions using simple Euler integration
    points = points + [velX; velY] * dt;
    
    % --- Logic to break lines that are too long (e.g., when splitting around cylinder) ---
    point_diffs = diff(points, 1, 2);
    distances = sqrt(sum(point_diffs.^2, 1));
    break_indices = find(distances > max_line_distance);
    
    plotX = points(1,:);
    plotY = points(2,:);
    for idx = fliplr(break_indices)
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
disp('Flow front animation complete.');
%% --- End of Script ---
disp('End of script.');

%% --- Local Functions ---

function bcMatrix = parabolicInlet(region, U_avg, H)
  % Calculates psi and w for a fully developed parabolic flow.
  y = region.y;
  
  % Stream function: psi(y) = integral of u(y) dy
  % where u(y) = 4 * U_avg * y .* (H - y) / H^2
  psi = U_avg * y.^2 .* (2*H - (4/3)*y) / H^2;
  
  % Vorticity: w(y) = -du/dy
  w = -4 * U_avg * (H - 2*y) / H^2;
  
  bcMatrix = [psi; w];
end

function psi = parabolicStreamFunction(y, U_avg, H)
    % Helper function to calculate just the stream function value for a
    % given y-position in the parabolic flow profile.
    psi = U_avg * y.^2 .* (2*H - (4/3)*y) / H^2;
end

function fMatrix = fCoeffFunc(region, state)
  % This function computes the 'f' vector for the PDE system.
  % f = [ -w; -(u*dw/dx + v*dw/dy) ]
  f_row1 = -state.u(2,:); % f1 = -w
  
  % Convection term for vorticity transport equation
  u_vel = state.uy(1, :);  % u = d(psi)/dy
  v_vel = -state.ux(1, :); % v = -d(psi)/dx
  dw_dx = state.ux(2, :);
  dw_dy = state.uy(2, :);
  f_row2 = -(u_vel .* dw_dx + v_vel .* dw_dy);

  fMatrix = [f_row1; f_row2];
end