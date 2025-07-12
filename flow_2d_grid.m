% MATLAB Script for 2D Fluid Flow Around a Grid of Cylinders using PDE Toolbox
%
% This script solves the 2D incompressible Navier-Stokes equations using the
% stream-function/vorticity formulation with NO-SLIP top and bottom walls.
% It simulates flow around a user-defined m x n grid of cylinders.
% It includes a progress indicator for the solver.

% Clear workspace and close figures
clear all;
close all;
clc;

disp('Setting up the fluid dynamics model...');

%% 1. --- Model Parameters ---
Re = 1;                 % Reynolds number (e.g., 100 for vortex shedding)
avgInletVelocity = 0.2;   % Average inlet velocity for the parabolic profile
simulationTime = 50;      % Total simulation time
timeStep = 1;           % Time step size for saving results

% --- Cylinder Grid Parameters ---
cylinderRadius = 0.5;     % Radius of each cylinder in the grid
gridRows = 10;             % Number of rows in the cylinder grid (m)
gridCols = 10;             % Number of columns in the cylinder grid (n)
xSpacing = 4;            % Center-to-center spacing in the x-direction
ySpacing = 4;            % Center-to-center spacing in the y-direction

channelHeight = (gridRows+1) * ySpacing;
channelLength = (gridCols+1) * xSpacing;
porosity = pi * cylinderRadius^2 / (xSpacing * ySpacing);

% Calculate the starting position to center the grid
gridWidth = (gridCols - 1) * xSpacing;
gridHeight = (gridRows - 1) * ySpacing;
gridStartX = (channelLength - gridWidth) / 2;
gridStartY = (channelHeight - gridHeight) / 2;


%% 2. --- PDE Model Setup ---
numberOfPDE = 2;
model = createpde(numberOfPDE);

%% 3. --- Geometry Definition ---
% Define the outer rectangular channel
R1 = [3, 4, 0, channelLength, channelLength, 0, 0, 0, channelHeight, channelHeight]';

% Initialize geometry description matrices
gdm = R1;
ns = 'R001';
sf = 'R001'; % Start the set formula with the main rectangle

% Loop to create the grid of cylinders
cylinderCounter = 1;
for i = 1:gridRows
    for j = 1:gridCols
        % Calculate center coordinates for the current cylinder
        cylCenterX = gridStartX + (j - 1) * xSpacing;
        cylCenterY = gridStartY + (i - 1) * ySpacing;

        % Create the cylinder geometry definition
        C_current = [1, cylCenterX, cylCenterY, cylinderRadius]';
        C_current = [C_current; zeros(length(R1) - length(C_current), 1)]; % Pad with zeros

        % Add the cylinder to the geometry matrix
        gdm = [gdm, C_current];

        % Create a unique name for the cylinder (e.g., C1, C2, ...)
        cylName = ['C', sprintf('%03d', cylinderCounter)];
        ns = [ns; cylName];
        
        % Subtract the cylinder from the set formula
        sf = [sf, '-', cylName];
        
        cylinderCounter = cylinderCounter + 1;
    end
end
ns = ns'; % Transpose ns to the correct format

% Create the geometry from the defined shapes
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

% Inlet (Edge 4) - Dirichlet for both psi and w
applyBoundaryCondition(model, 'dirichlet', 'Edge', 4, 'u', @(region,state) inlet_bc_handle(region, avgInletVelocity, channelHeight), 'Vectorized', 'on');

% Outlet (Edge 2) - Neumann for both psi and w (allows flow to exit freely)
applyBoundaryCondition(model, 'neumann', 'Edge', 2, 'g', [0; 0]);

% --- NO-SLIP Walls ---
h_mat = [1 0; 0 0]; q_mat = [0 0; 0 0]; g_vec = [0; 0];

% Bottom Wall (Edge 1): psi = 0
applyBoundaryCondition(model, 'mixed', 'Edge', 1, 'h', h_mat, 'r', [0;0], 'g', g_vec, 'q', q_mat);

% Top Wall (Edge 3): psi = constant, value corresponds to total flow rate
psi_top = (2/3) * avgInletVelocity * channelHeight;
r_top = [psi_top; 0]; 
applyBoundaryCondition(model, 'mixed', 'Edge', 3, 'h', h_mat, 'r', r_top, 'g', g_vec, 'q', q_mat);

% --- NO-SLIP Cylinder Walls ---
disp('Applying boundary conditions to cylinder grid...');
numCylinders = gridRows * gridCols;
% The first 4 edges belong to the outer rectangle.
% Each subsequent cylinder adds 4 edges to the model.
firstCylinderEdge = 5;
edgesPerCylinder = 4;

cylinderCounter = 1;
for i = 1:gridRows
    for j = 1:gridCols
        % Calculate y-center for the current cylinder
        cylCenterY = gridStartY + (i - 1) * ySpacing;
        
        % Calculate the stream function value at the cylinder's center
        psi_cylinder_value = parabolicStreamFunction(cylCenterY, avgInletVelocity, channelHeight);
        r_cyl = [psi_cylinder_value; 0];

        % Identify the edges for the current cylinder
        startEdge = firstCylinderEdge + (cylinderCounter - 1) * edgesPerCylinder;
        endEdge = startEdge + edgesPerCylinder - 1;
        cylinderEdges = startEdge:endEdge;
        
        % Apply the boundary condition
        applyBoundaryCondition(model, 'mixed', 'Edge', cylinderEdges, 'h', h_mat, 'r', r_cyl, 'g', g_vec, 'q', q_mat);

        cylinderCounter = cylinderCounter + 1;
    end
end


%% 6. --- Initial Conditions ---
initial_cond_handle = @(location) inlet_bc_handle(location, avgInletVelocity, channelHeight);
setInitialConditions(model, initial_cond_handle);

%% 7. --- Mesh Generation ---
disp('Generating mesh...');
generateMesh(model, 'Hmax', 0.5, 'GeometricOrder', 'linear'); % Note: Hmax may need adjustment for complex grids

%% 8. --- Solve the PDE System ---
disp('Solving the time-dependent PDE...');

tlist = 0:timeStep:simulationTime;
numSteps = length(tlist) - 1;
numNodes = size(model.Mesh.Nodes, 2);
fullSolution = zeros(numNodes, numberOfPDE, length(tlist));

location.x = model.Mesh.Nodes(1,:);
location.y = model.Mesh.Nodes(2,:);
fullSolution(:,:,1) = initial_cond_handle(location)';
model.SolverOptions.RelativeTolerance = 1e-3;
model.SolverOptions.AbsoluteTolerance = 1e-5;

tic; 
for i = 1:numSteps
    tspan = [tlist(i), tlist(i+1)];
    percentComplete = (i / numSteps) * 100;
    elapsedTime = toc;
    estimatedTotalTime = (elapsedTime / i) * numSteps;
    remainingTime = estimatedTotalTime - elapsedTime;
    
    fprintf('Solver Progress: %6.2f%% | Elapsed: %s | Remaining (est.): %s\n', ...
            percentComplete, ...
            datestr(seconds(elapsedTime), 'HH:MM:SS'), ...
            datestr(seconds(remainingTime), 'HH:MM:SS'));
    
    results_step = solvepde(model, tspan);
    fullSolution(:, :, i+1) = results_step.NodalSolution(:, :, end);
    setInitialConditions(model, results_step);
end
disp('Solution complete.');

results.NodalSolution = fullSolution;
results.Tlist = tlist;


%% 9. --- Post-processing and Visualization ---

% Get the solution at the final time step
finalSolution = results.NodalSolution(:, :, end)';
psi = finalSolution(1, :); % This is a 1xNp row vector
w = finalSolution(2, :);   % This is a 1xNp row vector

% Create a grid for interpolation and plotting
x_grid = linspace(0, channelLength, 400);
y_grid = linspace(0, channelHeight, 200);
[X, Y] = meshgrid(x_grid, y_grid);

% Interpolate the stream function (psi) from the unstructured mesh nodes
x_nodes = model.Mesh.Nodes(1,:)';
y_nodes = model.Mesh.Nodes(2,:)';
psi_grid = griddata(x_nodes, y_nodes, psi', X, Y, 'natural');

% --- Simplified and Corrected Velocity Calculation ---
hx = x_grid(2) - x_grid(1);
hy = y_grid(2) - y_grid(1);

% Calculate the gradient of the gridded stream function.
% [d(psi)/dx, d(psi)/dy]
[psi_gradx_grid, psi_grady_grid] = gradient(psi_grid, hx, hy);

% Define velocity components from the stream function definition
% u = d(psi)/dy
% v = -d(psi)/dx
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

% Draw all cylinders on the plot
for i = 1:gridRows
    for j = 1:gridCols
        cylCenterX = gridStartX + (j - 1) * xSpacing;
        cylCenterY = gridStartY + (i - 1) * ySpacing;
        rectangle('Position', [cylCenterX-cylinderRadius, cylCenterY-cylinderRadius, 2*cylinderRadius, 2*cylinderRadius], 'Curvature', [1 1], 'FaceColor', 'w');
    end
end

colorbar;
title(sprintf('Velocity Magnitude at t = %.1f s (Re = %d)', simulationTime, Re));
xlabel('x'); ylabel('y');
axis equal; axis([0 channelLength 0 channelHeight]);
hold off;

%% --- ADDED: Flow Front Animation ---
disp('Starting flow front animation...');
figure('Units','normalized','OuterPosition',[0 0 1 1]);
hold on;
title('Flow Front Animation');
xlabel('x'); ylabel('y');

% Draw all cylinders for the animation background
for i = 1:gridRows
    for j = 1:gridCols
        cylCenterX = gridStartX + (j - 1) * xSpacing;
        cylCenterY = gridStartY + (i - 1) * ySpacing;
        rectangle('Position', [cylCenterX-cylinderRadius, cylCenterY-cylinderRadius, 2*cylinderRadius, 2*cylinderRadius], 'Curvature', [1 1], 'FaceColor', 'w');
    end
end
axis equal; axis([0 channelLength 0 channelHeight]);

% --- Video Setup ---
% Create a dynamic filename, e.g., "flow_animation_3x5_grid.mp4"
videoFilename = sprintf('flow_animation_%dx%d_grid.mp4', gridRows, gridCols);
disp(['Preparing video file: ' videoFilename]);
video = VideoWriter(videoFilename, 'MPEG-4');
video.FrameRate = 60;  % You can adjust the frame rate (frames per second)
open(video);

% --- Animation Parameters ---
dt = 1;
max_time = 1e5; 
num_points = 1e4; 
points = [zeros(1, num_points); linspace(0, channelHeight, num_points)];
initial_spacing = channelHeight / (num_points - 1);
max_line_distance = 10 * initial_spacing;
h_front = plot(points(1,:), points(2,:), 'k.-', 'MarkerSize', 2);
maxX = max(points(1,:));
t = 0;

% --- Animation Loop ---
while maxX < channelLength && t < max_time
    % Interpolate velocity and update particle positions
    velX = interp2(X, Y, Ux_grid, points(1,:), points(2,:), 'linear', 0);
    velY = interp2(X, Y, Uy_grid, points(1,:), points(2,:), 'linear', 0);
    points = points + [velX; velY] * dt;
    
    % Logic to break lines that stretch too far
    point_diffs = diff(points, 1, 2);
    distances = sqrt(sum(point_diffs.^2, 1));
    break_indices = find(distances > max_line_distance);
    
    plotX = points(1,:);
    plotY = points(2,:);
    for idx = fliplr(break_indices)
        plotX = [plotX(1:idx), NaN, plotX(idx+1:end)];
        plotY = [plotY(1:idx), NaN, plotY(idx+1:end)];
    end
    
    % Update plot and draw
    set(h_front, 'XData', plotX, 'YData', plotY);
    drawnow;
    
    % Capture Frame for Video
    frame = getframe(gcf);
    writeVideo(video, frame);
    
    % Update loop conditions
    maxX = max(points(1,:));
    t = t + dt;
end

hold off;

% --- Finalize Video ---
close(video);
disp(['Flow front animation complete and saved to ' videoFilename]);

%% --- End of Script ---
disp('End of script.');

%% --- Local Functions ---

function bcMatrix = parabolicInlet(region, U_avg, H)
  y = region.y;
  psi = U_avg * y.^2 .* (2*H - (4/3)*y) / H^2;
  w = -4 * U_avg * (H - 2*y) / H^2;
  bcMatrix = [psi; w];
end

function psi = parabolicStreamFunction(y, U_avg, H)
    psi = U_avg * y.^2 .* (2*H - (4/3)*y) / H^2;
end

function fMatrix = fCoeffFunc(region, state)
  f_row1 = -state.u(2,:);
  u_vel = state.uy(1, :);
  v_vel = -state.ux(1, :);
  dw_dx = state.ux(2, :);
  dw_dy = state.uy(2, :);
  f_row2 = -(u_vel .* dw_dx + v_vel .* dw_dy);
  fMatrix = [f_row1; f_row2];
end