% =========================================================================
% 1. Problem Setup & Initialization
% =========================================================================
clear; clc; close all;
fprintf('--- Starting Simulation Setup ---\n');
% --- Physical and Numerical Parameters for Stokes Flow---
rho = 1.0;            % Fluid density [kg/m^3]
mu = 100;             % High viscosity for Stokes Flow
P_max = 1;          % Pressure for slow, creeping flow
pressure_ramptime = 1; % Time in seconds to reach max pressure
H = 5;                % Half-height of the channel
dt = 0.1;             % Time step [s]
t_end = 10.0;         % Total simulation time [s]
alpha_p = 0.1;        % Pressure relaxation factor
fprintf('Parameters loaded.\n');

% --- Geometry and Mesh ---
model = createpde(2);
R1 = [3, 4, -10, 10, 10, -10, -H, -H, H, H]';
C1 = [1, 0, 0, 1]';
C1 = [C1; zeros(length(R1) - length(C1), 1)];
geom = decsg([R1, C1], 'R1-C1', ['R1';'C1']');
geometryFromEdges(model, geom);
fprintf('Geometry created.\n');
inletEdge = 4;
outletEdge = 2;
cylinderEdges = [5, 6, 7, 8];
wallEdges = [1, 3];
h_cylinder = 0.3;
edge_pairs = {cylinderEdges, h_cylinder};
generateMesh(model, 'Hmax', 0.6, 'Hedge', edge_pairs, 'GeometricOrder', 'linear');
fprintf('Mesh generated successfully.\n');

% --- Boundary Conditions ---
applyBoundaryCondition(model, 'dirichlet', 'Edge', [wallEdges, cylinderEdges], 'u', [0, 0]);
applyBoundaryCondition(model, 'neumann', 'Edge', [inletEdge, outletEdge], 'g', [0;0]);
fprintf('Boundary conditions applied.\n');

% --- Initialization ---
nodes = model.Mesh.Nodes;
N = size(nodes, 2);
u_prev = zeros(N, 1);
v_prev = zeros(N, 1);
p_prev = zeros(N, 1);
inlet_x_coordinate = -10;
outlet_x_coordinate = 10;
tolerance = 1e-6;
inlet_node_indices = find(abs(model.Mesh.Nodes(1,:) - inlet_x_coordinate) < tolerance);
outlet_node_indices = find(abs(model.Mesh.Nodes(1,:) - outlet_x_coordinate) < tolerance);
fprintf('Fields initialized.\n');

% --- Storage for Animation ---
anim_frames = {};
figure('Position', [100, 100, 1200, 500]);

% =========================================================================
% 2. Time-Stepping Loop (SIMPLE Method)
% =========================================================================
fprintf('\n--- Starting Transient Simulation ---\n');
time_steps = 0:dt:t_end;
for t = time_steps
    fprintf('Solving for time t = %.2f s / %.2f s\n', t, t_end);
    
    P_inlet_t = P_max * 0.5 * (1 + tanh((t - pressure_ramptime/2) / (pressure_ramptime/4)));
    if t >= pressure_ramptime
        P_inlet_t = P_max;
    end
    
    [grad_p_x, grad_p_y] = evaluateGradient(model.Mesh, p_prev);
    f_handle = @(loc,s) predictor_f_function_Stokes(loc, s, rho, dt, u_prev, v_prev, grad_p_x, grad_p_y, nodes);
    specifyCoefficients(model, 'm', 0, 'd', 0, 'c', mu, 'a', rho/dt, 'f', f_handle);
    results_predictor = solvepde(model);
    u_star = results_predictor.NodalSolution(:, 1);
    v_star = results_predictor.NodalSolution(:, 2);
    
    [du_star_dx, ~] = evaluateGradient(model.Mesh, u_star);
    [~, dv_star_dy] = evaluateGradient(model.Mesh, v_star);
    div_u_star = du_star_dx + dv_star_dy;
    
    model_p = createpde(1);
    model_p.Geometry = model.Geometry;
    model_p.Mesh = model.Mesh;
    
    % --- Pressure Correction BCs ---
    % Inlet is Dirichlet to drive the flow
    p_prime_inlet_target = (P_inlet_t - p_prev(inlet_node_indices)) / alpha_p;
    inlet_nodes_y = model.Mesh.Nodes(2, inlet_node_indices);
    [sorted_inlet_y, sort_idx_in] = sort(inlet_nodes_y);
    inlet_bc_func = @(loc,s) interp1(sorted_inlet_y, p_prime_inlet_target(sort_idx_in), loc.y, 'linear', 'extrap');
    
    % MODIFIED: Outlet is now Dirichlet (fixed pressure) again.
    % This is stable for the Stokes flow model.
    p_prime_outlet_target = (0 - p_prev(outlet_node_indices)) / alpha_p;
    outlet_nodes_y = model.Mesh.Nodes(2, outlet_node_indices);
    [sorted_outlet_y, sort_idx_out] = sort(outlet_nodes_y);
    outlet_bc_func = @(loc,s) interp1(sorted_outlet_y, p_prime_outlet_target(sort_idx_out), loc.y, 'linear', 'extrap');
    
    applyBoundaryCondition(model_p, 'dirichlet', 'Edge', inletEdge, 'u', inlet_bc_func);
    applyBoundaryCondition(model_p, 'dirichlet', 'Edge', outletEdge, 'u', outlet_bc_func);
    applyBoundaryCondition(model_p, 'neumann', 'Edge', [wallEdges, cylinderEdges], 'g', 0);
    
    specifyCoefficients(model_p, 'm', 0, 'd', 0, 'c', 1, 'a', 0, 'f', @(loc,s) pressure_f_function(loc, s, rho, dt, div_u_star, nodes));
    results_p_prime = solvepde(model_p);
    p_prime = results_p_prime.NodalSolution;
    
    % --- Step 3: Correct Velocity and Pressure ---
    [grad_p_prime_x, grad_p_prime_y] = evaluateGradient(model.Mesh, p_prime);
    u_new = u_star - (dt / rho) * grad_p_prime_x;
    v_new = v_star - (dt / rho) * grad_p_prime_y;
    p_new = p_prev + alpha_p * p_prime;
    
    % NOTE: The pressure pinning step has been removed.
    
    % --- Step 4: Update for Next Time Step & Visualize ---
    u_prev = u_new;
    v_prev = v_new;
    p_prev = p_new;
    
    subplot(1, 2, 1);
    pdeplot(model, 'XYData', sqrt(u_new.^2 + v_new.^2), 'Contour', 'on');
    title(sprintf('Velocity Magnitude (m/s) at t = %.2f s', t));
    axis equal;
    colormap(gca, 'jet');
    caxis([0, max(max(sqrt(u_new.^2+v_new.^2)), 0.1)]); % Adjusted color axis for visibility
    colorbar;
    
    subplot(1, 2, 2);
    pdeplot(model, 'XYData', p_new, 'Contour', 'on', 'ColorMap', 'jet');
    title(sprintf('Pressure (Pa) at t = %.2f s', t));
    axis equal;
    colorbar;
    
    drawnow;
    anim_frames{end+1} = getframe(gcf);
end
fprintf('\n--- Simulation Finished ---\n');

% =========================================================================
% 3. Save Animation & Helper Functions
% =========================================================================
try
    fprintf('Saving animation...\n');
    writerObj = VideoWriter('transient_flow_around_cylinder.mp4', 'MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);
    for i = 1:length(anim_frames)
        writeVideo(writerObj, anim_frames{i});
    end
    close(writerObj);
    fprintf('Animation saved to "transient_flow_around_cylinder.mp4"\n');
catch e
    fprintf('Could not save animation. Error: %s\n', e.message);
end

% --- HELPER FUNCTION FOR STOKES FLOW ---
function f = predictor_f_function_Stokes(location, state, rho, dt, u_prev, v_prev, grad_px, grad_py, nodes)
    F_u = scatteredInterpolant(nodes(1,:)', nodes(2,:)', u_prev);
    F_v = scatteredInterpolant(nodes(1,:)', nodes(2,:)', v_prev);
    F_gpx = scatteredInterpolant(nodes(1,:)', nodes(2,:)', grad_px);
    F_gpy = scatteredInterpolant(nodes(1,:)', nodes(2,:)', grad_py);

    u = F_u(location.x, location.y);
    v = F_v(location.x, location.y);
    gpx = F_gpx(location.x, location.y);
    gpy = F_gpy(location.x, location.y);
    
    nr = length(location.x);
    f = zeros(2, nr);
    
    f(1, :) = (rho/dt) * u.' - gpx.';
    f(2, :) = (rho/dt) * v.' - gpy.';
end

function f = pressure_f_function(location, state, rho, dt, div_u_star, nodes)
    F_div = scatteredInterpolant(nodes(1,:)', nodes(2,:)', div_u_star);
    div_loc = F_div(location.x, location.y);
    f = (rho/dt) * div_loc;
end

function [gradx, grady] = evaluateGradient(mesh, data)
    [gradx_c, grady_c] = pdegrad(mesh.Nodes, mesh.Elements, data);
    gradx = pdeprtni(mesh.Nodes, mesh.Elements, gradx_c);
    grady = pdeprtni(mesh.Nodes, mesh.Elements, grady_c);
end