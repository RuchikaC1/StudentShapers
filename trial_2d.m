
%% 2D Darcy flow around circular preform particles
% Domain: 10 mm × 4 mm rectangle  
% Particles: 6 circles of radius 0.5 mm
% Inlet (left): P = Pin; Outlet (right): P = 0; No‐flow on top/bottom and particle walls.

close all; clear;

Pin   = 1e5;      % inlet pressure [Pa]
mu    = 100;      % viscosity [Pa·s]
pf = 0; %pre-set pressure at flow-front

Vf_target = 0.4; % The desired fiber volume fraction (60%)
W = 4e-3;       % Total width of the preform
H = 3e-3;        % Total height of the preform
r = 0.1e-3;      % Radius of one circular fiber

kx    = ((1-Vf_target)^3) / (4.2*(6/(2*r))^2 *Vf_target^2);    % permeability [m^2]
ky = kx;
total_area = W * H;
area_one_fiber = pi * r^2;

% Calculate the theoretical number of circles
nCirc_float = (Vf_target * total_area) / area_one_fiber;

% You can't have a fraction of a circle, so round to the nearest whole number
nCirc = round(nCirc_float);

%% 1) Create PDE model
model = createpde(1);    % scalar PDE (Darcy potential P)
%% 2) Construct geometry: rectangle minus uniformly spaced circles

W = 4.5e-3;   H = 3e-3;   % domain [m]
R1 = [3;4; 0;W;W;0; 0;0;H;H];
gd = R1;
ns = 'R1';
sf = 'R1';
% 
% % --- Call the new random placement function ---
% % 'nCirc' is calculated earlier in your script based on Vf_target
% centers = placeFibersRandomly(nCirc, W, H, r);
% 
% % --- Build the geometry description from the returned centers ---
% % Update nCirc in case the function returned fewer circles than requested
% nCirc_placed = size(centers, 1); 
% 
% for i = 1:nCirc_placed
%     x0 = centers(i,1);
%     y0 = centers(i,2);
% 
%     Ci = [1; x0; y0; r; zeros(6,1)];
%     gd = [gd, Ci];
%     ns = char(ns, sprintf('C%d',i));      
%     sf = [sf, sprintf('-C%d',i)];       
% end

% %% -- USER-DEFINED PARAMETERS -- %%
r     = r;      % circle radius [m]
nRows = 10;           % number of rows
nCols = 15;           % number of columns

% Define manual spacing and offsets [m]
x_offset  = 0.5e-3;    % distance from left edge to first circle's center
y_offset  = 0.27e-3;     % distance from bottom edge to first circle's center
x_spacing = 0.25e-3;    % horizontal spacing between circle centers
y_spacing = 0.27e-3;     % vertical spacing between circle centers
% %% ---------------------------- %%

nCirc = nRows * nCols; % total number of circles
centers = zeros(nCirc,2);
count = 0;

for i = 1:nRows
    for j = 1:nCols
        count = count + 1;
        % Calculate positions based on manual offsets and spacing
        x0 = x_offset + (j-1) * x_spacing;
        y0 = y_offset + (i-1) * y_spacing;

        % --- The rest of the code is the same ---
        centers(count,:) = [x0, y0];
        Ci = [1; x0; y0; r; zeros(6,1)];
        gd = [gd, Ci];
        ns = char(ns, sprintf('C%d',count));      
        sf = [sf, sprintf('-C%d',count)];       
    end
end

ns = ns';
[g, bt] = decsg(gd, sf, ns);
geometryFromEdges(model, g); % This line would then use the generated geometry

%% 3) Mesh once to classify edges
mesh = generateMesh(model,'Hmax',0.1e-3, 'GeometricOrder','linear');

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
nodes =mesh.Nodes';
elements = mesh.Elements' ; 

% Plot the mesh
figure;
pdemesh(model);
title('Finite Element Mesh');
axis equal;

%% 5) flow-cell development

Nnodes = size(nodes,1);
Ne= size(elements,1);

elemIdx = repmat((1:Ne)',3,1);
nodeIdx = elements(:); %all node IDS, element by element

elemList = accumarray(nodeIdx, elemIdx, [Nnodes,1], @(x){x}); % this tells us for each ode, the list of elements that includes it

%defining local and global flow cells
t = 1; %unit thickness

globalVol = zeros(Nnodes,1);

for e= 1:Ne
    vid = elements(e,:); %three node IDS 
    v = nodes(vid, :); %coordinates
    cen(e,:) = mean(v, 1); %compute triangle centroid

    %compute midpoints of each edge
    m12 = (v(1,:) + v(2,:)) / 2;
    m23 = (v(2,:) + v(3,:)) / 2;
    m31 = (v(3,:) + v(1,:)) / 2;

    % local polygon around node 1: [v1 → m12 → cen → m31]
    poly1 = [ v(1,:) ; m12 ; cen(e,:) ; m31 ];
    A1     = polyarea(poly1(:,1), poly1(:,2));

    % similarly for node 2: [v2 → m23 → cen → m12]
    poly2 = [ v(2,:) ; m23 ; cen(e,:) ; m12 ];
    A2    = polyarea(poly2(:,1), poly2(:,2));

    % and node 3: [v3 → m31 → cen → m23]
    poly3 = [ v(3,:) ; m31 ; cen(e,:) ; m23 ];
    A3    = polyarea(poly3(:,1), poly3(:,2));

    % sanity check: these should sum to the triangle area
    Atri = polyarea(v(:,1), v(:,2));
    assert( abs((A1+A2+A3) - Atri) < 1e-12 );

    % accumulate into global volumes (area × thickness)
    globalVol(vid(1)) = globalVol(vid(1)) + A1*t;
    globalVol(vid(2)) = globalVol(vid(2)) + A2*t;
    globalVol(vid(3)) = globalVol(vid(3)) + A3*t;

    %get the normal vectors of each of the local cell segment
    % Define segments from centroid to midpoints
    vec_oa(e,:) = m23 - cen(e,:); % Segment 'oa' from text (centroid to midpoint opposite node 1)
    vec_ob(e,:) = m31 - cen(e,:); % Segment 'ob' from text (centroid to midpoint opposite node 2)
    vec_oc(e,:) = m12 - cen(e,:); % Segment 'oc' from text (centroid to midpoint opposite node 3)
    
    % Calculate lengths
    len_oa(e,:) = norm(vec_oa(e,:));
    len_ob(e,:) = norm(vec_ob(e,:));
    len_oc(e,:) = norm(vec_oc(e,:));

    % Define OUTWARD normal for each segment (away from centroid).
    % Consistent counter-clockwise rotation of the segment vector.
    normal_oa(e,:) = [-vec_oa(e,2), vec_oa(e,1)]; normal_oa(e,:) = normal_oa(e,:) / norm(normal_oa(e,:));
    normal_ob(e,:) = [-vec_ob(e,2), vec_ob(e,1)]; normal_ob(e,:) = normal_ob(e,:) / norm(normal_ob(e,:));
    normal_oc(e,:) = [-vec_oc(e,2), vec_oc(e,1)]; normal_oc(e,:) = normal_oc(e,:) / norm(normal_oc(e,:));

end

%% 4) standard boundary conditions
applyBoundaryCondition(model,'dirichlet','Edge',edgesLeft,  'u',Pin);
applyBoundaryCondition(model,'dirichlet','Edge',edgesRight, 'u',0);
applyBoundaryCondition(model,'neumann',  'Edge',edgesTopBottom,'g',0,'q',0);
%% Flow Analysis Network
% 5) PDE coefficients
% 'c' is now a vector for anisotropic permeability
c_coeff = [-kx/mu; -ky/mu];
specifyCoefficients(model,'m',0,'d',0,'c',c_coeff,'a',0,'f',0);
% specifyCoefficients(model,'m',0,'d',0,'c',-k/mu,'a',0,'f',0);

%6) Initialise all arrays 

%Store the final fill time for each node
%Initialize to infinty to signify "not yet filled"
node_fill_time = inf(Nnodes,1);

% A flag to track which nodes are completely full
is_node_full = false(Nnodes, 1);

%fill factor each node's global flow cell
f= zeros(Nnodes, 1);

current_time = 0; %start simulation time

%find the nodes at x=0
inlet_nodes = find(nodes(:,1) == 0);
f(inlet_nodes)=1;
is_node_full(inlet_nodes)=true;

%Initialise the first nodes 
first_nodes = find(nodes(:,1)<= 0.00046);
% Remove any nodes from the front list that are also inlet nodes.
first_nodes = setdiff(first_nodes, inlet_nodes);
f(first_nodes)= 1;
is_node_full(first_nodes)= true;   

%%initial flow pressure etc
max_iterations = 5000; % Safety break to prevent infinite loops
iteration = 0;
flow_history = {}; % Cell array to store results for animation
 %calculating the local and global flux
Gf = zeros(Nnodes,1);

%% 4) MAIN SIMULATION SUPER-LOOP (HANDLES TRAPPED ISLANDS)
fprintf('Starting main simulation loop...\n');
super_loop_count = 0;
while any(~is_node_full)
    super_loop_count = super_loop_count + 1;
    fprintf('\n--- Starting Filling Phase %d. Filled nodes: %d/%d ---\n', super_loop_count, sum(is_node_full), Nnodes);
    
    % --- Inner Simulation Loop (runs one phase until stall) ---
    while iteration < max_iterations
        iteration = iteration + 1;
        
        % DYNAMIC FRONT DEFINITION: The front is all partially filled nodes.
        % This is the most robust way to define the front.
        % 1. Find all nodes that are currently full.
        full_nodes = find(is_node_full);
    
        % 2. Find all elements connected to these full nodes.
        candidate_elements = unique(cat(1, elemList{full_nodes}));
    
        % 3. Get all nodes belonging to these candidate elements.
        candidate_nodes = unique(elements(candidate_elements, :));
    
        % 4. The flow front is the set of candidate nodes that are NOT yet full.
        flow_front_nodes_idx = setdiff(candidate_nodes, full_nodes);

         % The full active domain for the pressure solve
        active_nodes = union(full_nodes, flow_front_nodes_idx);
        
        % STALL CHECK: If no nodes are partially filled, the front has vanished.
        if isempty(flow_front_nodes_idx)
            disp('Stall detected: No partially filled nodes found. Breaking to check for islands..., no front');
            break; % Exit inner loop to let the super-loop handle it.
        end
        
        % % --- SOLVE FOR PRESSURE ---
        % fem_matrices = assembleFEMatrices(model);
        % K = fem_matrices.K;
        % F = fem_matrices.F;
        % 
        % % Step 2: Manually enforce all pressure BC for the flow front
        % % Enforce INLET condition
        % for i = 1:length(inlet_nodes)
        %     node_id = inlet_nodes(i);
        %     K(node_id, :) = 0;
        %     K(node_id, node_id) = 1;
        %     F(node_id) = Pin;
        % end
        % % For each node on the front, we will modify the K and F matrices.
        % for i = 1:length(flow_front_nodes_idx)
        %     node_id = flow_front_nodes_idx(i);
        % 
        %     % Get rid of the old equation for this node.
        %     K(node_id, :) = 0;
        % 
        %     % Set the new equation to be: 1 * P(node_id) = pf
        %     K(node_id, node_id) = 1;
        %     F(node_id) = pf;
        % end
        % P = K \ F;
        % 
        % % --- CALCULATE FLUXES ---
        % R = createPDEResults(model, P);
        % 
        % [gradPx_cen, gradPy_cen] = evaluateGradient(R, cen(:,1), cen(:,2));
        % ve = [-k/mu * gradPx_cen, -k/mu * gradPy_cen];
        % Gf = zeros(Nnodes, 1);
        % I1=zeros(Ne,1); I2=zeros(Ne,1); I3=zeros(Ne,1);
        % Q_1=zeros(Ne,1); Q_2=zeros(Ne,1); Q_3=zeros(Ne,1);
        % 
        % for i = 1:Ne
        %     vid = elements(i,:);
        %     I1(i) = dot(ve(i,:), normal_oc(i,:)) * len_oc(i) * t;
        %     I2(i) = dot(ve(i,:), normal_oa(i,:)) * len_oa(i) * t;
        %     I3(i) = dot(ve(i,:), normal_ob(i,:)) * len_ob(i) * t;
        %     Q_1(i) = f(vid(1))*I1(i) - f(vid(2))*I2(i);
        %     Q_2(i) = f(vid(2))*I2(i) - f(vid(3))*I3(i);
        %     Q_3(i) = f(vid(3))*I3(i) - f(vid(1))*I1(i);
        %     Gf(vid(1)) = Gf(vid(1)) + Q_1(i);
        %     Gf(vid(2)) = Gf(vid(2)) + Q_2(i);
        %     Gf(vid(3)) = Gf(vid(3)) + Q_3(i);
        % end

            % --- 2. SOLVE FOR PRESSURE IN THE ACTIVE DOMAIN (THE CORRECT WAY) ---
    fem_matrices = assembleFEMatrices(model);
    K_glob = fem_matrices.K;
    F_glob = fem_matrices.F;
    P = zeros(Nnodes, 1); % Initialize global pressure vector

    % Create a smaller sub-system for only the active nodes
    K_sub = K_glob(active_nodes, active_nodes);
    F_sub = F_glob(active_nodes);

    % Find the local indices within the sub-system for applying BCs
    [~, inlet_nodes_local_idx] = ismember(inlet_nodes, active_nodes);
    [~, front_nodes_local_idx] = ismember(flow_front_nodes_idx, active_nodes);
    inlet_nodes_local_idx(inlet_nodes_local_idx == 0) = []; % remove non-active inlets if any
    
    % Apply Dirichlet BCs to the SUB-SYSTEM
    for i = 1:length(inlet_nodes_local_idx)
        idx = inlet_nodes_local_idx(i);
        K_sub(idx, :) = 0;
        K_sub(idx, idx) = 1;
        F_sub(idx) = Pin;
    end
    for i = 1:length(front_nodes_local_idx)
        idx = front_nodes_local_idx(i);
        K_sub(idx, :) = 0;
        K_sub(idx, idx) = 1;
        F_sub(idx) = pf;
    end
    
    % Solve the smaller, correct system and map solution back to global P
    P_sub = K_sub \ F_sub;
    P(active_nodes) = P_sub;

    % --- 3. CALCULATE FLUXES (CORRECTED FVM LOGIC) ---
    R = createPDEResults(model, P);
    [gradPx_cen, gradPy_cen] = evaluateGradient(R, cen(:,1), cen(:,2));
    ve = [-kx/mu * gradPx_cen, -ky/mu * gradPy_cen];
    
    Gf = zeros(Nnodes, 1); % Reset global flux
    for i = 1:Ne
        vid = elements(i,:);
        P_nodes = P(vid);
        f_nodes = f(vid);
        
        % Fluxes across the three internal faces of the element's control volumes

        I_m12 = dot(ve(i,:), normal_oc(i,:)) *len_oc(i)*t; % ve and normal are for each element 
        I_m23 = dot(ve(i,:), normal_oa(i,:))*len_oa(i)*t;
        I_m31 = dot(ve(i,:), normal_ob(i,:))*len_ob(i)*t;
        
        % Upwinding: The flux carries the fill factor of the upstream node
        f12 = (P_nodes(1) > P_nodes(2)) * f_nodes(1) + (P_nodes(2) >= P_nodes(1)) * f_nodes(2);
        f23 = (P_nodes(2) > P_nodes(3)) * f_nodes(2) + (P_nodes(3) >= P_nodes(2)) * f_nodes(3);
        f31 = (P_nodes(3) > P_nodes(1)) * f_nodes(3) + (P_nodes(1) >= P_nodes(3)) * f_nodes(1);

        % Accumulate net flux for each node from this element's contribution
        Gf(vid(1)) = Gf(vid(1)) + f31*I_m31 - f12*I_m12;
        Gf(vid(2)) = Gf(vid(2)) + f12*I_m12 - f23*I_m23;
        Gf(vid(3)) = Gf(vid(3)) + f23*I_m23 - f31*I_m31;

        % Q_1(i)= f(vid(1))*I1(i) - f(vid(2))*I2(i);
        % Q_2(i) = f(vid(2))*I2(i) - f(vid(3))*I3(i);
        % Q_3(i) = f(vid(3))*I3(i) - f(vid(1))*I1(i);
        % 
        %global flux
        % Gf(vid(1)) = Gf(vid(1)) + Q_1(i);
        % Gf(vid(2)) = Gf(vid(2)) + Q_2(i);
        % Gf(vid(3)) = Gf(vid(3)) + Q_3(i);
    end
        
        % --- ADVANCE TIME STEP ---
        vol_to_fill = globalVol .* (1 - f);
        delta_t_array = vol_to_fill ./ Gf;
        delta_t_array(delta_t_array <= 1e-9) = inf;
        [delta_t, node_that_should_fill_idx] = min(delta_t_array);

        % --- FIX: Add a cap to prevent huge time jumps ---
        % max_delta_t = 0.5; % Set a reasonable maximum time step (e.g., 0.1 seconds)
        % delta_t = min(delta_t);
        
        if isinf(delta_t) || isnan(delta_t)
            disp('Stall detected: No positive flow. Breaking to check for islands...');
            break; % Exit inner loop
        end
        
        current_time = current_time + delta_t;
        f = f + (Gf * delta_t) ./ globalVol;
        
        % FORCE FILL: Guarantees progress against floating-point errors.
        if ~isempty(node_that_should_fill_idx)
            f(node_that_should_fill_idx) = 1.0;
        end
        
        % --- UPDATE NODE STATUS ---
        f(is_node_full) = 1.0; % Lock already full nodes
        f = max(f, 0.0);       % Prevent f from going negative
        newly_filled_nodes = find(f >= 0.999 & ~is_node_full);
        if ~isempty(newly_filled_nodes)
            is_node_full(newly_filled_nodes) = true;
            node_fill_time(newly_filled_nodes) = current_time;
        end
        
        % --- STORE HISTORY ---
        if mod(iteration, 1) == 0 || ~isempty(newly_filled_nodes)
            history_entry=struct(); 
            history_entry.time=current_time; 
            history_entry.fill_fraction=f; 
            history_entry.pressure = P; 
            flow_history{end+1}=history_entry;
            fprintf('Iter: %d, Time: %.4f s, Filled: %d/%d\n', iteration, current_time, sum(is_node_full), Nnodes);
        end
    end % --- End of inner simulation loop ---

    % --- ISLAND DETECTION AND HANDLING LOGIC (after a stall) ---

    if sum(is_node_full) >= Nnodes || iteration >= max_iterations
        disp('Mold is full or max iterations reached. Exiting.');
        break; % Exit the super-loop
    end
    
    fprintf('Finding largest unfilled island...\n');
    unfilled_nodes = find(~is_node_full);
    
    adj = sparse(elements(:,[1 2 2 3 3 1]), elements(:,[2 1 3 2 1 3]), 1, Nnodes, Nnodes);
    adj_unfilled = adj(unfilled_nodes, unfilled_nodes);
    G_unfilled = graph(adj_unfilled, 'omitselfloops');
    
    % Get the connected components (islands)
    [bins, ~] = conncomp(G_unfilled);

    if isempty(bins)
        fprintf('No more islands to fill. Finishing.\n');
        f(~is_node_full) = 1.0; is_node_full(~is_node_full) = true;
        continue;
    end
    
    % --- FIX: Handle the case of a single island ---
    % Find the largest island by counting occurrences in the 'bins' vector
    unique_bins = unique(bins);
    bin_counts = histc(bins, unique_bins);
    [~, max_idx] = max(bin_counts);
    largest_island_bin_id = unique_bins(max_idx);
    
    % Get the node indices for the largest island
    nodes_in_largest_island = unfilled_nodes(bins == largest_island_bin_id);
    
    boundary_of_island = [];
    for node_id = nodes_in_largest_island'
        neighbors = find(adj(node_id, :));
        if any(is_node_full(neighbors))
            boundary_of_island(end+1) = node_id;
        end
    end
    
    if isempty(boundary_of_island)
        fprintf('Could not identify a new boundary. Cannot proceed.\n');
        break;
    end
    
    fprintf('Identified new front with %d nodes. Restarting fill...\n', length(boundary_of_island));
    f(unique(boundary_of_island)) = 1e-5; % Kickstart the new front
    
end % --- End of outer "super-loop" ---
disp('Simulation Finished.');

%% 5) SIMPLER ANIMATION USING A SCATTER PLOT
% load('simulation_results_3.mat')
% disp('Creating simple scatter plot animation...');
% 
% % --- 1. Set up the VideoWriter object ---
% video_filename = 'flow_animation.mp4';
% v = VideoWriter(video_filename, 'MPEG-4');
% v.FrameRate = 10; % Adjust frame rate (frames per second) as needed
% open(v); % Open the file for writing
% 
% figure('Name', 'Simple Flow Animation');
% nodes = model.Mesh.Nodes'; % Get node coordinates (Nx2 matrix)
% Nnodes = size(nodes, 1);
% 
% for n = 1:length(flow_history)
% 
%     % --- CORRECTED CHECK ---
%     % Plot if the frame number 'n' is a multiple of 800, OR if it's the last frame.
%     if mod(n, 10) == 0 || n == length(flow_history)
% 
%         current_f = flow_history{n}.fill_fraction;
%         current_time = flow_history{n}.time;
%         current_f(~isfinite(current_f)) = 0; % Clean data as before
% 
%         % --- Define colors based on fill status ---
%         node_colors = repmat([0.8, 0.8, 0.8], Nnodes, 1); 
%         full_nodes_idx = (current_f >= 0.999);
%         node_colors(full_nodes_idx, :) = repmat([0, 0.4, 0.8], sum(full_nodes_idx), 1);
%         % front_nodes_idx = (current_f > 1e-6 & current_f < 0.999);
%         % node_colors(front_nodes_idx, :) = repmat([1, 0, 0], sum(front_nodes_idx), 1);
% 
%         % --- Create the Plot ---
%         clf; % Clear the frame
%         scatter(nodes(:,1), nodes(:,2), 10, node_colors, 'filled');
% 
%         axis equal;
%         title(sprintf('Flow Front Motion at Time = %.4f s (Frame %d)', current_time, n));
%         xlabel('Domain Length (m)');
%         ylabel('Domain Height (m)');
% 
%         % Add a legend
%         hold on;
%         p1 = plot(NaN,NaN,'o','MarkerFaceColor',[0.8, 0.8, 0.8],'MarkerEdgeColor','k');
%         p2 = plot(NaN,NaN,'o','MarkerFaceColor',[0, 0.4, 0.8],'MarkerEdgeColor','k');
%         p3 = plot(NaN,NaN,'o','MarkerFaceColor',[1, 0, 0],'MarkerEdgeColor','k');
%         legend([p1,p2], {'Empty','Filled'}, 'Location', 'southeast');
%         hold off;
% 
%         drawnow;
%          % --- 2. Capture the current figure as a video frame ---
%         frame = getframe(gcf); % gcf gets the current figure handle
%         writeVideo(v, frame); % Write the frame to the video file
% 
%     end % --- END OF THE IF STATEMENT ---
% end
% 
% close(v);
% 
% disp('Animation finished.');

disp('Creating contour plot animation...');

% --- 1. Get mesh data (do this once before the loop) ---
nodes = model.Mesh.Nodes';      % Node coordinates (Nx2 matrix)
elements = model.Mesh.Elements'; % Element connectivity (Mx3 matrix)

% --- 2. Set up the VideoWriter object ---
video_filename = 'flow_contour_animation.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 15; % Increased frame rate for smoother animation
open(v);

% --- 3. Loop through history and create frames ---
figure('Name', 'Flow Front Contour Animation');
for n = 1:length(flow_history)
    
    % Plot every 10th frame to speed up video generation
    if mod(n, 10) == 0 || n == length(flow_history)
    
        % Get the data for the current time step
        current_f = flow_history{n}.fill_fraction;
        current_time = flow_history{n}.time;
        
        % --- Create the Plot using patch() ---
        clf; % Clear the current figure
        
        patch('Faces', elements, 'Vertices', nodes, ...
              'FaceVertexCData', current_f, ... % Data at each vertex
              'FaceColor', 'interp', ...       % Interpolate color across the face
              'EdgeColor', 'none');            % Hide the mesh lines for a smooth look
        
        % --- Formatting ---
        axis equal;
        caxis([0 1]); % Set consistent color limits (0=unfilled, 1=filled)
        h_bar = colorbar; % Add a color bar
        ylabel(h_bar, 'Fill Fraction (f)'); % Add a label to the color bar
        
        title(sprintf('Flow Front Motion at Time = %.4f s', current_time));
        xlabel('Domain Length (m)');
        ylabel('Domain Height (m)');
        
        drawnow;
        
        % --- 4. Capture the frame for the video ---
        frame = getframe(gcf);
        writeVideo(v, frame);
    
    end
end

close(v); % Finalize and save the video
fprintf('Animation saved to %s\n', video_filename);
%% 
save('simulation_results.mat', 'flow_history', 'model');

function centers = placeFibersRandomly(nCirc, W, H, r)
    % Places nCirc circular fibers of radius r randomly within a WxH domain.
    % Ensures that no fibers overlap with each other or the boundaries.
    %
    % INPUTS:
    %   nCirc - Number of circles (fibers) to place.
    %   W     - Width of the rectangular domain [m].
    %   H     - Height of the rectangular domain [m].
    %   r     - Radius of each circle [m].
    %
    % OUTPUT:
    %   centers - An [nCirc x 2] matrix of the [x, y] coordinates for the centers.

    fprintf('Attempting to randomly place %d fibers...\n', nCirc);
    
    centers = zeros(nCirc, 2);
    max_attempts_per_fiber = 5000; % Safety break

    for i = 1:nCirc
        is_placed = false;
        attempts = 0;
        
        while ~is_placed && attempts < max_attempts_per_fiber
            attempts = attempts + 1;
            
            % Generate a random center coordinate.
            % The valid area for centers is inset by 'r' from the main boundary.
            x_candidate = r + (W - 2*r) * rand();
            y_candidate = r + (H - 2*r) * rand();
            
            % Assume the position is valid until we find an overlap.
            is_valid_position = true;
            
            % Check for overlap with previously placed circles.
            % We only need to check up to the circles we've already placed (i-1).
            for j = 1:(i-1)
                dist_sq = (x_candidate - centers(j,1))^2 + (y_candidate - centers(j,2))^2;
                
                % If the squared distance is less than the squared diameter, they overlap.
                if dist_sq < (2*r)^2
                    is_valid_position = false;
                    break; % No need to check other circles, move to a new attempt.
                end
            end
            
            % If the position is still valid after all checks, place the circle.
            if is_valid_position
                centers(i, :) = [x_candidate, y_candidate];
                is_placed = true;
            end
        end
        
        % If a fiber couldn't be placed, warn the user and exit.
        if ~is_placed
            warning('Could not place fiber #%d after %d attempts. The domain may be too full. Returning %d placed fibers.', i, max_attempts_per_fiber, i-1);
            centers = centers(1:i-1, :); % Return only the successfully placed fibers
            return;
        end
    end
    fprintf('Successfully placed all %d fibers.\n', nCirc);
end