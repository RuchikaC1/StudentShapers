
%% 2D Darcy flow around circular preform particles
% Domain: 10 mm × 4 mm rectangle  
% Particles: 6 circles of radius 0.5 mm
% Inlet (left): P = Pin; Outlet (right): P = 0; No‐flow on top/bottom and particle walls.

close all; clear;

Pin   = 1e5;      % inlet pressure [Pa]
k    = 1e-12;    % permeability [m^2]
mu    = 0.1;      % viscosity [Pa·s]
pf = 0; %pre-set pressure at flow-front

%% 1) Create PDE model
model = createpde(1);    % scalar PDE (Darcy potential P)
%% 2) Construct geometry: rectangle minus uniformly spaced circles
W = 10e-3;   H = 5e-3;   % domain [m]
R1 = [3;4; 0;W;W;0; 0;0;H;H];

gd = R1;
ns = 'R1';
sf = 'R1';
nCirc = 2;          % total number of circles
r     = 0.5e-3;     % circle radius

% specify grid dimensions: 2 rows × 3 cols = 6 circles
nRows = 2;
nCols = 1;

% preallocate
centers = zeros(nCirc,2);

count = 0;
for i = 1:nRows
    for j = 1:nCols
        count = count + 1;
        % evenly spaced in x: avoid edges by offsetting half‐cell
        x0 = (j/(nCols+1)) * W;
        % evenly spaced in y: two rows at H/3 and 2H/3
        y0 = (i/(nRows+1)) * H;
        centers(count,:) = [x0, y0];

        % build the circle description Ci = [1; x0; y0; r; zeros(6,1)];
        Ci = [1; x0; y0; r; zeros(6,1)];
        gd = [gd, Ci];

        % append name and set‐formula
        ns = char(ns, sprintf('C%d',count));      
        sf = [sf, sprintf('-C%d',count)];       
    end
end

% turn ns into char array, decompose and apply
ns = ns';
[g, bt] = decsg(gd, sf, ns);
geometryFromEdges(model, g);

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
specifyCoefficients(model,'m',0,'d',0,'c',-k/mu,'a',0,'f',0);

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
first_nodes = find(nodes(:,1)<= 0.00026);
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

while any(~is_node_full) && iteration < max_iterations % Loop until all nodes are full
    
    iteration = iteration + 1;
    %set the boundary conditions 

    % % 1. Get the x-coordinates for ONLY the nodes in your initial set.
    % x_coords = model.Mesh.Nodes(1, first_nodes);
    % 
    % % 2. Find the maximum x-coordinate among these nodes.
    % max_x = max(x_coords);
    % 
    % % 3. Find the indices of all nodes in the 'first_nodes' list
    % %    that have this maximum x-coordinate. A small tolerance is used
    % %    for a safe comparison with floating-point numbers.
    % front_thickness = 1e-7;
    % is_rightmost_node = (x_coords >= (max_x +- front_thickness));
    % 
    % % 4. The flow front is the list of node IDs from 'first_nodes'
    % %    that are on the rightmost edge.
    % flow_front_nodes_idx = first_nodes(is_rightmost_node);

    % PLACE THIS LOGIC INSIDE YOUR WHILE LOOP

    % 1. Find all nodes that are currently full.
    full_nodes = find(is_node_full);
    
    % 2. Find all elements connected to these full nodes.
    candidate_elements = unique(cat(1, elemList{full_nodes}));
    
    % 3. Get all nodes belonging to these candidate elements.
    candidate_nodes = unique(elements(candidate_elements, :));
    
    % 4. The flow front is the set of candidate nodes that are NOT yet full.
    flow_front_nodes_idx = setdiff(candidate_nodes, full_nodes);
    
    % Handle the case where the front is empty but mold isn't full (should not happen in a connected mesh)
    if isempty(flow_front_nodes_idx) && any(~is_node_full)
        disp('Warning: No flow front found, but mold is not full. Exiting.');
        break;
    end

    % figure
    % scatter(nodes(flow_front_nodes_idx,1), nodes(flow_front_nodes_idx,2))

    % Apply the Dirichlet boundary condition to the identified flow front nodes
    % --- Step 1: Assemble the raw system matrices ---
    % This creates the matrices for the PDE without any boundary conditions applied.
    fem_matrices = assembleFEMatrices(model);
    K = fem_matrices.K;
    F = fem_matrices.F;
    
    % Step 2: Manually enforce all pressure BC for the flow front
    % Enforce INLET condition
    for i = 1:length(inlet_nodes)
        node_id = inlet_nodes(i);
        K(node_id, :) = 0;
        K(node_id, node_id) = 1;
        F(node_id) = Pin;
    end
    % For each node on the front, we will modify the K and F matrices.
    for i = 1:length(flow_front_nodes_idx)
        node_id = flow_front_nodes_idx(i);
        
        % Get rid of the old equation for this node.
        K(node_id, :) = 0;
        
        % Set the new equation to be: 1 * P(node_id) = pf
        K(node_id, node_id) = 1;
        F(node_id) = pf;
    end
    
    % --- Step 3: Solve the modified linear system ---
    P = K \ F;
    R = createPDEResults(model, P);
    [gradPx,gradPy] = evaluateGradient(R);
    [gradPx_cen, gradPy_cen] = evaluateGradient(R, cen(:,1),cen(:,2) );
    vx = -k/mu * gradPx;
    vy = -k/mu * gradPy;
    vx_cen = -k/mu * gradPx_cen;
    vy_cen= -k/mu * gradPy_cen;

    ve= [vx_cen,vy_cen];
   
    for i= 1:Ne
        vid = elements(i,:); %three node IDS 
        v = nodes(vid, :); %coordinates

        %local volumetric flow rate
        I1(i) = dot(ve(i,:), normal_oc(i,:)) *len_oc(i)*t; % ve and normal are for each element 
        I2(i) = dot(ve(i,:), normal_oa(i,:))*len_oa(i)*t;
        I3(i) = dot(ve(i,:), normal_ob(i,:))*len_ob(i)*t;
    
        Q_1(i)= f(vid(1))*I1(1,i) - f(vid(2))*I2(1,i);
        Q_2(i) = f(vid(2))*I2(1,i) - f(vid(3))*I3(1,i);
        Q_3(i) = f(vid(3))*I3(1,i) - f(vid(1))*I1(1,i);

        %global flux
        Gf(vid(1)) = Gf(vid(1)) + Q_1(i);
        Gf(vid(2)) = Gf(vid(2)) + Q_2(i);
        Gf(vid(3)) = Gf(vid(3)) + Q_3(i);
    end

    delta_t_array = globalVol .* (1- f) ./ Gf;
    delta_t_array(delta_t_array <=0) = inf; % Ignore non-positive time steps

    delta_t = min(delta_t_array);
    indices = find(delta_t_array ~= inf);

    % Check for invalid time step (can happen if flow stagnates)
    if isinf(delta_t) || isnan(delta_t)
        disp('Simulation stalled. No positive flow into any cell.');
        break;
    end

    f = f + (Gf *delta_t)./globalVol; 
    f = min(f, 1.0); % Cap fill fraction at 1

    % 4. Update the status of nodes that just became full.
    %    Find nodes that are now full but weren't before.
    newly_filled_nodes = find(f >= 1 & ~is_node_full);

    if ~isempty(newly_filled_nodes)
    is_node_full(newly_filled_nodes) = true;
    % Store the fill time only once
    node_fill_time(newly_filled_nodes) = current_time; 
    end

    % Store the state for this time step
    history_entry = struct();
    history_entry.time = current_time;
    history_entry.fill_fraction = f;
    history_entry.pressure = P; % store pressure field
    flow_history{end+1} = history_entry;

end
