% Start fresh
close all; clear all; clc;

%% Define the geometry of the spoke
r = 2.0; % The initial apex height of the spoke, mm
a = 0.0; % The thickness distribution coefficient
L = 75; % The vertical height of the spoke, mm
t0 = 2.0; % The initial horizontal distance between the left and right boundary, mm
D = 60; % Depth of the spoke, mm

%% Define the mesh space
nx = 3; % The number of nodes in the x-direction
ny = 39; % The number of nodes in the y-direction
nz = 16; % Number of nodes in the z-direction
num_nodes = nx * ny * nz; % Total number of nodes in 3D
num_elems = (nx-1) * (ny-1) * (nz-1) * 1; % Total number of hexahedral elements in 3D

%% Define the coordinate of the nodes
x0 = zeros(1,ny);
t = zeros(1,ny);
x_l = zeros(1,ny);
x_r = zeros(1,ny);
X = zeros(ny,nx);
y0 = linspace(0, L, ny); % The partition of nodes in y-direction
for i = 1:ny
    x0(i) = -r*(cos(2*pi*y0(i)/L)-1); % The backbone shape function
    t(i) = 2*a*t0*(cos(2*pi*y0(i)/L))^2 + t0*(1-a); % The thickness distribution function
    x_l(i) = x0(i) - t(i)/2; % The left boundary of the spoke
    x_r(i) = x0(i) + t(i)/2; % The right boundary of the spoke
    X(i,:) = linspace(x_l(i), x_r(i), nx);
end

%% Generate the node coordinates for 3D
z0 = linspace(0, D, nz); % Z-coordinates for the layers
nodes = zeros(num_nodes, 3);
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            index = ((k-1)*nx*ny) + (j-1)*nx + i;
            nodes(index, :) = [X(j, i), y0(j), z0(k)];
        end
    end
end

%% Rotate first spoke by v degrees
nodes_spoke1 = rotateSpoke(nodes, 10);

%% Translate the second spoke in Z-direction with 10mm interval
interval = 10; % Interval between spokes in mm
max_z = max(nodes(:,3)); % Get the maximum Z value from the first spoke
translation_z = max_z + interval; % Calculate the translation distance in Z

% Duplicate the nodes for the second spoke and translate in Z-direction
nodes_spoke2 = nodes;
nodes_spoke2(:,3) = nodes_spoke2(:,3) + translation_z;

% Rotate second spoke by 20 degrees
nodes_spoke2 = rotateSpoke(nodes_spoke2, -10);

% Combine node sets
nodes_combined = [nodes_spoke1; nodes_spoke2];

%% Adjust element connectivity for the second spoke
elem_connect = zeros(num_elems, 8); % Preallocate for efficiency
elem_count = 0;
for k = 1:nz-1
    for j = 1:ny-1
        for i = 1:nx-1
            n1 = (k-1)*nx*ny + (j-1)*nx + i;
            n2 = n1 + 1;
            n3 = n1 + nx + 1;
            n4 = n1 + nx;
            n5 = n1 + nx*ny;
            n6 = n2 + nx*ny;
            n7 = n3 + nx*ny;
            n8 = n4 + nx*ny;
            elem_count = elem_count + 1;
            elem_connect(elem_count, :) = [n1, n2, n3, n4, n5, n6, n7, n8];
        end
    end
end

% Adjust element connectivity for the second spoke
elem_connect_spoke2 = elem_connect + num_nodes;

% Combine element connectivity
elem_connect_combined = [elem_connect; elem_connect_spoke2];

%% Rotate the entire model around the X-axis by 90 degrees
nodes_combined = rotateModelX(nodes_combined, 90);

%% Translate the entire model in Z and Y direction
nodes_combined(:, 3) = nodes_combined(:, 3) - 362; % Translate in Z-direction
nodes_combined(:, 2) = nodes_combined(:, 2) + 65; % Translate in Y-direction

%% Generate and rotate the spokes
total_spokes = 1;
angle_between_spokes = 5;
all_nodes = [];
all_elem_connect = [];
node_offset = 0;

for i = 0:total_spokes-1
    current_angle = i * angle_between_spokes;
    rotated_nodes = rotateNodesY(nodes_combined, current_angle);
    
    % Append the new nodes
    all_nodes = [all_nodes; rotated_nodes];
    
    % Update the element connectivity with the current node offset
    current_elem_connect = elem_connect_combined + node_offset;
    all_elem_connect = [all_elem_connect; current_elem_connect];
    
    % Update the node offset for the next iteration
    node_offset = node_offset + size(nodes_combined, 1);
end

%% Write the combined INP file
inp_filename = 'one_pair_r=2,b=5.inp';
fileID = fopen(inp_filename, 'w');
fprintf(fileID, '*Heading\n');
fprintf(fileID, '** 3D Spokes Model with Translation\n');
fprintf(fileID, '*Node\n');
for i = 1:size(all_nodes, 1)
    fprintf(fileID, '%d, %f, %f, %f\n', i, all_nodes(i, 1), all_nodes(i, 2), all_nodes(i, 3));
end
fprintf(fileID, '*Element, type=C3D8R\n');
for i = 1:size(all_elem_connect, 1)
    fprintf(fileID, '%d, %d, %d, %d, %d, %d, %d, %d, %d\n', i, all_elem_connect(i, :));
end
fclose(fileID);

disp('3D 72 spokes model with translation generated and written to INP file.');

%% Function definitions
function rotated_nodes = rotateSpoke(nodes, angle_degrees)
    % Function to rotate nodes of a spoke around the Y-axis
    angle_radians = deg2rad(angle_degrees); % Convert degrees to radians
    % Define the rotation matrix around the Y-axis
    Ry = [cos(angle_radians), 0, sin(angle_radians);
          0, 1, 0;
         -sin(angle_radians), 0, cos(angle_radians)];
    % Translate nodes to the origin
    nodes_translated = nodes - mean(nodes);
    % Apply the rotation to all nodes
    rotated_nodes_translated = (Ry * nodes_translated')'; % Transpose back to original dimensions
    % Translate nodes back to their original location
    rotated_nodes = rotated_nodes_translated + mean(nodes);
end

function rotated_nodes = rotateModelX(nodes, angle_degrees)
    % Function to rotate the entire model around the X-axis
    angle_radians = deg2rad(angle_degrees); % Convert degrees to radians
    % Define the rotation matrix around the X-axis
    Rx = [1, 0, 0;
          0, cos(angle_radians), -sin(angle_radians);
          0, sin(angle_radians), cos(angle_radians)];
    % Apply the rotation to all nodes
    rotated_nodes = (Rx * nodes')'; % Transpose back to original dimensions
end

function rotated_nodes = rotateNodesY(nodes, angle_degrees)
    % Function to rotate nodes around the Y-axis by given angle
    angle_radians = deg2rad(angle_degrees);
    Ry = [cos(angle_radians), 0, sin(angle_radians);
          0, 1, 0;
         -sin(angle_radians), 0, cos(angle_radians)];
    rotated_nodes = (Ry * nodes')';
end
