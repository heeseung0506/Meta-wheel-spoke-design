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
ny = 76; % The number of nodes in the y-direction
nz = 25; % Number of nodes in the z-direction
num_nodes = nx * ny * nz; % Total number of nodes in 3D
num_elems = (nx-1) * (ny-1) * (nz-1); % Total number of hexahedral elements in 3D

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

%% Initialize element connectivity array
elem_connect = zeros(num_elems, 8);
elem_count = 0;
for k = 1:nz-1
    for j = 1:ny-1
        for i = 1:nx-1
            n1 = (k-1)*nx*ny + (j-1)*nx + i + 1;
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

%% Define parameters for circular layout
num_spokes = 72; % Number of spokes
circle_radius = 362; % Radius of the circle on which spokes will be placed, mm
angle_increment = 360 / num_spokes; % Angle increment per spoke in degrees

%% Initialize storage for all spokes nodes and element connectivity
all_spokes_nodes = [];
all_spokes_elements = [];
base_node_index = 0;

%% Create and transform spokes
for i = 1:num_spokes
    angle_degrees = (i-1) * angle_increment; % Current angle for the spoke
    rotated_nodes = rotateSpoke(nodes, angle_degrees); % Rotate nodes
    % Translate nodes
    translated_nodes = rotated_nodes;
    translated_nodes(:,1) = translated_nodes(:,1) + circle_radius * cosd(angle_degrees);
    translated_nodes(:,3) = translated_nodes(:,3) + circle_radius * sind(angle_degrees);
    % Append to all spokes nodes
    all_spokes_nodes = [all_spokes_nodes; translated_nodes];
    % Adjust and append element connectivity
    num_nodes_per_spoke = size(nodes, 1);
    elem_connect_spoke = elem_connect + base_node_index;
    all_spokes_elements = [all_spokes_elements; elem_connect_spoke];
    base_node_index = base_node_index + num_nodes_per_spoke;
end

%% Write the combined INP file
inp_filename = '72spokes_circular_array.inp';
% Specify an absolute path if necessary, e.g., 'C:/Users/YourUsername/Documents/MATLAB/72spokes_circular_array.inp'
fileID = fopen(inp_filename, 'w');

% Check if the file was opened successfully
if fileID == -1
    error('Failed to open file for writing. Check permissions or path validity.');
end

fprintf(fileID, '*Heading\n');
fprintf(fileID, '** 3D Model of 72 Spokes Arranged in a Circular Pattern\n');
fprintf(fileID, '*Node\n');
for i = 1:size(all_spokes_nodes, 1)
    fprintf(fileID, '%d, %f, %f, %f\n', i, all_spokes_nodes(i, 1), all_spokes_nodes(i, 2), all_spokes_nodes(i, 3));
end
fprintf(fileID, '*Element, type=C3D8R\n');
for i = 1:size(all_spokes_elements, 1)
    fprintf(fileID, '%d, %d, %d, %d, %d, %d, %d, %d, %d\n', i, all_spokes_elements(i, :));
end
fclose(fileID);

disp('3D model of 72 spokes arranged in a circular pattern generated and written to INP file.');

%% Function definitions
function rotated_nodes = rotateSpoke(nodes, angle_degrees)
    % Function to rotate nodes of a spoke around the Y-axis
    angle_radians = deg2rad(angle_degrees); % Convert degrees to radians
    % Define the rotation matrix around the Y-axis
    Ry = [cos(angle_radians) 0 sin(angle_radians);
          0 1 0;
         -sin(angle_radians) 0 cos(angle_radians)];
    % Apply the rotation to all nodes
    rotated_nodes = nodes;
    rotated_nodes(:,1) = nodes(:,1) * cos(angle_radians) - nodes(:,3) * sin(angle_radians);
    rotated_nodes(:,3) = nodes(:,1) * sin(angle_radians) + nodes(:,3) * cos(angle_radians);
end
