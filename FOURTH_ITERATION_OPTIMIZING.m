%% Opens and reads the .stl file
brain_mesh = stlread("utah_teapot.stl");

%% Function to pick a random node to become infected
%activatedNode = randi([1 size(brain_mesh.Points,1)]); % Picks a random node to set as active
patient_zero = 5;
%% Define out variables

iteration_number = 40; % Choose how many iterations/timesteps/intervals you would like to calculate the model for
node_quantity = size(brain_mesh.Points,1); % Determines how many nodes are in the model
node_values = zeros(node_quantity,iteration_number); % Matrix to be filled with values representing how active the node is (from 0 - 1)
activated_node_log = zeros(node_quantity,iteration_number); % Matrix to be filled with the row numbers each node represents
adjacent_node_list = zeros(0, 3);
%activatedNodesExpDec = activatedNodes; % Matrix to store the node values are exponential decay function applied to them

%% Sets the randomly activated node into the matrix
node_values(patient_zero,1) = 1;
activated_node_log(1,1) = patient_zero;

%% NODE PROPAGATION: Simulates the activation of neurons throughout a CSD episode. Based on the idea that from one group of neurons in the visual cortex, neuron activity propagates, spreading across the surface of the brain.
% Calculates the adjacent nodes to the activated nodes
tic
for i1 = 2:iteration_number
    progress = i1
    for j = 1:node_quantity
        if (activated_node_log(j,i1-1) ~= 0) %Checks if the node is activated before checking for adjacent nodes
            [row,col] = find(brain_mesh.ConnectivityList == activated_node_log(j,i1-1)); % The find function determines the row and column of the infected nodes in the Connectivity List
            unique_rows = unique(row);
            all_adjacent_nodes = brain_mesh.ConnectivityList(unique_rows(:,1),:); % Determines what triangulation any activated node is part of
            unique_adjacent_nodes = unique(all_adjacent_nodes);
            adjacent_node_list = [adjacent_node_list; unique_adjacent_nodes];
        end
    end

    %Puts the newly determined activated nodes into the matrix log
    %and the node value matrix.
    unique_adjacent_nodes = unique(adjacent_node_list);
    for k = 1:size(unique_adjacent_nodes)
        activated_node_log(k,i1) = unique_adjacent_nodes(k,1);
        node_values(unique_adjacent_nodes(:,1), i1) = 1;
    end
end
toc
%% Matrix to store when a node becomes activated and how many iterations have passed
%In a CSD episode after a neuron activates, its BOLD signal will drop of
% exponentially. Hence this position tracker keeps a log of what iteration
% a node becomes activated.
position_tracker = node_values;
for i2 = 1:node_quantity
    for j = 1:iteration_number
        if (position_tracker(i2,j) ~= 0)
            position_tracker(i2,j+1) = position_tracker(i2,j) + 1;
        end
    end
end

%% EXPONENTIAL DECAY: Calculates the exponential decay of neuron activity after the wave has passed over.
% Combines the position tracker matrix with the activated nodes matrix to determine the value of a node over time
% Exponential decay function
limit = 1; % The lower the limit, the more gradual the decay
tau = linspace(0, limit, iteration_number); % Creates equally spaced constants for our exponential function

for i2 = 1:node_quantity
    for j = 1:iteration_number
        if (node_values(i2,j) ~= 0)
            node_values(i2,j) = exp(-tau(1,position_tracker(i2,j)));
        end
    end
end

%% BACKGROUND NOISE: Generates artificial noise to simulate what an MRI scanner would detect in a hospital.
%Generates a matrix of "random static" to represent background noise.
%Random values picked based on normal distribution and reduced by a damping
%factor
damping = 0.2;
background_noise = damping * randn(node_quantity, iteration_number);

%Adds the background noise to the node value matrix
%If the node has a value greater than one, it is set to one. If the node has a value less
%than zero, it is set to zero
node_values = background_noise + node_values;
node_values(node_values > 1) = 1;
node_values(node_values < 0) = 0;

%% Creates and changes the settings of a figure
f = figure;
axis vis3d
axis equal

%% Plots and colours the stl file
colourMap = patch('faces', brain_mesh.ConnectivityList, 'Vertices', brain_mesh.Points);
colourMap.FaceColor = [0 0.40 0.13];
colourMap.EdgeColor = [0 0 0];
hold on

%% Animates the propagation of activated nodes by placing red markers on the infected nodes each iteration

for j = 1:iteration_number
    for i2 = 1:node_quantity
        plot3(brain_mesh.Points(i2,1), brain_mesh.Points(i2,2), brain_mesh.Points(i2,3), '.', 'markersize', 20, 'color', ([node_values(i2,j) 0 0]));
        hold on
    end
    pause(0.05); %Quickest it can go
end
