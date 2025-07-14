%% Fresh start
clear
close all
%% Opens and reads the .stl file
brain_mesh = stlread("stl/dans_brain.stl");
[success_chime, Fs] = audioread("sound/UI_CHIME_SOUND.mp3");

%% Define our variables

iteration_number = 60; % Choose how many iterations/timesteps/intervals you would like to calculate the model for
node_quantity = size(brain_mesh.Points,1); % Determines how many nodes are in the model
node_values = zeros(node_quantity,iteration_number); % Matrix to be filled with values representing how active the node is (from 0 - 1)
activated_node_log = zeros(node_quantity,iteration_number); % Matrix to be filled with the row numbers each node represents
adjacent_node_list = zeros(0, 3);
start_node = 1531; %V3A
end_node = 5715; %V5
node_values(start_node,1) = 1;
activated_node_log(1,1) = start_node;

%% NODE PROPAGATION: Simulates the activation of neurons throughout a CSD episode. Based on the idea that from one group of neurons in the visual cortex, neuron activity propagates, spreading across the surface of the brain.
% Calculates the adjacent nodes to the activated nodes
tic
for i = 2:iteration_number
    propagation_progress = i
    tic
    for j = 1:node_quantity
        if (activated_node_log(j,i-1) ~= 0) %Checks if the node is activated before checking for adjacent nodes
            [row,col] = find(brain_mesh.ConnectivityList == activated_node_log(j,i-1)); % The find function determines the row and column of the infected nodes in the Connectivity List
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
        activated_node_log(k,i) = unique_adjacent_nodes(k,1);
        node_values(unique_adjacent_nodes(:,1), i) = 1;
    end
    %Breaks if reached V5 area of cortex
    if (node_values(end_node,i) ~= 0)
        sound(success_chime, Fs);
        return;
    end
    toc
    %sound(success_chime, Fs);
end
toc

clear activated_node_log adjacent_node_list all_adjacent_nodes col i j k patient_zero row unique_adjacent_nodes unique_rows
%% Matrix to store when a node becomes activated and how many iterations have passed
%In a CSD episode after a neuron activates, its BOLD signal will drop
% exponentially. Hence this position tracker keeps a log of what iteration
% a node becomes activated.
position_tracker = node_values;
for i = 1:node_quantity
    for j = 1:iteration_number
        if (position_tracker(i,j) ~= 0)
            position_tracker(i,j+1) = position_tracker(i,j) + 1;
        end
    end
end
clear i j
%% EXPONENTIAL DECAY: Calculates the exponential decay of neuron activity after the wave has passed over.
% Combines the position tracker matrix with the activated nodes matrix to determine the value of a node over time
% Exponential decay function
time = 45; %How long the exp decay is modelled for
decay_rate = -0.1; %Negative to be a decaying model, higher to be sharper rate of decay
decay = zeros(1,time); % Creates a matrix with exp values, updates /min (45 mins)

%Creates the data for the exponential decay
for i = 6:time
    decay(1,i) = exp(decay_rate*(i-5)); 
end
decay(1,1:5) = 1; %Sets the first five minutes of the value to stay at one
decay(1,time) = 0; %Sets the final minute to be zero

%Limits the position tracker to be less than the time in the exponential
%decay data.

for i = 1:node_quantity
    for j = 1:iteration_number
        %Limits the position tracker as anything over 45 is not needed
        if (position_tracker(i,j) > time)
            position_tracker(i,j) = time;
        end
        %Sets the node value to its decayed value
        if (position_tracker(i,j) ~= 0)
            node_values(i,j) = decay(1,(position_tracker(i,j)));
        end
    end
end
check = "success"
clear i j limit tau

%% Creates and plots the brain
f = figure;

colorbar
colourMap = patch('faces', brain_mesh.ConnectivityList, 'Vertices', brain_mesh.Points);
colourMap.FaceColor = [0.53 0.55 0.64];
colourMap.EdgeColor = [0 0 0];
hold on

%% Animates the propagation of activated nodes by placing red markers on the infected nodes each iteration

for j = 1:iteration_number
    plot_progress = j
    for i = 1:node_quantity
        if (node_values(i,j) ~= 0)
            plot3(brain_mesh.Points(i,1), brain_mesh.Points(i,2), brain_mesh.Points(i,3), '.', 'markersize', 20, 'color', ([node_values(i,j) 0 0]));
            hold on
        end
    end
    pause(0.05); %Quickest it can go
end

%% Plots the brain activty as patches
% Create a colormap based on the node values

for j = 1:iteration_number
    plot_progress = j
    
    % Create a matrix of color values based on the node values
    node_colors = zeros(node_quantity, 3);
    node_colors(:,2) = 1; %Green base color
    for i = 1:node_quantity
        if (node_values(i,j) ~= 0)
            bal = node_values(i,j);
            node_colors(i,2) = 1 - bal;
            node_colors(i,1) = bal; %Changes to red when activated
        end
    end
    
    % Create a patch object and set the FaceVertexCData property to node_colors
    colourMap = patch('faces', brain_mesh.ConnectivityList, 'Vertices', brain_mesh.Points,'FaceVertexCData', node_colors, 'FaceColor', 'interp');
    colorbar
    
    % Set the axis equal, add a title, and pause
    title(sprintf('Iteration: %d', j))
    pause(0.05);
    
    % Delete the patch object to clear the plot for the next iteration
    delete(colourMap);
end
%% Plots a random set of markers to find an appropriate start point for CSD

start = round(linspace(1, 10000, 99));

clf
colourMap = patch('faces', brain_mesh.ConnectivityList, 'Vertices', brain_mesh.Points);
axis vis3d
axis equal
colourMap.FaceColor = [0 0.40 0.13];
colourMap.EdgeColor = [0 0 0];
hold on

for i = 2:size(start,2)
    text(brain_mesh.Points(start(1,i),1),brain_mesh.Points(start(1,i),2),brain_mesh.Points(start(1,i),3),num2str(start(1,i)), 'Color','r','FontSize',20);
end

for i = 1:node_quantity/2
    plot3(brain_mesh.Points(i,1), brain_mesh.Points(i,2), brain_mesh.Points(i,3), '.', 'markersize', 20);
end


%% Plots a graph of the amount of activated nodes over time
total = sum(node_values);
quantity = zeros(1,iteration_number);

for i = 1:iteration_number
    quantity(1,i) = size(find(node_values(:,i) ~= 0),1);
end

plot(1:iteration_number, quantity(1,1:iteration_number), 'LineWidth',3);
ylabel('No. of activated nodes.', 'FontSize',30);
xlabel('Model time: /min','FontSize',30);
title('Quantity of activated nodes for the duration of the simulated model.','FontSize',30);