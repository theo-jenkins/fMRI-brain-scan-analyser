%% Clear and reset varaibles
clear
close all

%% Loads our voxel maps
load saved_data/voxel_maps/147_brain_2of3.mat;
coordinate_image = imread("pictures/square.png");
[success_chime, Fs] = audioread("sound/UI_CHIME_SOUND.mp3");
voxel_maps = voxel_maps_noisy;
clear voxel_maps_noisy

%% Loads Dr Gallichans functions
%addpath E:\Users\Theo Jenkins\Documents\University\Engineering\YEAR 3\PROJECT\PERSONAL MATLAB TESTING\WAVE_CHECK_SANDBOX/utils;

%% Sets a grid to be full of data
x_size = size(voxel_maps,1);
y_size = size(voxel_maps,2);
z_size = size(voxel_maps,3);
timesteps = size(voxel_maps,4);

%% Distributes markers across the 2D plane
% Set the number of markers PER PLANE, so true amount will be cubed (num_markers^3)!
num_markers = 20;

% Sets the border for where the markers can be placed (-5% from the edge)
cap = 0.05;
x_max_limit = x_size - (cap * x_size); x_min_limit = (cap * x_size);
x_limit = linspace(x_min_limit, x_max_limit, num_markers);

y_max_limit = y_size - (cap * y_size); y_min_limit = (cap * y_size);
y_limit = linspace(y_min_limit, y_max_limit, num_markers);

z_max_limit = z_size - (cap * z_size); z_min_limit = (cap * z_size);
z_limit = linspace(z_min_limit, z_max_limit, num_markers);

% Sets the cartesian coordinates for the markers
[Ym, Xm, Zm] = ndgrid(y_limit, x_limit, z_limit);
marker_xyz = zeros(3,(num_markers^3), timesteps);

for i = 1:(num_markers^3)
    [y_id, x_id, z_id] = ind2sub(size(Xm),i);
    marker_xyz(1,i,1) = Xm(x_id, y_id, z_id);
    marker_xyz(2,i,1) = Ym(x_id, y_id, z_id);
    marker_xyz(3,i,1) = Zm(x_id, y_id, z_id);
end

clear x_max_limit x_min_limit x_limit y_max_limit y_min_limit y_limit z_max_limit z_min_limit z_limit x_id y_id z_id
%% Gradient of the markers:
%Creates a matrix to store all the [x,y,z] gradients for the markers over
%time
gradients_xyz = zeros(3, num_markers, timesteps);

clear smooth_factor Kernel

%% Moves each marker in the direction of the local minima

%Sets a 3D matrix of the grid point cartesian points
[X, Y, Z] = meshgrid(1:x_size, 1:y_size, 1:z_size);

for i = 2:timesteps
    timestep_progress = i
    for j = 1:(num_markers^3)
        %Loads the cartesian coordinates of a marker
        current_X = marker_xyz(1, j, (i-1));
        current_Y = marker_xyz(2, j, (i-1));
        current_Z = marker_xyz(3, j, (i-1));

        %Calculates the gradient matrices in the x/y/z directions
        [gx, gy, gz] = gradient(voxel_maps(:,:,:,i));

        %Estimates values for the gradient at the marker using 3D
        %interpolation
        gradient_X = interp3(X,Y,Z,gx,current_X,current_Y,current_Z);
        gradient_Y = interp3(X,Y,Z,gy,current_X,current_Y,current_Z);
        gradient_Z = interp3(X,Y,Z,gz,current_X,current_Y,current_Z);

        %Determines the direction for the marker to move in
        dir_X = -(gradient_X);
        dir_Y = -(gradient_Y);
        dir_Z = -(gradient_Z);

        %Calculates the markers new position
        step_size = 2.5;
        new_X = current_X + (step_size * dir_X);
        new_Y = current_Y + (step_size * dir_Y);
        new_Z = current_Z + (step_size * dir_Z);

        %Loads the new location into the matrix
        marker_xyz(1,j,i) = new_X;
        marker_xyz(2,j,i) = new_Y;
        marker_xyz(3,j,i) = new_Z;
    end
    sound(success_chime, Fs);
end

clear current_X current_Y current_Z gx gy gz gradient_X gradient_Y gradient_Z dir_X dir_Y dir_Z step_size new_X new_Y new_Z x_size y_size z_size timestep_progress

%% Calculates the total displacement of the markers
%Uses a for loop to determine the total distance moved from the first to
%the last voxel for every marker

% Calculates the displacement for marker_xyz
marker_displacement = zeros(num_markers^3, 1);
for i = 1:(num_markers^3)
    displacement = sqrt(((marker_xyz(1,i,timesteps) - marker_xyz(1,i,1))^2) + ((marker_xyz(2,i,timesteps) - marker_xyz(2,i,1))^2) + ((marker_xyz(3,i,timesteps) - marker_xyz(3,i,1))^2));
    marker_displacement(i,1) = displacement;
end

%Checks for NaN values and sets to zero
marker_displacement(isnan(marker_displacement)) = 0;

%Sorts the markers into order of greatest displacement, into the matrix
%'chosen_markers'
[~, chosen_markers] = sort(marker_displacement, 'descend');
marker_limit =1000;

%% Plots histograms for patient_markers and simulated_markers
% Calculates the displacement for simulated/patient_markers
simulated_displacement = zeros(num_markers^2,1);
patient_displacement = zeros(num_markers^2,1);
for i = 1:(num_markers^3)
    simulated_d = sqrt(((simulated_markers(1,i,timesteps) - simulated_markers(1,i,1))^2) + ((simulated_markers(2,i,timesteps) - simulated_markers(2,i,1))^2) + ((simulated_markers(3,i,timesteps) - simulated_markers(3,i,1))^2));
    patient_d = sqrt(((patient_markers(1,i,timesteps) - patient_markers(1,i,1))^2) + ((patient_markers(2,i,timesteps) - patient_markers(2,i,1))^2) + ((patient_markers(3,i,timesteps) - patient_markers(3,i,1))^2));
    simulated_displacement(i,1) = simulated_d;
    patient_displacement(i,1) = patient_d;
end

% Sets the data between 0 - 1 and removes 'NaN'
normalised_patient = normalize(patient_displacement, 'range', [0 1]);
normalised_patient(isnan(normalised_patient)) = 0;
normalised_simulated = normalize(simulated_displacement, 'range', [0 1]);
normalised_simulated(isnan(normalised_simulated)) = 0;

% Plots histograms
figure;
%h_p = histogram(patient,100);
h_p = histogram(normalize(patient, 'range', [0 1]),100);
h_p.FaceColor = [0.8500 0.3250 0.0980]
xlabel('Displacement (mm)', 'fontsize', 25)
ylabel('Frequency of markers', 'FontSize', 25)

figure;
%h_e = histogram(empty, 100);
h_e = histogram(normalize(empty, 'range', [0 1]),100);
h_e.FaceColor = [0.6350 0.0780 0.1840]
xlabel('Displacement (mm)', 'fontsize', 25)
ylabel('Frequency of markers', 'FontSize', 25)

figure;
%h_s = histogram(simulated,100);
h_s = histogram(normalize(simulated, 'range', [0 1]),100);
xlabel('Displacement (mm)', 'fontsize', 25)
ylabel('Frequency of markers', 'FontSize', 25)


%% Graphs the displacement of the markers

figure;
hold on;
h1 = histogram(patient, 80);
figure;
h2 = histogram(simulated, 80);
h = histogram(marker_displacement(1,:), 8000);
sum = zeros(1,8000);
for i = 2 :size(marker_displacement,1)
    sum(1,i) = sum(1,i-1) + marker_displacement(i-1,1);
end
b = bar(1:num_markers^3, sum(1,:), 1);
b = bar(1:num_markers^3, marker_displacement(1:num_markers^3,1), 1);
xlabel('Index of marker', 'FontSize',25);
ylabel('Displacement (mm)', 'FontSize',25);
b.FaceColor = 'flat';

for i = 1:length(marker_displacement)
    progress = i;
    if marker_displacement(i,1) < 1 %Lower limit for change of color
        color = [0.8 0 0];%[0.98 0.13 0.08];
    else
        color = [0 0.1 0.8];%[0.2 0.35 0.9];
    end
    b.CData(i,:) = color;
end

hold off;

%% Plots the markers which have moved the most

f = figure;
for i = 1:timesteps
    plot_progress = i
    plot3(marker_xyz(1,chosen_markers(1:marker_limit),i), marker_xyz(2,chosen_markers(1:marker_limit),i), marker_xyz(3,chosen_markers(1:marker_limit),i), 'x', 'Color', 'r');
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis ([0 100 0 100 0 100]);
    pause(0.1);
end

%% Plots all markers no matter if they moved or not
f = figure;
for i = 1:timesteps
    plot_progress = i
    plot3(marker_xyz(1,:,i), marker_xyz(2,:,i), marker_xyz(3,:,i), 'x', 'Color', 'r');
    pause(0.2);
end

%% Plots the voxels of the brain scan with the markers movement overlayed
f = figure;
max_c = max(voxel_maps(:));
min_c = min(voxel_maps(:));
c_limit = 15;

%Create subplot 1: (top down view)
subplot(2,3,4);
v1 = imagesc(NaN);
hold on
m1 = plot(NaN, NaN, '.', 'Color', 'w', 'LineWidth',10,'markersize',10);
title('Coronal View', 'FontSize',15);
xlabel('Z','FontSize',20);
ylabel('Y','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight
caxis([(c_limit*min_c) (c_limit*max_c)])

%Create subplot 1 (back view):
subplot(2,3,5);
v2 = imagesc(NaN);
hold on
m2 = plot(NaN, NaN, '.', 'Color', 'w', 'LineWidth',10,'markersize',10);
title('Horizontal View', 'FontSize',15);
xlabel('Z','FontSize',20);
ylabel('X','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight
caxis([(c_limit*min_c) (c_limit*max_c)])

%Create subplot 3 (side view [infection side]):
subplot(2,3,6);
v3 = imagesc(NaN);
hold on
m3 = plot(NaN, NaN, '.', 'Color', 'w', 'LineWidth',10,'markersize',10);
title('Sagittal View', 'FontSize', 15);
xlabel('X','FontSize',20);
ylabel('Y','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight
caxis([(c_limit*min_c) (c_limit*max_c)])

for i = 1:timesteps
    %Progress check
    plot_progress = i

    %Plots the ZY plane: A view of the voxels (Xv) along the X axis, 
    subplot(2,3,4);
    normalised_x = squeeze(sum(voxel_maps(:,:,:,i),2));
    normalised_x = flip(rot90(normalised_x));
    set(v1, 'CData', normalised_x);
    %Markers:
    Xm_view = marker_xyz([2 3], chosen_markers(1:marker_limit), i);
    m1.XData = Xm_view(1,:); m1.YData = Xm_view(2,:);

    %Plots the XZ plane: A view of the voxels (Yv) along the Y axis,
    %overlayed with the movement of the markers (Ym)
    subplot(2,3,5);
    normalised_y = squeeze(sum(voxel_maps(:,:,:,i),3));
    normalised_y = flip(rot90(normalised_y));
    set(v2, 'CData', normalised_y); 
    %Markers:
    Ym_view = marker_xyz([1 2], chosen_markers(1:marker_limit), i);
    m2.XData = Ym_view(2,:); m2.YData = Ym_view(1,:);

    %Plots the XY plane: A view of the voxels (Zv) along the Z axis,
    %overlayed with the movement of the markers (Zm)
    subplot(2,3,6);
    normalised_z = squeeze(sum(voxel_maps(:,:,:,i),1));
    normalised_z = rot90((normalised_z),3);
    set(v3, 'CData', normalised_z);
    %Markers:
    Zm_view = 30 - (marker_xyz([1 3], chosen_markers(1:marker_limit), i));
    m3.XData = Zm_view(1,:); m3.YData = Zm_view(2,:);
    %Pause
    
    %Creates subplot 4 (coordinate key):
    subplot(2,3,[1,2,3]);
    title('Volume: ');
    axis equal
    plot(i,1,'o','Color','b');
    axis([0 timesteps 0 2])
    hold on

    %Pause
    pause(0.1);
end

clear displacement i j plot_progress X Xm Xm_view Xv_view Y Ym Ym_view Yv_view Z Zm Zm_view Zv_view corrected_Xv_view corrected_Yv_view corrected_Zm_view corrected_Zv_view

%% Plots the voxels of the brain: WITH NO MARKERS
f = figure;

%Create subplot 1: (top down view)
subplot(2,3,4);
v1 = imagesc(NaN);
hold on
m1 = plot(NaN, NaN, 'x', 'Color', 'w');
title('Coronal View', 'FontSize',15);
xlabel('Z','FontSize',20);
ylabel('Y','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight

%Create subplot 1 (back view):
subplot(2,3,5);
v2 = imagesc(NaN);
hold on
m2 = plot(NaN, NaN, 'x', 'Color', 'w');
title('Horizontal View', 'FontSize',15);
xlabel('Z','FontSize',20);
ylabel('X','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight

%Create subplot 3 (side view [infection side]):
subplot(2,3,6);
v3 = imagesc(NaN);
hold on
m3 = plot(NaN, NaN, 'x', 'Color', 'r');
title('Sagittal View', 'FontSize', 15);
xlabel('X','FontSize',20);
ylabel('Y','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight

for i = 1:timesteps
    %Progress check
    plot_progress = i

    %Plots the ZY plane: A view of the voxels (Xv) along the X axis, 
    subplot(2,3,4);
    normalised_x = squeeze(sum(voxel_maps(:,:,:,i),2));
    normalised_x = flip(rot90(normalised_x));
    set(v1, 'CData', normalised_x); 

    %Plots the XZ plane: A view of the voxels (Yv) along the Y axis,
    %overlayed with the movement of the markers (Ym)
    subplot(2,3,5);
    normalised_y = squeeze(sum(voxel_maps(:,:,:,i),3));
    normalised_y = flip(rot90(normalised_y));
    set(v2, 'CData', normalised_y); 

    %Plots the XY plane: A view of the voxels (Zv) along the Z axis,
    %overlayed with the movement of the markers (Zm)
    subplot(2,3,6);
    normalised_z = squeeze(sum(voxel_maps(:,:,:,i),1));
    normalised_z = rot90((normalised_z),3);
    set(v3, 'CData', normalised_z); 
    %Pause
    
    %Creates subplot 4 (coordinate key):
    subplot(2,3,[1,2,3]);
    title('Volume: ');
    axis equal
    plot(i,1,'o','Color','b');
    axis([0 timesteps 0 2])
    hold on

    %Pause
    pause(0.1);
end