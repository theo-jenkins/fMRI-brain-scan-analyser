%% Clear and reset varaibles
clear
close all

%% 'Patient names' available to load into voxel viewer
%Eleven subjects from eleven datasets
names = {'alex','brian','chris','dug','ella','finn','george','harry','isla','josh','kirsty'};

%% Loads our voxel map
patient = niftiread("openneuro_data\ella.nii.gz");
voxel_maps = double(patient);
voxel_maps = voxel_maps - mean(voxel_maps,4);

[success_chime, Fs] = audioread("sound/UI_CHIME_SOUND.mp3");
%% Smooths the voxel data

%Sets the parameters for the Gaussian filter
smooth_factor = 2;
sigma = 1; %std of gf
kernel_size = ceil(sigma * smooth_factor)+1;

%Creates and applies the gaussian filter kernel
gaussian_kernel = fspecial3('gaussian', kernel_size, sigma);
voxel_maps = imfilter(voxel_maps, gaussian_kernel, 'symmetric');
%% Sets a grid to be full of data
x_size = size(voxel_maps,1);
y_size = size(voxel_maps,2);
z_size = size(voxel_maps,3);
timesteps = size(voxel_maps,4);
timesteps = 110;

%% Distributes markers across the 2D plane
% Set the number of markers PER PLANE, so true amount will be cubed (num_markers^3)!
num_markers = 20;

%Marker offset from voxel boundary
f = 0.15; %f<0.5

% Sets the border for where the markers can be placed
x_max_limit = x_size - (f * x_size); x_min_limit = (f * x_size);
x_limit = linspace(x_min_limit, x_max_limit, num_markers);

y_max_limit = y_size - (f * y_size); y_min_limit = (f * y_size);
y_limit = linspace(y_min_limit, y_max_limit, num_markers);

z_max_limit = z_size - (f * z_size); z_min_limit = (f * z_size);
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
        [gx, gy, gz] = gradient(voxel_maps(:,:,:,1));

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
        step_size = 0.25;
        new_X = current_X + (step_size * dir_X);
        new_Y = current_Y + (step_size * dir_Y);
        new_Z = current_Z + (step_size * dir_Z);

        %Keeps the markers jumping too far
        cap = 0.05;
        if abs(current_X - new_X) > (cap * x_size)
            new_X = current_X + (cap * x_size);
        end
        if abs(current_Y - new_Y) > (cap * y_size)
            new_Y = current_Y + (cap * y_size);
        end
        if abs(current_Z - new_Z) > (cap * z_size)
            new_Z = current_Z + (cap * z_size);
        end

        %Keeps the markers in the bounds of the box
        if new_X > x_size
            new_X = x_size - step_size;
        end
        if new_X < 0
            new_X = 0 + step_size;
        end
        
        if new_Y > y_size
            new_Y = y_size - step_size;
        end
        if new_Y < 0
            new_Y = 0 + step_size;
        end
        
        if new_Z > z_size
            new_Z = z_size - step_size;
        end
        if new_Z < 0
            new_Z = 0 + step_size;
        end

        %Loads the new location into the matrix
        marker_xyz(1,j,i) = new_X;
        marker_xyz(2,j,i) = new_Y;
        marker_xyz(3,j,i) = new_Z;
    end
    sound(success_chime, Fs);
end

clear current_X current_Y current_Z gx gy gz gradient_X gradient_Y gradient_Z dir_X dir_Y dir_Z step_size new_X new_Y new_Z timestep_progress

%% Calculates the total displacement of the markers
%Uses a for loop to determine the total distance moved from the first to
%the last voxel for every marker
marker_displacement = zeros(1,(num_markers^3));
marker_limit = 500;
for i = 1:(num_markers^3)
    displacement = sqrt(((marker_xyz(1,i,timesteps) - marker_xyz(1,i,1))^2) + ((marker_xyz(2,i,timesteps) - marker_xyz(2,i,1))^2) + ((marker_xyz(3,i,timesteps) - marker_xyz(3,i,1))^2));
    marker_displacement(1,i) = displacement;
end

%Checks for NaN values and sets to zero
marker_displacement(isnan(marker_displacement)) = 0;

%Sorts the markers into order of greatest displacement, into the matrix
%'chosen_markers'
[~, chosen_markers] = sort(marker_displacement, 'descend');
if marker_limit > (num_markers^3)
    marker_limit = num_markers^3;
end
plot_markers = chosen_markers(1:marker_limit);

%% Plots the voxels of the brain: WITH NO MARKERS
f = figure;

%Create subplot 1: (top down view)
subplot(2,3,4);
v1 = imagesc(NaN);
hold on
m1 = plot(NaN, NaN, 'x', 'Color', 'w');
title('Coronal View', 'FontSize',15);
xlabel('Z');
ylabel('Y');
set(gca, 'YDir', 'normal');
axis equal tight

%Create subplot 1 (back view):
subplot(2,3,5);
v2 = imagesc(NaN);
hold on
m2 = plot(NaN, NaN, 'x', 'Color', 'w');
title('Horizontal View', 'FontSize',15);
xlabel('Z');
ylabel('X');
set(gca, 'YDir', 'normal');
axis equal tight

%Create subplot 3 (side view [infection side]):
subplot(2,3,6);
v3 = imagesc(NaN);
hold on
m3 = plot(NaN, NaN, 'x', 'Color', 'r');
title('Sagittal View', 'FontSize', 15);
xlabel('X');
ylabel('Y');
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

%% Plots the voxels of the brain scan with the markers movement overlayed

f = figure;
%Create subplot 1: (top down view)
subplot(2,3,4);
v1 = imagesc(NaN);
colorbar
hold on
m1 = plot(NaN, NaN, '.', 'Color', 'w', 'LineWidth',10,'markersize',10);
title('Coronal View', 'FontSize',15);
xlabel('Z','FontSize',20);
ylabel('Y','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight

%Create subplot 2 (back view):
subplot(2,3,5);
v2 = imagesc(NaN);
colorbar
hold on
m2 = plot(NaN, NaN, '.', 'Color', 'w', 'LineWidth',10,'markersize',10);
title('Horizontal View', 'FontSize',15);
xlabel('Z','FontSize',20);
ylabel('X','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight

%Create subplot 3 (side view [infection side]):
subplot(2,3,6);
v3 = imagesc(NaN);
colorbar
hold on
m3 = plot(NaN, NaN, '.', 'Color', 'w', 'LineWidth',10,'markersize',10);
title('Sagittal View', 'FontSize', 15);
xlabel('X','FontSize',20);
ylabel('Y','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal

%Creates subplot 4 (coordinate key):
subplot(2,3,[1,3])
title('Volume:');
axis([0 timesteps 0 2]);

for i = 1:timesteps
    %Progress check
    plot_progress = i

    %Plots the ZY plane: A view of the voxels (Xv) along the X axis, 
    subplot(2,3,4);
    Xv_view = squeeze(sum(voxel_maps(:,:,:,i),2));
    Xv_view = flip(rot90(Xv_view));
    Xm_view = marker_xyz([2 3], plot_markers, i);
    set(v1, 'CData', Xv_view); 
    m1.XData = Xm_view(1,:); m1.YData = Xm_view(2,:);

    %Plots the XZ plane: A view of the voxels (Yv) along the Y axis,
    %overlayed with the movement of the markers (Ym)
    subplot(2,3,5);
    Yv_view = squeeze(sum(voxel_maps(:,:,:,i),3));
    Yv_view = flip(rot90(Yv_view));
    Ym_view = marker_xyz([1 2], plot_markers, i);
    set(v2, 'CData', Yv_view);
    m2.XData = Ym_view(2,:); m2.YData = Ym_view(1,:);

    %Plots the XY plane: A view of the voxels (Zv) along the Z axis,
    %overlayed with the movement of the markers (Zm)
    subplot(2,3,6);
    axis tight
    Zv_view = squeeze(sum(voxel_maps(:,:,:,i,1)));
    Zv_view = rot90(Zv_view,3);
    Zm_view = (marker_xyz([1 3], plot_markers, i));
    set(v3, 'CData', Zv_view);
    m3.XData = Zm_view(1,:); m3.YData = Zm_view(2,:);

    %Plot 4:
    subplot(2,3,[1,3])
    plot(i, 1, 'o', 'color', 'b');
    hold on
    %Pause
    pause(0.3);
end

clear displacement i j plot_progress X Xm Xm_view Xv_view Y Ym Ym_view Yv_view Z Zm Zm_view Zv_view corrected_Xv_view corrected_Yv_view corrected_Zm_view corrected_Zv_view