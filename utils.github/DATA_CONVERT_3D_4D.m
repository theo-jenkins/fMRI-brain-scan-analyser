%% Fresh start
clear
close all

%% Load saved data
brain_mesh = stlread("stl/dans_brain.stl");
%brain_mesh = stlread("utah_teapot.stl");

[success_chime, Fs] = audioread("sound/UI_CHIME_SOUND.mp3");

node_quantity = size(node_values(:,1),1);
iteration_number = size(node_values(1,:),2);

%% Defines the resolution of the voxel data
resolution = 30;

Nx = resolution;
Ny = resolution;
Nz = resolution;

%% Sets a 3 coordinate system of voxels surrounding the brain mesh
xrange = linspace(min(brain_mesh.Points(:,1)), max(brain_mesh.Points(:,1)), Nx);
yrange = linspace(min(brain_mesh.Points(:,2)), max(brain_mesh.Points(:,2)), Ny);
zrange = linspace(min(brain_mesh.Points(:,3)), max(brain_mesh.Points(:,3)), Nz);

%% Creates a 4D matrix to be filled with volumes
voxel_maps = zeros(Nx,Ny,Nz, iteration_number);

%% For each volume/iteration, the program determines what voxel it inhabits
tic
for vol = 1:iteration_number
    voxel_progress = vol
    for node = 1:node_quantity
        if (node_values(node,vol) ~= 0)
            nx = brain_mesh.Points(node,1);
            ny = brain_mesh.Points(node,2);
            nz = brain_mesh.Points(node,3);

            voxel_x = find(xrange>nx,1);
            voxel_y = find(yrange>ny,1);
            voxel_z = find(zrange>nz,1);

            voxel_maps(voxel_x, voxel_y, voxel_z, vol) = voxel_maps(voxel_x, voxel_y, voxel_z, vol) + node_values(node,vol);

            
        end
    end
    sound(success_chime, Fs);
end
toc
clear node nx Nx ny Ny nz Nz resolution vol voxel_x voxel_y voxel_z xrange yrange zrange voxel_progress

%% Adds noise to the voxel data
noise_scaling = 0.15;
for i = 1:iteration_number
    noise = noise_scaling * randn(size(voxel_maps,1), size(voxel_maps,2), size(voxel_maps,3));
    voxel_maps_noisy(:,:,:,i) = voxel_maps(:,:,:,i) + noise;
end
%% Smooths the voxel data

%Sets the parameters for the Gaussian filter
smooth_factor = 2;
sigma = 2;
kernel_size = ceil(sigma * smooth_factor) + 1;

%Creates and applies the gaussian filter kernel
gaussian_kernel = fspecial3('gaussian', kernel_size, sigma);
voxel_maps_noisy = imfilter(voxel_maps_noisy, gaussian_kernel, 'symmetric');

%% Plots the voxels of the brain: NO NOISE
f = figure;

%Create subplot 1: (top down view)
subplot(2,3,4);
v1 = imagesc(NaN);
hold on
title('Coronal View', 'FontSize',15);
xlabel('Z','FontSize',20);
ylabel('Y','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight

%Create subplot 1 (back view):
subplot(2,3,5);
v2 = imagesc(NaN);
hold on
title('Horizontal View', 'FontSize',15);
xlabel('Z','FontSize',20);
ylabel('X','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight

%Create subplot 3 (side view [infection side]):
subplot(2,3,6);
v3 = imagesc(NaN);
hold on
title('Sagittal View', 'FontSize', 15);
xlabel('X','FontSize',20);
ylabel('Y','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight

for i = 1:iteration_number
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
    axis([0 iteration_number 0 2])
    hold on

    %Pause
    pause(0.1);
end
%% Plots the voxel_maps_noisy data
f = figure;

%Create subplot 1: (top down view)
subplot(2,3,4);
v1 = imagesc(NaN);
hold on
title('Coronal View', 'FontSize',15);
xlabel('Z','FontSize',20);
ylabel('Y','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight

%Create subplot 1 (back view):
subplot(2,3,5);
v2 = imagesc(NaN);
hold on

title('Horizontal View', 'FontSize',15);
xlabel('Z','FontSize',20);
ylabel('X','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight

%Create subplot 3 (side view [infection side]):
subplot(2,3,6);
v3 = imagesc(NaN);
hold on
title('Sagittal View', 'FontSize', 15);
xlabel('X','FontSize',20);
ylabel('Y','FontSize',20);
set(gca, 'YDir', 'normal');
axis equal tight


for i = 1:iteration_number
    %Progress check
    plot_progress = i

    %Plots the ZY plane: A view of the voxels (Xv) along the X axis, 
    subplot(2,3,4);
    normalised_x = squeeze(sum(voxel_maps_noisy(:,:,:,i),2));
    %normalised_x = squeeze(sum(voxel_maps(:,:,:,i),2));
    normalised_x = flip(rot90(normalised_x));
    set(v1, 'CData', normalised_x); 

    %Plots the XZ plane: A view of the voxels (Yv) along the Y axis,
    %overlayed with the movement of the markers (Ym)
    subplot(2,3,5);
    normalised_y = squeeze(sum(voxel_maps_noisy(:,:,:,i),3));
    %normalised_y = squeeze(sum(voxel_maps(:,:,:,i),3));
    normalised_y = flip(rot90(normalised_y));
    set(v2, 'CData', normalised_y); 

    %Plots the XY plane: A view of the voxels (Zv) along the Z axis,
    %overlayed with the movement of the markers (Zm)
    subplot(2,3,6);
    normalised_z = squeeze(sum(voxel_maps_noisy(:,:,:,i),1));
    %normalised_z = squeeze(sum(voxel_maps(:,:,:,i),1));
    normalised_z = rot90((normalised_z),3);
    set(v3, 'CData', normalised_z); 
    %Pause
    
    %Creates subplot 4 (coordinate key):
    subplot(2,3,[1,2,3]);
    title('Volume: ');
    axis equal
    plot(i,1,'o','Color','b');
    axis([0 iteration_number 0 2])
    hold on

    %Pause
    pause(0.1);
end