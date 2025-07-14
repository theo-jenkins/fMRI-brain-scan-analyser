%% Clear and reset varaibles
clear
close all

%% Loads 2D Image
load('diss_algorithm\perlin_image.mat');
image = perlin_image;

%% Sets size
x_size = size(image,1);
y_size = size(image,2);
timesteps = 100;

%% Smooths the voxel data

%Smooths the data
smooth_factor = 5;
Kernel = 1/(smooth_factor^2) * ones(smooth_factor);
image = convn(image, Kernel, 'same');

%% Distributes markers across the 2D plane
% Set the number of markers PER PLANE, so true amount will be cubed (num_markers^3)!
num_markers = 14;
cap = 0.05;

% Sets the border for where the markers can be placed (-15% from the edge)
x_max_limit = x_size - (cap * x_size); x_min_limit = (cap * x_size);
x_limit = linspace(x_min_limit, x_max_limit, num_markers);

y_max_limit = y_size - (cap * y_size); y_min_limit = (cap * y_size);
y_limit = linspace(y_min_limit, y_max_limit, num_markers);

% Sets the cartesian coordinates for the markers
[Ym, Xm] = ndgrid(y_limit, x_limit);
marker_xy = zeros(2,(num_markers^2), timesteps);

for i = 1:(num_markers^2)
    [y_id, x_id] = ind2sub(size(Xm),i);
    marker_xy(1,i,1) = Xm(x_id, y_id);
    marker_xy(2,i,1) = Ym(x_id, y_id);
end

%% update marker positions

%Creates a matrix to store all the [x,y] gradients for the markers over
%time
gradients_of_markers = zeros(2, num_markers^2, timesteps);

%% Moves each marker in the direction of the local minima
for i = 2:timesteps
    for j = 1:num_markers^2
        %Loads the cartesian coordinates of a marker
        current_X = marker_xy(1,j,i-1);
        current_Y = marker_xy(2,j,i-1);
    
        %Loads the gradient matrices for the x and y direction 
        [gx, gy] = gradient(image);
    
        %Creates a cartesian coordinate grid for the interpolation
        [X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
    
        %Loads a value for the gradient at the coordinates of the marker
        grad_X_at_current = interp2(X, Y, gx, current_X, current_Y);
        grad_Y_at_current = interp2(X, Y, gy, current_X, current_Y);
    
        %Sets the intended gradient of the marker (negative as finding local
        %minima)
        dir_X = (grad_X_at_current);
        dir_Y = (grad_Y_at_current);
    
        %Calculates markers new position
        step_size = 0.4;
        new_X = current_X + (step_size * dir_X);
        new_Y = current_Y + (step_size * dir_Y);

        %Keeps the markers jumping too far
        cap = 0.05;
        if abs(current_X - new_X) > (cap * x_size)
            new_X = current_X + (cap * x_size);
        end
        if abs(current_Y - new_Y) > (cap * y_size)
            new_Y = current_Y + (cap * y_size);
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

    
        %Loads the new location into the intended matrix
        marker_xy(1,j,i) = new_X;
        marker_xy(2,j,i) = new_Y;
    end
end

%% Plots the animation of the volume
clims = [-20 20];

for i = 1:timesteps
    output = i
    %Plots the next timestep
    imagesc(image, clims);
    hold on
    %set(sandbox, 'CData', sandbox_with_time(:,:,i));

    %For loop to plot each marker and there changing position
    imagesc(image);
    hold on
    for j = 1:num_markers^2
        plot(marker_xy(1,j,i), marker_xy(2,j,i), '.', 'Color', 'w', 'LineWidth',5,'MarkerSize',20);
        hold on
    end

    pause(0.1);
end

