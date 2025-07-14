%% Defines the resolution of the voxel data
resolution = 40;

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
    voxelprogress = vol
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
end
toc

%% Animates and plots the voxel data
for frame = 1:iteration_number
    plotprogress = frame
    subplot(1,3,1)
    title("x axis")
    imagesc(sum(voxel_maps(:,:,:,frame),3))
    hold on

    subplot(1,3,2)
    title("y axis")
    imagesc(squeeze(sum(voxel_maps(:,:,:,frame),2)))
    hold on

    subplot(1,3,3)
    title("z axis")
    imagesc(squeeze(sum(voxel_maps(:,:,:,frame),1)))
    hold on

    pause(0.2)
end