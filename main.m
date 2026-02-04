%% AUSM

Lx = 1;
Ly = 1;
dt = 0.01;

Nx = 100; % Number of grid points in x-direction
Ny = 100; % Number of grid points in y-direction
dx = Lx / (Nx - 1); % Grid spacing in x-direction
dy = Ly / (Ny - 1); % Grid spacing in y-direction

MESH = zeros(Ny, Nx);