%% AUSM
clc
clear all;
addpath("SRC\FUNCTIONS\");
addpath("SRC\FUNCTIONS\FLUX\");
addpath("SRC\FUNCTIONS\SPLIT\");

%% Propiedades fisicas

gamma = 1.4;
c = 300;

%% Discretizacion matematica
Lx = 1;
Ly = 1;
dt = 0.0001;
N_steps = 1e3;

Nx = 50; % Number of grid points in x-direction
Ny = 50; % Number of grid points in y-direction
dx = Lx / (Nx - 1); % Grid spacing in x-direction
dy = Ly / (Ny - 1); % Grid spacing in y-direction

CFL = c*dt/dx;

% Create a mesh grid for the computational domain
mesh_x = zeros(Ny, Nx);
mesh_y = zeros(Ny, Nx);

for i = 1:Nx
    mesh_x(:, i) = Lx/Nx*(i);
end

for i = 1:Ny
    mesh_y(i, :) = Ly/Ny*(i);
end

%% Initialize variables

ro_0 = 1;
u_0 = 1;
v_0 = 0.2;
P_0 = 1;
E_0 = P_0/((gamma - 1)*ro_0) + (u_0^2 + v_0^2)/2;

U(1, :, :) = ones(Ny, Nx)*ro_0;               % Densidad constante en toda la malla
U(2, :, :) = ones(Ny, Nx)*ro_0 * u_0;         % Momento x constante
U(3, :, :) = ones(Ny, Nx)*ro_0 * v_0;         % Momento y constante
U(4, :, :) = ones(Ny, Nx)*ro_0 * E_0;        % Energía total (rho * Et) constante


% Inicializacion animacion

figure
h = pcolor(mesh_x, mesh_y, zeros(size(mesh_x)));
shading flat
axis equal tight
colormap(turbo)
colorbar
caxis([-5 5])
xlabel('x')
ylabel('y')
title('Vorticidad \omega')
[Ny_loc, Nx_loc] = size(mesh_x);
vort = zeros(Ny_loc, Nx_loc);
drawnow

for n = 1:N_steps
    U = Euler_explicit(U, mesh_x, mesh_y, dt);
    
    rho = squeeze(U(1,:,:));
    u   = squeeze(U(2,:,:)) ./ rho;
    v   = squeeze(U(3,:,:)) ./ rho;
    
    vort(2:Ny_loc-1, 2:Nx_loc-1) = ...
    (v(2:Ny_loc-1, 3:Nx_loc) - v(2:Ny_loc-1, 1:Nx_loc-2))/(2*dx) ...
  - (u(3:Ny_loc, 2:Nx_loc-1) - u(1:Ny_loc-2, 2:Nx_loc-1))/(2*dy);
title(sprintf('Vorticidad \\omega  |  Paso %d', n))
drawnow

end

