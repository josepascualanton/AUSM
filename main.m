%% AUSM
clc
clear all;
close all;

addpath("SRC\FUNCTIONS\");
addpath("SRC\FUNCTIONS\FLUX\");
addpath("SRC\FUNCTIONS\SPLIT\");

%% Propiedades fisicas

gamma = 1.4;
c = 300;

%% Discretizacion matematica
Lx = 1;
Ly = 1;
dt = 0.000001;
N_steps = 5000;

Nx = 100; % Number of grid points in x-direction
Ny = 100; % Number of grid points in y-direction
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
u_0 = 10;
v_0 = 20;
P_0 = 100;
E_0 = P_0/((gamma - 1)*ro_0) + (u_0^2 + v_0^2)/2;

epsilon = 0.1;
U(1, :, :) = ones(Ny, Nx)*ro_0;               % Densidad constante en toda la malla
U(2,:,:) = ro_0 * (u_0 + epsilon*randn(Ny,Nx));
U(3,:,:) = ro_0 * (v_0 + epsilon*randn(Ny,Nx));
U(4, :, :) = ones(Ny, Nx)*ro_0 * E_0;        % Energía total (rho * Et) constante


% Inicializacion animacion
fig = figure;


subplot(1,2,1)
h1 = pcolor(mesh_x, mesh_y, zeros(size(mesh_x)));
shading flat
axis equal tight
colormap(turbo)
colorbar
caxis([-5 5])
xlabel('x')
ylabel('y')
title('Vorticidad \omega')

subplot(1,2,2)
h2 = pcolor(mesh_x, mesh_y, zeros(size(mesh_x)));
shading flat
axis equal tight
colormap(turbo)
colorbar
xlabel('x')
ylabel('y')
title('Presión p')

[Ny_loc, Nx_loc] = size(mesh_x);
vort = zeros(Ny_loc, Nx_loc);
drawnow


for n = 1:N_steps
    U = Euler_explicit(U, mesh_x, mesh_y, dt);
    
    rho = squeeze(U(1,:,:));
    u2 = squeeze(U(2,:,:));
    u3 = squeeze(U(3,:,:));
    u4 = squeeze(U(4,:,:));

    u = u2 ./ rho;
    v = u3 ./ rho;

    p = (gamma-1) .* (u4 - 0.5 .* rho .* (u.^2 + v.^2));
    p = max(p, eps);
    
    % --- VORTICIDAD ---
vort(2:Ny_loc-1, 2:Nx_loc-1) = ...
    (v(2:Ny_loc-1, 3:Nx_loc) - v(2:Ny_loc-1, 1:Nx_loc-2))/(2*dx) ...
  - (u(3:Ny_loc, 2:Nx_loc-1) - u(1:Ny_loc-2, 2:Nx_loc-1))/(2*dy);

    set(h1, 'CData', vort)
    set(h2, 'CData', p)
    
    drawnow limitrate
end

