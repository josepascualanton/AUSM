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
dt = 0.1;
N_steps = 1e3;

Nx = 10; % Number of grid points in x-direction
Ny = 10; % Number of grid points in y-direction
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

[x, y] = meshgrid(0:dx:Lx, 0:dy:Ly);

%% Initialize variables

ro_0 = 1;
u_0 = 1;
v_0 = 0.2;
P_0 = 1;
E_0 = P_0/((gamma - 1)*ro_0) + (u_0^2 + v_0^2);

U(1, :, :) = ones(Ny, Nx)*ro_0;               % Densidad constante en toda la malla
U(2, :, :) = ones(Ny, Nx)*ro_0 * u_0;         % Momento x constante
U(3, :, :) = ones(Ny, Nx)*ro_0 * v_0;         % Momento y constante
U(4, :, :) = ones(Ny, Nx)*ro_0 * E_0;        % Energía total (rho * Et) constante


% Inicializacion animacion

figure;
h = imagesc(x(1,:), y(:,1), squeeze(U(1,:,:)));
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Densidad \rho')
caxis([0.8 1.2])   % ajusta según tu problema
drawnow


for n = 1:N_steps
    U = Euler_explicit(U, mesh_x, mesh_y, dt);

    if mod(n,5) == 0
        set(h, 'CData', squeeze(U(2,:,:)./U(1,:,:);
        title(sprintf('Densidad \\rho  |  Paso %d', n))
        drawnow
    end

end

