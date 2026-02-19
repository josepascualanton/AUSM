%% AUSM
clc
clear all;
close all;

addpath("SRC\FUNCTIONS\");
addpath("SRC\FUNCTIONS\FLUX\");
addpath("SRC\FUNCTIONS\SPLIT\");

%% Propiedades fisicas

gamma = 1.4;

%% Discretizacion matematica
Lx = 1;
Ly = 1;
dt = 0.000001;
N_steps = 5000;

Nx = 50; % Number of grid points in x-direction
Ny = 50; % Number of grid points in y-direction
dx = Lx / (Nx - 1); % Grid spacing in x-direction
dy = Ly / (Ny - 1); % Grid spacing in y-direction


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

ro_0 = 1.225;
u_0 = 5;
v_0 = 2;
P_0 = 101325;
E_0 = P_0/((gamma - 1)*ro_0) + (u_0^2 + v_0^2)/2;

epsilon = 0.1;
U(1, :, :) = ones(Ny, Nx)*ro_0;               % Densidad constante en toda la malla
U(2,:,:) = ro_0 * (u_0 + 0*epsilon*randn(Ny,Nx));
U(3,:,:) = ro_0 * (v_0 + 0*epsilon*randn(Ny,Nx));
U(4, :, :) = ones(Ny, Nx)*ro_0 * E_0;        % Energía total (rho * Et) constante


% Inicializacion animacion
fig = figure('Units','pixels','Position',[200 100 1500 800]);


tiledlayout(1, 2,'TileSpacing','compact','Padding','compact')

% --- VORTICIDAD + VELOCIDAD ---
nexttile
h1 = pcolor(mesh_x, mesh_y, zeros(size(mesh_x)));
shading flat
hold on
hQ = quiver(mesh_x, mesh_y, ...
            zeros(size(mesh_x)), ...
            zeros(size(mesh_y)), ...
            'k');
hold off
axis equal tight
colormap(turbo)
colorbar
caxis([-5 5])
xlabel('x')
ylabel('y')
title('Vorticidad + campo de velocidades')

% --- PRESIÓN ---
nexttile
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
skip = 4;
set(hQ, 'AutoScale', 'on', 'AutoScaleFactor', 1)


hIter = annotation('textbox',[0.45 0.92 0.2 0.05], ...
    'String','Iter = 0','FitBoxToText','on', ...
    'BackgroundColor','w','EdgeColor','k', ...
    'FontSize',12,'FontWeight','bold');


for n = 1:N_steps

    a0 = sqrt(gamma*P_0/ro_0);
    dt = 0.05 * min(dx,dy) / (abs(u_0) + a0);
    
    U(2, Ny/2, Nx/2) = 0;

    % Izquierda / derecha
U(:, :, 1)   = U(:, :, 2);
U(:, :, end) = U(:, :, end-1);

% Abajo / arriba
U(:, 1, :)   = U(:, 2, :);
U(:, end, :) = U(:, end-1, :);

    U = Euler_explicit(U, mesh_x, mesh_y, dt);
    
    rho = squeeze(U(1,:,:));
    u2 = squeeze(U(2,:,:));
    u3 = squeeze(U(3,:,:));
    u4 = squeeze(U(4,:,:));

    u = u2 ./ rho;
    v = u3 ./ rho;

    % --- CAMPO DE VELOCIDADES ---
% --- CAMPO DE VELOCIDADES (sobre vorticidad) ---
set(hQ, ...
    'XData', mesh_x(1:skip:end, 1:skip:end), ...
    'YData', mesh_y(1:skip:end, 1:skip:end), ...
    'UData', u(1:skip:end, 1:skip:end), ...
    'VData', v(1:skip:end, 1:skip:end));

set(hIter,'String',sprintf('Iter = %d',n));


    p = (gamma-1) .* (u4 - 0.5 .* rho .* (u.^2 + v.^2));
    p = max(p, eps);
    rho = max(rho, 1e-6);
    
    if any(rho(:) <= 0)
    [iy, ix] = find(rho <= 0, 1);
    fprintf('rho negativa en (%d,%d), iter %d\n', iy, ix, n);
    fprintf('rho = %e\n', rho(iy,ix));
    keyboard
end

    % --- VORTICIDAD ---
vort(2:Ny_loc-1, 2:Nx_loc-1) = ...
    (v(2:Ny_loc-1, 3:Nx_loc) - v(2:Ny_loc-1, 1:Nx_loc-2))/(2*dx) ...
  - (u(3:Ny_loc, 2:Nx_loc-1) - u(1:Ny_loc-2, 2:Nx_loc-1))/(2*dy);

    set(h1, 'CData', vort)
    set(h2, 'CData', p)
    
    drawnow limitrate
end

