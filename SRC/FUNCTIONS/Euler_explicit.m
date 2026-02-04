function U = Euler_explicit(U, mesh_x, mesh_y, dt)
    
    dx = mesh_x(1,2) - mesh_x(1,1);
    dy = mesh_y(2,1) - mesh_y(1,1);

    vol = dx * dy;
    U = U - dt*FLUX(U, mesh_x, mesh_y)/vol;
end