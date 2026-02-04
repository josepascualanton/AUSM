function U = Euler_explicit(U, mesh_x, mesh_y, dt)
    

    U = U + dt*FLUX(U, mesh_x, mesh_y);
end