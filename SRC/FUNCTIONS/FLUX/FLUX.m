function F = FLUX(U, mesh_x, mesh_y)
    
    dx = mesh_x(1,2) - mesh_x(1,1);
    dy = mesh_y(2,1) - mesh_y(1,1);

    Fc = convective(U);
    M = mach(U);
    P = pressure(U);
    F = dy*F_UP(Fc, M, P) + dy*F_DOWN(Fc, M, P) + dx*F_RIGHT(Fc, M, P) + dx*F_LEFT(Fc, M, P);
end