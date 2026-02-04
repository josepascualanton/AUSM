function F = FLUX(U, mesh_x, mesh_y)
    
    dx = mesh_x(1,2) - mesh_x(1,1);
    dy = mesh_y(2,1) - mesh_y(1,1);

    Fc = convective(U)
    [Mu, Mv] = mach(U);
    P = pressure(U);

    F = dx*F_UP(Fc, Mv, P) + dx*F_DOWN(Fc, Mv, P) + dy*F_RIGHT(Fc, Mu, P) + dy*F_LEFT(Fc, Mu, P);

end