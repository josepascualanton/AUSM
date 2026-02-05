function F = FLUX(U, mesh_x, mesh_y)
    
    dx = mesh_x(1,2) - mesh_x(1,1);
    dy = mesh_y(2,1) - mesh_y(1,1);

    [Fcx, Fcy] = convective(U);
   
    [Mu, Mv] = mach(U);
    P = pressure(U);

    F = dx*F_UP(Fcy, Mv, P) + dx*F_DOWN(Fcy, Mv, P) + dy*F_RIGHT(Fcx, Mu, P) + dy*F_LEFT(Fcx, Mu, P);

end