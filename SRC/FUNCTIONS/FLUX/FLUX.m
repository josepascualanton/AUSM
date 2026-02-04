function F = FLUX(U, mesh_x, mesh_y)
    
    Fc = convective(U);
    M = mach(U);
    P = pressure(U);
    F = F_UP(Fc, M, P) + F_DOWN(Fc, M, P) + F_RIGHT(Fc, M, P) + F_LEFT(Fc, M, P);
    F = F*0.01;
end