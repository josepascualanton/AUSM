function [Fcx, Fcy] = convective(U)
    [~, Ny, Nx] = size(U);
    gamma = 1.4;
    Fcx = zeros(4, Ny, Nx);
    Fcy = zeros(4, Ny, Nx);

    rho(:, :) = squeeze(U(1,:,:));
    u2(:, :) = squeeze(U(2, :, :));
    u3(:, :) = squeeze(U(3, :, :));
    u4(:, :) = squeeze(U(4, :, :));

    u = u2 ./ rho;
    v = u3 ./ rho;
    
    p   = (gamma-1) .* (u4 - 0.5 .* rho .* (u.^2 + v.^2));
    p = max(p, 0.01);

    a   = sqrt(gamma .* p ./ rho);
    H   = (u4 + p) ./ rho;
    
    Fcx(1, :, :) = rho.*a;
    Fcx(2, :, :) = rho.*a.*u;
    Fcx(3, :, :) = rho.*a.*v;
    Fcx(4, :, :) = rho.*a.*H;

    Fcy(1, :, :) = rho.*a;
    Fcy(2, :, :) = rho.*a.*u;
    Fcy(3, :, :) = rho.*a.*v;
    Fcy(4, :, :) = rho.*a.*H;
    
    
end