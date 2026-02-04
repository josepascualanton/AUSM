function M = mach(U)

[~, Ny, Nx] = size(U);
    
    M = zeros(Ny, Nx);
    gamma = 1.4;

    rho(:, :) = U(1, :, :);
    
    u2(:, :) = U(2, :, :);
    u3(:, :) = U(3, :, :);
    
    ux = u2 ./ rho;
    uy = u3 ./ rho;

    c = 300;

    M = sqrt(ux.^2 + uy.^2)/c;
end