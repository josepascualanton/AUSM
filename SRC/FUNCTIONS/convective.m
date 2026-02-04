function Fc = convective(U)
    [~, Ny, Nx] = size(U);
    gamma = 1.4;
    Fc = zeros(4, Ny, Nx);

    rho(:, :) = U(1,:,:);
    u2(:, :) = U(2, :, :);
    u3(:, :) = U(3, :, :);
    u4(:, :) = U(4, :, :);

    u = u2 ./ rho;
    v = u3 ./ rho;
    
    p   = (gamma-1) .* (u4 - 0.5 .* rho .* (u.^2 + v.^2));
    a   = sqrt(gamma .* p ./ rho);
    H   = (u4 + p) ./ rho;
    M   = u ./ a;

    Fc(1, :, :) = rho.*a;
    Fc(2, :, :) = rho.*a.*u;
    Fc(3, :, :) = rho.*a.*v;
    Fc(4, :, :) = rho.*a.*H;
    
end