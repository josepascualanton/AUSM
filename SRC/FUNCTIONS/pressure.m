function p = pressure(U)
    gamma = 1.4;
    
    rho(:, :) = U(1,:,:);
    u2(:, :) = U(2, :, :);
    u3(:, :) = U(3, :, :);
    u4(:, :) = U(4, :, :);

    u = u2 ./ rho;
    v = u3 ./ rho;
    p   = (gamma-1) .* (u4 - 0.5 .* rho .* (u.^2 + v.^2));

end