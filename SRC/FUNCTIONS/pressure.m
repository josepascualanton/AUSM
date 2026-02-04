function p = pressure(U)

    gamma = 1.4;

    rho = squeeze(U(1,:,:));
    u2  = squeeze(U(2,:,:));
    u3  = squeeze(U(3,:,:));
    u4  = squeeze(U(4,:,:));

    u = u2 ./ rho;
    v = u3 ./ rho;

    p = (gamma-1) .* (u4 - 0.5 .* rho .* (u.^2 + v.^2));
    p = max(p, 1);

end
