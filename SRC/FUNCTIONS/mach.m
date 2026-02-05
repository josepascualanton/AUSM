function [Mu, Mv] = mach(U)

    gamma = 1.4;

    rho = squeeze(U(1,:,:));
    u2 = squeeze(U(2,:,:));
    u3 = squeeze(U(3,:,:));
    u4 = squeeze(U(4,:,:));

    ux = u2 ./ rho;
    uy = u3 ./ rho;

    p = (gamma-1) .* (u4 - 0.5 .* rho .* (ux.^2 + uy.^2));
    p = max(p, eps);

    a = sqrt(gamma .* p ./ rho);
    
    Mu = ux./a;
    Mv = uy./a;
end