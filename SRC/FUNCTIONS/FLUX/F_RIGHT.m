function F = F_RIGHT(Fc, Mu, P)

[Ny, Nx] = size(Mu);
F = zeros(4, Ny, Nx);

for j = 2:Ny-1
    for i = 2:Nx-1

        Mij = Msplit_sum(Mu(j, i), Mu(j, i+1));
        Pij = Psplit_sum(P(j, i), P(j, i+1), Mu(j, i), Mu(j, i+1));

        F(:, j, i) = ...
            0.5*Mij*(Fc(:, j, i) + Fc(:, j, i+1)) ...
          - 0.5*abs(Mij)*(Fc(:, j, i+1) - Fc(:, j, i)) ...
          + Pij*[0; 1; 0; 0];
    end
end
end
