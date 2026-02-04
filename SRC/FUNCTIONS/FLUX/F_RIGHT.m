function F = F_RIGHT(Fc, M, P)

[Ny, Nx] = size(M);
F = zeros(4, Ny, Nx);
for i = 2:(Nx - 1)
    for j = 2:(Ny - 1)
        Mij = Msplit_sum(M(i, j), M(i + 1, j));
        Pij = Psplit_sum(P(i, j), P(i + 1, j), M(i, j), M(i + 1, j));
        F(:, i, j) = 1/2*Mij*(Fc(:, i, j) + Fc(:, i + 1, j)) - 1/2*abs(Mij)*(Fc(:, i + 1, j) - Fc(:, i, j)) + Pij*[0; 1; 0; 0];
    end
end
end