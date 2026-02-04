function F = F_UP(Fc, M, P)

[Ny, Nx] = size(M);
F = zeros(4, Ny, Nx);

for i = 2: (Nx - 1)
    for j = 2: (Ny - 1)

        Mij = Msplit_sum(M(i, j), M(i, j + 1));
        Pij = Psplit_sum(P(i, j), P(i, j + 1), M(i, j), M(i, j +1));
        
        F(:, i, j) = 1/2*Mij*(Fc(:, i, j) + Fc(:, i, j + 1)) - 1/2*abs(Mij)*(Fc(:, i, j + 1) - Fc(:, i, j)) + Pij*[0; 0; 1; 0];
    end
end
end