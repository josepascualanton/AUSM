function F = F_LEFT(Fc, M, P)

[Ny, Nx] = size(M);
F = zeros(4, Ny, Nx);

for i = 2:(Nx - 1)
    for j = 2:(Ny - 1)
        Mij = Msplit_sum(M(i - 1, j), M(i, j));
        Pij = Psplit_sum(P(i - 1, j), P(i, j), M(i, j), M(i + 1, j));
        
        F(:, i, j) = 1/2*Mij*(Fc(:, i - 1, j) + Fc(:, i, j)) - 1/2*abs(Mij)*(Fc(:, i, j) - Fc(:, i - 1, j)) + Pij*[0; -1; 0; 0];
    end
end
end