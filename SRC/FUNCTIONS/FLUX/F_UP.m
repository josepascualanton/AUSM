function F = F_UP(Fc, Mv, P)

[Ny, Nx] = size(Mv);
F = zeros(4, Ny, Nx);

for j = 2:Ny-1
    for i = 2:Nx-1

        Mij = Msplit_sum(Mv(j, i), Mv(j+1, i));
        Pij = Psplit_sum(P(j, i), P(j+1, i), Mv(j, i), Mv(j+1, i));

        F(:, j, i) = ...
            0.5*Mij*(Fc(:, j, i) + Fc(:, j+1, i)) ...
          - 0.5*abs(Mij)*(Fc(:, j+1, i) - Fc(:, j, i)) ...
          + Pij*[0; 0; 1; 0];
    end
end
end
