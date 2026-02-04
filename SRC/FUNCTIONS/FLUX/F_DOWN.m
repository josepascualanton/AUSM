function F = F_DOWN(Fc, Mv, P)

[Ny, Nx] = size(Mv);
F = zeros(4, Ny, Nx);

for j = 2:Ny-1
    for i = 2:Nx-1

        Mij = Msplit_sum(Mv(j-1, i), Mv(j, i));
        Pij = Psplit_sum(P(j-1, i), P(j, i), Mv(j-1, i), Mv(j, i));

        F(:, j, i) = ...
            0.5*Mij*(Fc(:, j-1, i) + Fc(:, j, i)) ...
          - 0.5*abs(Mij)*(Fc(:, j, i) - Fc(:, j-1, i)) ...
          + Pij*[0; 0; -1; 0];
    end
end
end
