function F = F_RIGH(U, MESH)

size = size(MESH);

for i = 2:size(1) - 1
    for j = 2:size(2) - 1
        M = Msplit_sum(Mr, Ml);
        F(i, j) = 1/2*M*(F(i, j) + F(i + 1, j)) - 1/2*abs(M)*(F(i + 1, j) - F(i, j)) + P*[0; 1; 0; 0];
    end
end
end