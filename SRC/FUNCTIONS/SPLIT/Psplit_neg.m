function P = Psplit_neg(P, M)
    
    if abs(M) < 1
        P = P*((M - 1)^2)*(2 + M)/4;
    else    
        P = P*(M - abs(M))/M;
    end
end