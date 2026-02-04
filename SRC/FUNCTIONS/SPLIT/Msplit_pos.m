function M = Msplit_pos(M)
    
    if abs(M) < 1
        M = ((M+1)^2)/4;
    else    
        M = (M + abs(M))/2;
    end
end