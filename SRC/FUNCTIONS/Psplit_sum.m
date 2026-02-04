function P = Psplit_sum(Pr, Pl, Ml, Mr)
    
    P = Psplit_pos(Pl, Ml) + Psplit_neg(Pr, Mr);
end