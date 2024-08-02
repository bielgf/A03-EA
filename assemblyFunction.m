function [K_global,F_global] = assemblyFunction(ndof,nel,nne,ni,Td,Kel,fel)
    F_global = zeros(ndof,1);

    for e = 1:nel
        for i = 1:(nne*ni)
            F_global(Td(e,i)) = F_global(Td(e,i)) + fel(i,e);
        end
    end 
    
    s.Kel  = Kel;
    s.Td   = Td;
    s.nel  = nel;
    s.nne  = nne;
    s.ni   = ni;
    s.ndof = ndof;
    gsmc = GlobalStiffnessMatrixComputer(s);
    gsmc.compute();
    K_global = gsmc.KGlobal;
end