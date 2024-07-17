function [K_global,F_global] = assemblyFunction(data,Td,Kel,fel)

    F_global = zeros(data.ndof,1);

    for e = 1:data.nel
        for i = 1:(data.nne*data.ni)
            F_global(Td(e,i)) = F_global(Td(e,i)) + fel(i,e);
        end
    end 

    gsmc = GlobalStiffnessMatrixComputer(Kel,Td,data);
    gsmc.compute();
    K_global = gsmc.KGlobal;
end