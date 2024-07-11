classdef GlobalStiffnessMatrixComputer

    properties
        K_Global
    end

    methods (Access = public)

        function obj = assemblyFunction(obj,data,Td,Kel)
            
            obj.K_Global = zeros(data.ndof,data.ndof);
            for e = 1:data.nel
                for i = 1:(data.nne*data.ni)
                    for j = 1:(data.nne*data.ni)
                        obj.K_Global(Td(e,i),Td(e,j)) = obj.K_Global(Td(e,i),Td(e,j)) + Kel(i,j,e);
                    end
                end
            end

        end

    end

end