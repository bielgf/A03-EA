classdef GlobalStiffnessMatrixComputer < handle

    properties (Access = private)

    end

    properties (Access = public)
        KGlobal
    end

    methods (Access = public)

        function obj = GlobalStiffnessMatrixComputer() % passar els 4 arguments
            obj.init(); % Inicialització

            % Calculs previs (en aquest cas no cal)
        end

        function compute(obj)
            obj.assemblyMatrix();
        end

    end

    methods (Access = private)
        function init(obj,...)
                obj.Kelem = ...;
                %...
        end

            function assemblyMatrix(obj)
                obj.KGlobal = zeros(data.ndof,data.ndof);
                for e = 1:data.nel
                    for i = 1:(data.nne*data.ni)
                        for j = 1:(data.nne*data.ni)
                            obj.KGlobal(Td(e,i),Td(e,j)) = obj.KGlobal(Td(e,i),Td(e,j)) + Kel(i,j,e);
                        end
                    end
                end
            end
    end

end