classdef GlobalStiffnessMatrixComputer < handle

    properties (Access = private)
        Kelem
        Tnod
        nel
        nne
        ni
        ndof
    end

    properties (Access = public)
        KGlobal
    end

    methods (Access = public)

        function obj = GlobalStiffnessMatrixComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.assemblyMatrix();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.Kelem = cParams.Kel;
            obj.Tnod  = cParams.Td;
            obj.nel   = cParams.nel;
            obj.nne   = cParams.nne;
            obj.ni    = cParams.ni;
            obj.ndof  = cParams.ndof;
        end

        function assemblyMatrix(obj)
            obj.KGlobal = zeros(obj.ndof,obj.ndof);
            for e = 1:obj.nel
                for i = 1:(obj.nne*obj.ni)
                    for j = 1:(obj.nne*obj.ni)
                        obj.KGlobal(obj.Tnod(e,i),obj.Tnod(e,j)) = obj.KGlobal(obj.Tnod(e,i),obj.Tnod(e,j)) + obj.Kelem(i,j,e);
                    end
                end
            end
        end

    end

end