classdef GlobalStiffnessMatrixComputer < handle

    properties (Access = private)
        Kelem
        Tnod
        data % treure
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
            obj.Tnod = cParams.Td;
            obj.data = cParams.data; % treure
            % afegir:
            % nel
            % nne
            % ni
            % ndof
        end

        function assemblyMatrix(obj)
            obj.KGlobal = zeros(obj.data.ndof,obj.data.ndof);
            for e = 1:obj.data.nel
                for i = 1:(obj.data.nne*obj.data.ni)
                    for j = 1:(obj.data.nne*obj.data.ni)
                        obj.KGlobal(obj.Tnod(e,i),obj.Tnod(e,j)) = obj.KGlobal(obj.Tnod(e,i),obj.Tnod(e,j)) + obj.Kelem(i,j,e);
                    end
                end
            end
        end

    end

end