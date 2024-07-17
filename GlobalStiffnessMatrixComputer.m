classdef GlobalStiffnessMatrixComputer < handle

    properties (Access = private)
        Kelem
        Tnod
        data
    end

    properties (Access = public)
        KGlobal
    end

    methods (Access = public)

        function obj = GlobalStiffnessMatrixComputer(Kelem,Tnod,data)
            obj.init(Kelem,Tnod,data);
        end

        function compute(obj)
            obj.assemblyMatrix();
        end

    end

    methods (Access = private)

        function init(obj,Kelem,Tnod,data)
            obj.Kelem = Kelem;
            obj.Tnod = Tnod;
            obj.data = data;
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