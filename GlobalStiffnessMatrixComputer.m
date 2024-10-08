classdef GlobalStiffnessMatrixComputer < handle

    properties (Access = private)
        beamParams
        Connec
        KElem
    end

    properties (Access = public)
        K
    end

    methods (Access = public)

        function obj = GlobalStiffnessMatrixComputer(cParams,Kel)
            obj.init(cParams,Kel);
        end

        function computeGlobalStiffnessMatrix(obj)
            obj.compute();
        end

    end

    methods (Access = private)

        function init(obj,cParams,Kel)
            obj.beamParams = cParams.beamP;
            obj.Connec     = cParams.con;
            obj.KElem      = Kel;
        end

        function compute(obj)
            ndof = obj.beamParams.totalDOFs;
            nel  = obj.beamParams.numElements;
            nne  = obj.beamParams.numNodesElem;
            ni   = obj.beamParams.numDOFsNode;
            Td = obj.Connec.dofsConnec;
            
            KGlobal = zeros(ndof,ndof);
            for e = 1:nel
                for i = 1:(nne*ni)
                    for j = 1:(nne*ni)
                        KGlobal(Td(e,i),Td(e,j)) = KGlobal(Td(e,i),Td(e,j)) + obj.KElem(i,j,e);
                    end
                end
            end
            obj.K = KGlobal;
        end

    end

end