classdef AssemblyFunctionClass < handle

    properties (Access = private)
        beamParams
        connec
        KElem
        fElem
    end

    properties (Access = public)
        K
        F
    end

    methods (Access = public)
        
        function obj = AssemblyFunctionClass(cParams)
            obj.init(cParams);
        end

        function computeAssembly(obj)
            obj.compute;
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.beamParams = cParams.beamP;
            obj.connec     = cParams.con;
            obj.KElem      = cParams.Kel;
            obj.fElem      = cParams.fel;
        end

        function compute(obj)
            ndof = obj.beamParams.totalDOFs;
            nel  = obj.beamParams.numElements;
            nne  = obj.beamParams.numNodesElem;
            ni   = obj.beamParams.numDOFsNode;
            Td   = obj.connec.dofsConnec;
            Kel  = obj.KElem;
            fel  = obj.fElem;

            obj.F = zeros(ndof,1);

            for e = 1:nel
                for i = 1:(nne*ni)
                    obj.F(Td(e,i)) = obj.F(Td(e,i)) + fel(i,e);
                end
            end % private function + cridar class global ext forces
            
            s.Kel  = Kel;
            s.Td   = Td;
            s.nel  = nel;
            s.nne  = nne;
            s.ni   = ni;
            s.ndof = ndof;
            gsmc = GlobalStiffnessMatrixComputer(s);
            gsmc.compute();
            obj.K = gsmc.KGlobal; % private function
        end

    end
    
end