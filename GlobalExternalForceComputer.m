classdef GlobalExternalForceComputer < handle
    
    properties (Access = private)
        beamParams
        Connec
        fElem
    end

    properties (Access = public)
        F
    end

    methods (Access = public)
        
        function obj = GlobalExternalForceComputer(cParams,fel)
            obj.init(cParams,fel);
        end

        function computeExternalForces(obj)
            obj.compute();
        end

    end

    methods (Access = private)
        
        function init(obj,cParams,fel)
            obj.beamParams = cParams.beamP;
            obj.Connec     = cParams.con;
            obj.fElem      = fel;
        end

        function compute(obj)
            ndof = obj.beamParams.totalDOFs;
            nel  = obj.beamParams.numElements;
            nne  = obj.beamParams.numNodesElem;
            ni   = obj.beamParams.numDOFsNode;
            Td   = obj.Connec.dofsConnec;
            
            GlobF = zeros(ndof,1);
            for e = 1:nel
                for i = 1:(nne*ni)
                    GlobF(Td(e,i)) = GlobF(Td(e,i)) + obj.fElem(i,e);
                end
            end
            obj.F = GlobF;
        end

    end

end