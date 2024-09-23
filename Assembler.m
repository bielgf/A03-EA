classdef Assembler < handle

    properties (Access = private)
        beamParams
        connec
        KElem
        fElem
    end

    methods (Access = public)
        
        function obj = Assembler(cParams)
            obj.init(cParams);
        end

        function [K,F] = computeAssembly(obj)
            [K,F] = obj.compute;
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.beamParams = cParams.beamP;
            obj.connec     = cParams.con;
            obj.KElem      = cParams.Kel;
            obj.fElem      = cParams.fel;
        end

        function [K,F] = compute(obj)
            s.beamP = obj.beamParams;
            s.con   = obj.connec;
            Kel     = obj.KElem;
            fel     = obj.fElem;

            gsmc = GlobalStiffnessMatrixComputer(s,Kel);
            gsmc.computeGlobalStiffnessMatrix();
            K = gsmc.K;

            gefc = GlobalExternalForceComputer(s,fel);
            gefc.computeExternalForces();
            F = gefc.F;
        end

    end
    
end