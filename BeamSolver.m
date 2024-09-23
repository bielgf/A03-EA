classdef BeamSolver < handle
    
    properties (Access = private)
        beamParams
        connec
        beamProp
        externalForce
    end

    methods (Access = public)
        
        function obj = BeamSolver(cParams)
            obj.init(cParams);
        end

        function [K,F,u,r] = compute(obj)
            s.beamP  = obj.beamParams;
            s.con    = obj.connec;
            s.beamPr = obj.beamProp;
            s.extF   = obj.externalForce;

            Kel = obj.computeElementalStiffnessMatrix();
            s.Kel = Kel;

            fF = ElementForceMatrixComputer(s);
            fF.computeForceVector();
            s.fel = fF.fel;

            aF = Assembler(s);
            aF.computeAssembly();
            s.K = aF.K;
            s.F = aF.F;

            bc = BoundaryConditionsApplier(s);
            bc.computeBC();
            s.up = bc.up;
            s.vp = bc.vp;

            pL = PointLoadsApplier(s);
            pL.computePointLoads();
            s.F = pL.F;

            sS = SystemSolver(s);
            sS.computeSystem();
            u = sS.u;
            r = sS.r;

            K = s.K;
            F = s.F;
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.beamParams    = cParams.beamParams;
            obj.connec        = cParams.connec;
            obj.beamProp      = cParams.beamProp;
            obj.externalForce = cParams.externalForce;
        end

        function Kel = computeElementalStiffnessMatrix(obj)
            s.beamP  = obj.beamParams;
            s.con    = obj.connec;
            s.beamPr = obj.beamProp;
            K        = ElemenetalStiffnessMatrixComputer(s);
            Kel      = K.computeStiffnessMatrix();
        end

    end

end