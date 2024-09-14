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

            sF = StiffnessFunctionClass(s);
            sF.computeStiffnessMatrix();
            s.Kel = sF.Kel;

            fF = ForceFunctionClass(s);
            fF.computeForceVector();
            s.fel = fF.fel;

            aF = AssemblyFunctionClass(s);
            aF.computeAssembly();
            s.K = aF.K;
            s.F = aF.F;

            bc = applyBoundaryConditionsClass(s);
            bc.computeBC();
            s.up = bc.up;
            s.vp = bc.vp;

            pL = pointLoadsClass(s);
            pL.computePointLoads();
            s.F = pL.F;

            sS = SolveSystemClass(s);
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

    end

end