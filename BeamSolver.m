classdef BeamSolver < handle
    
    properties (Access = private)
        beamParams
        connec
        externalForce
    end

    methods (Access = public)
        
        function obj = BeamSolver(cParams)
            obj.init(cParams);
        end

        function [K,F,u,r] = compute(obj)
            Kel     = obj.computeElementalStiffnessMatrix();
            fel     = obj.computeElementalForceVector();
            [K,F]   = obj.computeGlobalAssembly(Kel,fel);
            [up,vp] = obj.computeBoundaryConditions();
            F       = obj.computePointLoads(F);
            [u,r]   = obj.computeSystem(K,F,up,vp);
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.beamParams    = cParams.beamParams;
            obj.connec        = cParams.connec;
            obj.externalForce = cParams.externalForce;
        end

        function Kel = computeElementalStiffnessMatrix(obj)
            s.beamP  = obj.beamParams;
            s.con    = obj.connec;
            K        = ElemenetalStiffnessMatrixComputer(s);
            Kel      = K.computeStiffnessMatrix();
        end

        function fel = computeElementalForceVector(obj)
            s.beamP = obj.beamParams;
            s.con   = obj.connec;
            f       = ElementalForceVectorComputer(s);
            fel     = f.computeForceVector();
        end

        function [K,F] = computeGlobalAssembly(obj,Kel,fel)
            s.beamP = obj.beamParams;
            s.con   = obj.connec;
            s.Kel   = Kel;
            s.fel   = fel;
            globAssembly = Assembler(s);
            [K,F]   = globAssembly.computeAssembly();
        end

        function [up,vp] = computeBoundaryConditions(obj)
            s.beamP = obj.beamParams;
            s.con   = obj.connec;
            bc      = BoundaryConditionsApplier(s);
            [up,vp] = bc.computeBC();
        end

        function F = computePointLoads(obj,F)
            s.nDofsNode = obj.beamParams.numDOFsNode;
            s.extF      = obj.externalForce;
            s.F         = F;
            pL          = PointLoadsApplier(s);
            pL.computePointLoads();
            F = pL.F;
        end

        function [u,r] = computeSystem(obj,K,F,up,vp)
            s.nDofTotal = obj.beamParams.totalDOFs;
            s.K         = K;
            s.F         = F;
            s.up        = up;
            s.vp        = vp;
            solv  = SystemSolver(s);
            [u,r] = solv.computeSystem();
        end

    end

end