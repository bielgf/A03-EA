classdef FEMBeamComputer < handle

    properties (Access = private)
        geomParams
        beamParams
        aeroParams
        connec
        beamProp
        externalForce
    end

    properties (Access = public)
        K
        F
        u
        r
    end

    methods (Access = public)
        
        function obj = FEMBeamComputer(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.computeGeoDiscret();
            obj.computeSectionSolver();
            obj.computeBeamProp();
            obj.computeForceMomentElem();
            obj.computeExternalForce();
            obj.computeBeamSolver();
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.geomParams     = cParams.geomParams;
            obj.beamParams     = cParams.beamParams;
            obj.aeroParams     = cParams.aeroParams;
            obj.connec         = cParams.connec;
        end

        function computeGeoDiscret(obj)
            gd = obj.geomParams;
            geoDiscret   = GeometricDiscretizationSolver(gd);
            [xSection,materialSecConnec,nodalSecConnec] = geoDiscret.compute();
            obj.geomParams.xSection          = xSection;
            obj.geomParams.materialSecConnec = materialSecConnec;
            obj.geomParams.nodalSecConnec    = nodalSecConnec;
        end

        function computeSectionSolver(obj)
            s   = obj.geomParams;
            sec = SectionSolver(s);
            [xShearCenter,xSecInertia,J] = sec.compute();
            obj.geomParams.xShearCenter = xShearCenter;
            obj.geomParams.xSecInertia  = xSecInertia;
            obj.geomParams.J            = J;
        end

        function computeBeamProp(obj)
            E   = obj.geomParams.E;
            G   = obj.geomParams.G;
            Ixx = obj.geomParams.xSecInertia;
            J   = obj.geomParams.J;
            obj.beamProp = [E   G   Ixx    J];
        end

        % -----------------------------------------------

        function computeForceMomentElem(obj)
            fm.beamWidth    = obj.geomParams.beamWidth;
            fm.chiP         = obj.geomParams.chiP;
            fm.aeroCenter   = obj.geomParams.aeroCenter;
            fm.gravCenter   = obj.geomParams.gravCenter;
            fm.xShearCenter = obj.geomParams.xShearCenter;
            fm.wingspan     = obj.geomParams.wingspan;
            fm.chord        = obj.geomParams.chord;
            fm.rhoInf       = obj.aeroParams.rhoInf;
            fm.vInf         = obj.aeroParams.vInf;            
            fm.Cl           = obj.aeroParams.Cl;
            fm.g            = obj.aeroParams.g;
            fm.lambda       = obj.beamParams.lambda;
            fm.numDOFsNode  = obj.beamParams.numDOFsNode;
            fm.numElements  = obj.beamParams.numElements;
            forceMomentElem = ForceMomentElemCompute(fm);
            [xGlobal,forceElem,momentElem,totalDOFs] = forceMomentElem.compute();
            obj.beamParams.xGlobal    = xGlobal;
            obj.beamParams.forceElem  = forceElem;
            obj.beamParams.momentElem = momentElem;
            obj.beamParams.totalDOFs  = totalDOFs;
        end

        function computeExternalForce(obj)
            ef.geomParams = obj.geomParams;
            ef.aeroParams = obj.aeroParams;
            ef.beamParams = obj.beamParams;
            efAssembly = ExternalForceAssembly(ef);
            efAssembly.assembly();
            obj.externalForce = efAssembly.externalForce;
        end

        function computeBeamSolver(obj)
            bm.numNodesElem   = obj.beamParams.numNodesElem;
            bm.numDOFsNode    = obj.beamParams.numDOFsNode;
            bm.numElements    = obj.beamParams.numElements;
            bm.xGlobal        = obj.beamParams.xGlobal;
            bm.forceElem      = obj.beamParams.forceElem;
            bm.momentElem     = obj.beamParams.momentElem;
            bm.totalDOFs      = obj.beamParams.totalDOFs;
            bm.nodalConnec    = obj.connec.nodalConnec;
            bm.dofsConnec     = obj.connec.dofsConnec;
            bm.materialConnec = obj.connec.materialConnec;
            bm.fixedNodes     = obj.connec.fixedNodes;
            bm.externalForce  = obj.externalForce;
            bm.beamProp       = obj.beamProp;
            beam              = BeamSolver(bm);
            [obj.K,obj.F,obj.u,obj.r] = beam.compute();
        end

    end

end