classdef FEMBeamComputer < handle

    properties (Access = private)
        geomParams
        beamParams
        aeroParams
        connec
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
            obj.computeForceMomentElem();
            obj.computeBeamSolver();
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.geomParams = cParams.geomParams;
            obj.beamParams = cParams.beamParams;
            obj.aeroParams = cParams.aeroParams;
            obj.connec     = cParams.connec;
        end

        function computeGeoDiscret(obj)
            s = obj.geomParams;
            geoDiscret   = GeometricDiscretizationSolver(s);
            [xSection,materialSecConnec,nodalSecConnec] = geoDiscret.compute();
            obj.geomParams.xSection          = xSection;
            obj.geomParams.materialSecConnec = materialSecConnec;
            obj.geomParams.nodalSecConnec    = nodalSecConnec;
        end

        function computeSectionSolver(obj)
            s   = obj.geomParams;
            sec = SectionSolver(s);
            [xShearCenter,xSecInertia,J] = sec.compute();
            obj.geomParams.xShearCenter  = xShearCenter;
            obj.geomParams.xSecInertia   = xSecInertia;
            obj.geomParams.J             = J;
        end

        function computeForceMomentElem(obj)
            s.geoParams   = obj.geomParams;
            s.aeroParams  = obj.aeroParams;
            s.lambda      = obj.beamParams.lambda;
            s.numDOFsNode = obj.beamParams.numDOFsNode;
            s.numElements = obj.beamParams.numElements;
            forceMomentElem = ForceMomentElemCompute(s);
            [xGlobal,forceElem,momentElem,totalDOFs] = forceMomentElem.compute();
            obj.beamParams.xGlobal    = xGlobal;
            obj.beamParams.forceElem  = forceElem;
            obj.beamParams.momentElem = momentElem;
            obj.beamParams.totalDOFs  = totalDOFs;
        end

        function extF = computeExternalForce(obj)
            s.geomParams = obj.geomParams;
            s.aeroParams = obj.aeroParams;
            s.beamParams = obj.beamParams;
            efAssembly = ExternalForceAssembly(s);
            efAssembly.assembly();
            extF = efAssembly.externalForce;
        end

        function computeBeamSolver(obj)
            s.beamParams      = obj.beamParams;
            s.connec          = obj.connec;
            s.externalForce   = obj.computeExternalForce();
            s.beamParams.Prop = obj.computeBeamProperties;
            beam            = BeamSolver(s);
            [obj.K,obj.F,obj.u,obj.r] = beam.compute();
        end

        function beamProp = computeBeamProperties(obj)
            E   = obj.geomParams.E;
            G   = obj.geomParams.G;
            Ixx = obj.geomParams.xSecInertia;
            J   = obj.geomParams.J;
            beamProp = [E   G   Ixx    J];
        end

    end

end