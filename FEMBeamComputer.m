classdef FEMBeamComputer < handle
    
    % ToDo:
    %- Simplificar properties
    %- Veure comentaris específics en el codi

    properties (Access = private)
        geomParams
        beamParams
        aeroParams
        connec
        secStress
    end

    properties (Access = private)
        xShearCenter
        xSecInertia
        J
        beamProp
        xGlobal
        totalDOFs
        forceElem
        momentElem
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
            obj.secStress      = cParams.secStress;
        end

        function computeGeoDiscret(obj)
            gd = obj.geomParams;
            geoDiscret   = GeometricDiscretizationSolver(gd);
            [xSection,materialSecConnec,nodalSecConnec] = geoDiscret.compute();
            obj.geomParams.xSection = xSection;
            obj.geomParams.materialSecConnec = materialSecConnec;
            obj.geomParams.nodalSecConnec = nodalSecConnec;
        end

        function computeSectionSolver(obj)
            s.xSection          = obj.geomParams.xSection;
            s.nodalSecConnec    = obj.geomParams.nodalSecConnec;
            s.materialProp      = obj.geomParams.materialProp;
            s.materialSecConnec = obj.geomParams.materialSecConnec;
            s.open              = obj.geomParams.open;
            s.xBendMoment       = obj.secStress.xBendMoment;
            s.yBendMoment       = obj.secStress.yBendMoment;
            s.zBendMoment       = obj.secStress.zBendMoment;
            s.xShearForce       = obj.secStress.xShearForce;
            s.yShearForce       = obj.secStress.yShearForce;
            sec                 = SectionSolver(s);
            [~,~,~,~,~,~,obj.xShearCenter,obj.xSecInertia,obj.J] = sec.compute();
        end

        function computeBeamProp(obj)
            obj.beamProp = [obj.geomParams.E   obj.geomParams.G   obj.xSecInertia  obj.J];
        end

        function computeForceMomentElem(obj)
            fm.numDOFsNode  = obj.beamParams.numDOFsNode;
            fm.numElements  = obj.beamParams.numElements;
            fm.wingspan     = obj.geomParams.wingspan;
            fm.rhoInf       = obj.aeroParams.rhoInf;
            fm.vInf         = obj.aeroParams.vInf;
            fm.chord        = obj.geomParams.chord;
            fm.Cl           = obj.aeroParams.Cl;
            fm.lambda       = obj.beamParams.lambda;
            fm.g            = obj.aeroParams.g;
            fm.beamWidth    = obj.geomParams.beamWidth;
            fm.chiP         = obj.geomParams.chiP;
            fm.aeroCenter   = obj.geomParams.aeroCenter;
            fm.gravCenter   = obj.geomParams.gravCenter;
            fm.xShearCenter = obj.xShearCenter;
            forceMomentElem = ForceMomentElemCompute(fm);
            [obj.xGlobal,obj.forceElem,obj.momentElem,obj.totalDOFs] = forceMomentElem.compute();
        end

        function computeExternalForce(obj)
            xG = obj.xGlobal;
            xE = obj.beamParams.xEngine;
            eM = obj.beamParams.engineMass;
            % ...
            F11 = find(xG == xE);
            F21 = F11;
            % ...
            obj.externalForce = [F11, 1, -eM*obj.aeroParams.g;
                                 F21, 3, -eM*obj.aeroParams.g*(((obj.geomParams.beamWidth + obj.geomParams.chiP) - obj.xShearCenter) - obj.beamParams.zEngine)];

            % Create a class for external force
        end

        function computeBeamSolver(obj)
            bm.numNodesElem   = obj.beamParams.numNodesElem;
            bm.numDOFsNode    = obj.beamParams.numDOFsNode;
            bm.numElements    = obj.beamParams.numElements;
            bm.xGlobal        = obj.xGlobal;
            bm.nodalConnec    = obj.connec.nodalConnec;
            bm.dofsConnec     = obj.connec.dofsConnec;
            bm.materialConnec = obj.connec.materialConnec;
            bm.fixedNodes     = obj.connec.fixedNodes;
            bm.beamProp       = obj.beamProp;
            bm.forceElem      = obj.forceElem;
            bm.momentElem     = obj.momentElem;
            bm.externalForce  = obj.externalForce;
            bm.totalDOFs      = obj.totalDOFs;
            beam              = BeamSolver(bm);
            [obj.K,obj.F,obj.u,obj.r] = beam.compute();
        end

    end

end