classdef FEMBeamComputer < handle
    
    % ToDo:
    %- Simplificar properties
    %- Veure comentaris específics en el codi

    properties (Access = private)
        geomParams

        connec

        aeroParams
        
        materialProp
        open
        xShearCenter
        xSecInertia
        J
        xBendMoment
        yBendMoment
        zBendMoment
        xShearForce
        yShearForce
        beamProp
        E
        G
        xGlobal
        numDOFsNode
        numNodesElem
        numElements
        totalDOFs
        nodalConnec
        materialConnec
        dofsConnec
        forceElem
        momentElem
        externalForce
        fixedNodes
        xEngine
        engineMass
        g
        zEngine

        lambda
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
%             obj.connec         = cParams.connec;
            obj.aeroParams     = cParams.aeroParams;

            obj.materialProp   = cParams.materialProp;
            obj.xBendMoment    = cParams.xBendMoment;
            obj.yBendMoment    = cParams.yBendMoment;
            obj.zBendMoment    = cParams.zBendMoment;
            obj.xShearForce    = cParams.xShearForce;
            obj.yShearForce    = cParams.yShearForce;
            obj.E              = cParams.E;
            obj.G              = cParams.G;
            obj.numDOFsNode    = cParams.numDOFsNode;
            obj.numNodesElem   = cParams.numNodesElem;
            obj.numElements    = cParams.numElements;
            obj.lambda         = cParams.lambda;
            obj.nodalConnec    = cParams.nodalConnec;
            obj.materialConnec = cParams.materialConnec;
            obj.dofsConnec     = cParams.dofsConnec;
            obj.fixedNodes     = cParams.fixedNodes;
            obj.xEngine        = cParams.xEngine;
            obj.engineMass     = cParams.engineMass;
            obj.zEngine        = cParams.zEngine;
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
            s.materialProp      = obj.materialProp;
            s.materialSecConnec = obj.geomParams.materialSecConnec;
            s.open              = obj.open;
            s.xBendMoment       = obj.xBendMoment;
            s.yBendMoment       = obj.yBendMoment;
            s.zBendMoment       = obj.zBendMoment;
            s.xShearForce       = obj.xShearForce;
            s.yShearForce       = obj.yShearForce;
            sec                 = SectionSolver(s);
            [~,~,~,~,~,~,obj.xShearCenter,obj.xSecInertia,obj.J] = sec.compute();
        end

        function computeBeamProp(obj)
            obj.beamProp = [obj.E   obj.G   obj.xSecInertia  obj.J];
        end

        function computeForceMomentElem(obj)
            fm.numDOFsNode  = obj.numDOFsNode;
            fm.numElements  = obj.numElements;
            fm.wingspan     = obj.geomParams.wingspan;
            fm.rhoInf       = obj.aeroParams.rhoInf;
            fm.vInf         = obj.aeroParams.vInf;
            fm.chord        = obj.geomParams.chord;
            fm.Cl           = obj.aeroParams.Cl;
            fm.lambda       = obj.lambda;
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
            xE = obj.xEngine;
            % ...
            F11 = find(xG == xE);
            F21 = F11;
            % ...
            obj.externalForce = [F11, 1, -obj.engineMass*obj.aeroParams.g;
                                 F21, 3, -obj.engineMass*obj.aeroParams.g*(((obj.geomParams.beamWidth + obj.geomParams.chiP) - obj.xShearCenter) - obj.zEngine)];

            % Create a class for external force
        end

        function computeBeamSolver(obj)
            bm.numNodesElem   = obj.numNodesElem;
            bm.numDOFsNode    = obj.numDOFsNode;
            bm.numElements    = obj.numElements;
            bm.xGlobal        = obj.xGlobal;
            bm.nodalConnec    = obj.nodalConnec;
            bm.dofsConnec     = obj.dofsConnec;
            bm.beamProp       = obj.beamProp;
            bm.materialConnec = obj.materialConnec;
            bm.fixedNodes     = obj.fixedNodes;
            bm.forceElem      = obj.forceElem;
            bm.momentElem     = obj.momentElem;
            bm.externalForce  = obj.externalForce;
            bm.totalDOFs      = obj.totalDOFs;
            beam              = BeamSolver(bm);
            [obj.K,obj.F,obj.u,obj.r] = beam.compute();
        end



    end

end