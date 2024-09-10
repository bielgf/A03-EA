classdef FEMBeamComputer < handle
    
    % ToDo:
    %- Simplificar properties
    %- Eliminar properties publiques que no fem servir
    %- Veure comentaris específics en el codi

    properties (Access = private)
        geomParams
        connec

        h1
        h2
        N1
        N2
        N3
        
        xSection
        materialSecConnec
        nodalSecConnec
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
        wingspan
        rhoInf
        vInf
        chord
        Cl
        lambda
        aeroCenter
        centerOfMass
        chiP
    end

    properties (Access = public)
        sigma
        normalStress
        tau_s
        s_shear
        tau_t
        s_tor
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
            obj.h1             = cParams.h1; % delete
            obj.h2             = cParams.h2; % delete
            obj.N1             = cParams.N1; % delete
            obj.N2             = cParams.N2; % delete
            obj.N3             = cParams.N3; % delete
            obj.materialProp   = cParams.materialProp;
            obj.open           = cParams.open;
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
            obj.wingspan       = cParams.wingspan;
            obj.rhoInf         = cParams.rhoInf;
            obj.vInf           = cParams.vInf;
            obj.chord          = cParams.chord;
            obj.Cl             = cParams.Cl;
            obj.lambda         = cParams.lambda;
            obj.nodalConnec    = cParams.nodalConnec;
            obj.materialConnec = cParams.materialConnec;
            obj.dofsConnec     = cParams.dofsConnec;
            obj.fixedNodes     = cParams.fixedNodes;
            obj.xEngine        = cParams.xEngine;
            obj.engineMass     = cParams.engineMass;
            obj.g              = cParams.g;
            obj.zEngine        = cParams.zEngine;
            obj.aeroCenter     = cParams.aeroCenter;
            obj.centerOfMass   = cParams.centerOfMass;
            obj.chiP           = cParams.chiP;
        end

        function computeGeoDiscret(obj)
            gd.N1        = obj.N1;
            gd.N2        = obj.N2;
            gd.N3        = obj.N3;
            gd.h1        = obj.h1;
            gd.h2        = obj.h2;
            gd.beamWidth = obj.geomParams.beamWidth;
            geoDiscret   = GeometricDiscretizationSolver(gd);
            [obj.xSection,obj.materialSecConnec,obj.nodalSecConnec] = geoDiscret.compute();
        end

        function computeSectionSolver(obj)
            s.xSection          = obj.xSection;
            s.nodalSecConnec    = obj.nodalSecConnec;
            s.materialProp      = obj.materialProp;
            s.materialSecConnec = obj.materialSecConnec;
            s.open              = obj.open;
            s.xBendMoment       = obj.xBendMoment;
            s.yBendMoment       = obj.yBendMoment;
            s.zBendMoment       = obj.zBendMoment;
            s.xShearForce       = obj.xShearForce;
            s.yShearForce       = obj.yShearForce;
            sec                 = SectionSolver(s);
            [obj.sigma,obj.normalStress,obj.tau_s,obj.s_shear,obj.tau_t,obj.s_tor,obj.xShearCenter,obj.xSecInertia,obj.J] = sec.compute();
        end

        function computeBeamProp(obj)
            obj.beamProp = [obj.E   obj.G   obj.xSecInertia  obj.J];
        end

        function computeForceMomentElem(obj)
            fm.numDOFsNode  = obj.numDOFsNode;
            fm.numElements  = obj.numElements;
            fm.wingspan     = obj.wingspan;
            fm.rhoInf       = obj.rhoInf;
            fm.vInf         = obj.vInf;
            fm.chord        = obj.chord;
            fm.Cl           = obj.Cl;
            fm.lambda       = obj.lambda;
            fm.g            = obj.g;
            fm.beamWidth    = obj.geomParams.beamWidth;
            fm.chiP         = obj.chiP;
            fm.aeroCenter   = obj.aeroCenter;
            fm.centerOfMass = obj.centerOfMass;
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
            obj.externalForce = [F11, 1, -obj.engineMass*obj.g;
                F21, 3, -obj.engineMass*obj.g*(((obj.geomParams.beamWidth + obj.chiP) - obj.xShearCenter) - obj.zEngine)];

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