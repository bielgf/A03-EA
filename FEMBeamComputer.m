classdef FEMBeamComputer < handle
    
    properties (Access = private)
        beamWidth
        h1
        h2
        N1
        N2
        N3
        x_prim
        Tm
        Tn
        materialProp
        open
        x_0_prim
        x_s_prim
        I_xx_prim
        J
        xBendMoment
        yBendMoment
        zBendMoment
        xShearForce
        yShearForce
        mD2
        E
        G
        xnod
        numDOFsNode
        numNodesElem
        numElements
        ndof
        nodalConnec
        TmD2
        TdD2
        fe
        me
        pD2
        be
        Me
        g
        ze
        xi_p
        wingspan
        rhoInf
        vInf
        chord
        Cl
        lambda
        za
        zm
    end

    properties (Access = public)
        sigma
        s_norm
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

        function computeGeoDiscret(obj)
            gd.N1        = obj.N1;
            gd.N2        = obj.N2;
            gd.N3        = obj.N3;
            gd.h1        = obj.h1;
            gd.h2        = obj.h2;
            gd.beamWidth = obj.beamWidth;
            geoDiscret   = GeometricDiscretizationSolver(gd);
            [obj.x_prim,obj.Tm,obj.Tn] = geoDiscret.compute();
        end
        
        function computeSectionSolver(obj)
            s.x_prim   = obj.x_prim;
            s.Tn       = obj.Tn;
            s.materialProp      = obj.materialProp;
            s.Tm       = obj.Tm;
            s.open     = obj.open;
            s.xBendMoment = obj.xBendMoment;
            s.yBendMoment = obj.yBendMoment;
            s.zBendMoment = obj.zBendMoment;
            s.xShearForce = obj.xShearForce;
            s.yShearForce = obj.yShearForce;
            sec        = SectionSolver(s);
            [obj.sigma,obj.s_norm,obj.tau_s,obj.s_shear,obj.tau_t,obj.s_tor,obj.x_s_prim,obj.I_xx_prim,obj.J] = sec.compute();
        end

        function computemD2(obj)
            obj.mD2 = [obj.E   obj.G   obj.I_xx_prim  obj.J];
        end

        function computeForceMomentElem(obj)
            fm.numDOFsNode       = obj.numDOFsNode;
            fm.numElements      = obj.numElements;
            fm.wingspan        = obj.wingspan;
            fm.rhoInf   = obj.rhoInf;
            fm.vInf     = obj.vInf;
            fm.chord    = obj.chord;
            fm.Cl       = obj.Cl;
            fm.lambda   = obj.lambda;
            fm.g        = obj.g;
            fm.d        = obj.beamWidth;
            fm.xi_p     = obj.xi_p;
            fm.za       = obj.za;
            fm.zm       = obj.zm;
            fm.x_s_prim = obj.x_s_prim;
            forceMomentElem = ForceMomentElemCompute(fm);
            [obj.xnod,obj.fe,obj.me,obj.ndof] = forceMomentElem.compute();
        end

        function computeBeamSolver(obj)
            bm.numNodesElem      = obj.numNodesElem;
            bm.numDOFsNode       = obj.numDOFsNode;
            bm.numElements      = obj.numElements;
            bm.xnod     = obj.xnod;
            bm.nodalConnec     = obj.nodalConnec;
            bm.TdD2     = obj.TdD2;
            bm.mD2      = obj.mD2;
            bm.TmD2     = obj.TmD2;
            bm.pD2      = obj.pD2;
            bm.fe       = obj.fe;
            bm.me       = obj.me;
            bm.ndof     = obj.ndof;
            bm.be       = obj.be;
            bm.Me       = obj.Me;
            bm.g        = obj.g;
            bm.beamWidth        = obj.beamWidth;
            bm.xi_p     = obj.xi_p;
            bm.x_s_prim = obj.x_s_prim;
            bm.ze       = obj.ze;
            beam        = BeamSolver(bm);
            [obj.K,obj.F,obj.u,obj.r] = beam.compute();
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.beamWidth    = cParams.beamWidth;
            obj.h1           = cParams.h1;
            obj.h2           = cParams.h2;
            obj.N1           = cParams.N1;
            obj.N2           = cParams.N2;
            obj.N3           = cParams.N3;
            obj.materialProp = cParams.materialProp;
            obj.open         = cParams.open;
            obj.xBendMoment  = cParams.xBendMoment;
            obj.yBendMoment  = cParams.yBendMoment;
            obj.zBendMoment  = cParams.zBendMoment;
            obj.xShearForce  = cParams.xShearForce;
            obj.yShearForce  = cParams.yShearForce;
            obj.E            = cParams.E;
            obj.G            = cParams.G;
            obj.numDOFsNode  = cParams.numDOFsNode;
            obj.numNodesElem = cParams.numNodesElem;
            obj.numElements  = cParams.numElements;
            obj.wingspan     = cParams.wingspan;
            obj.rhoInf       = cParams.rhoInf;
            obj.vInf         = cParams.vInf;
            obj.chord        = cParams.chord;
            obj.Cl           = cParams.Cl;
            obj.lambda       = cParams.lambda;
            obj.nodalConnec  = cParams.nodalConnec;
            obj.TmD2     = cParams.materialConnec;
            obj.TdD2     = cParams.dofsConnec;
            obj.pD2      = cParams.fixedNodes;
            obj.be       = cParams.xEngine;
            obj.Me       = cParams.engineMass;
            obj.g        = cParams.g;
            obj.ze       = cParams.zEngine;
            obj.za       = cParams.aeroCenter;
            obj.zm       = cParams.centerOfMass;
            obj.xi_p     = cParams.chiP;
        end

    end

end