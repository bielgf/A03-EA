classdef FEMBeamComputer < handle
    
    properties (Access = private)
        d
        h1
        h2
        N1
        N2
        N3
        x_prim
        Tm
        Tn
        mD1
        open
        x_0_prim
        y_0_prim
        x_s_prim
        y_s_prim
        A_tot
        I_xx_prim
        I_yy_prim
        I_xy_prim
        J
        A_in
        M_x_prim
        M_y_prim
        M_z_prim
        S_x_prim
        S_y_prim
        mD2
        E
        G
        xnod
        ni
        nne
        nel
        ndof
        TnD2
        TmD2
        TdD2
        fe
        me
        pD2
        FD2
        be
        Me
        g
        ze
        xi_p
        b
        rhoinf
        vinf
        c
        cl
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
            gd.N1 = obj.N1;
            gd.N2 = obj.N2;
            gd.N3 = obj.N3;
            gd.d  = obj.d;
            gd.h1 = obj.h1;
            gd.h2 = obj.h2;
            geoDiscret = GeometricDiscretizationSolver(gd);
            [obj.x_prim,obj.Tm,obj.Tn] = geoDiscret.compute();
        end
        
        function computeSectionSolver(obj)
            s.x_prim   = obj.x_prim;
            s.Tn       = obj.Tn;
            s.mD1      = obj.mD1;
            s.Tm       = obj.Tm;
            s.open     = obj.open;
            s.E        = obj.E;
            s.G        = obj.G;
            s.ni       = obj.ni;
            s.nel      = obj.nel; 
            s.b        = obj.b;
            s.rhoinf   = obj.rhoinf;
            s.vinf     = obj.vinf;
            s.c        = obj.c;
            s.cl       = obj.cl;
            s.lambda   = obj.lambda;
            s.g        = obj.g;
            s.d        = obj.d;
            s.xi_p     = obj.xi_p;
            s.za       = obj.za;
            s.zm       = obj.zm;
            s.M_x_prim = obj.M_x_prim;
            s.M_y_prim = obj.M_y_prim;
            s.M_z_prim = obj.M_z_prim;
            s.S_x_prim = obj.S_x_prim;
            s.S_y_prim = obj.S_y_prim;
            sec        = SectionSolver(s);
            [obj.sigma,obj.s_norm,obj.tau_s,obj.s_shear,obj.tau_t,obj.s_tor,obj.x_s_prim,obj.I_xx_prim,obj.J] = sec.compute();
        end

        function computemD2(obj)
            obj.mD2 = [obj.E   obj.G   obj.I_xx_prim  obj.J];
        end

        function computeForceMomentElem(obj)
            fm.ni       = obj.ni;
            fm.nel      = obj.nel;
            fm.b        = obj.b;
            fm.rhoinf   = obj.rhoinf;
            fm.vinf     = obj.vinf;
            fm.c        = obj.c;
            fm.cl       = obj.cl;
            fm.lambda   = obj.lambda;
            fm.g        = obj.g;
            fm.d        = obj.d;
            fm.xi_p     = obj.xi_p;
            fm.za       = obj.za;
            fm.zm       = obj.zm;
            fm.x_s_prim = obj.x_s_prim;
            forceMomentElem = ForceMomentElemCompute(fm);
            [obj.xnod,obj.fe,obj.me,obj.ndof] = forceMomentElem.compute();
        end

        function computeBeamSolver(obj)
            bm.nne      = obj.nne;
            bm.ni       = obj.ni;
            bm.nel      = obj.nel;
            bm.xnod     = obj.xnod;
            bm.TnD2     = obj.TnD2;
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
            bm.d        = obj.d;
            bm.xi_p     = obj.xi_p;
            bm.x_s_prim = obj.x_s_prim;
            bm.ze       = obj.ze;
            beam        = BeamSolver(bm);
            [obj.K,obj.F,obj.u,obj.r] = beam.compute();
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.d        = cParams.beamWidth;
            obj.h1       = cParams.h1;
            obj.h2       = cParams.h2;
            obj.N1       = cParams.N1;
            obj.N2       = cParams.N2;
            obj.N3       = cParams.N3;
            obj.mD1      = cParams.materialProp;
            obj.open     = cParams.open;
            obj.M_x_prim = cParams.xBendMoment;
            obj.M_y_prim = cParams.yBendMoment;
            obj.M_z_prim = cParams.zBendMoment;
            obj.S_x_prim = cParams.xShearForce;
            obj.S_y_prim = cParams.yShearForce;
            obj.E        = cParams.E;
            obj.G        = cParams.G;
            obj.ni       = cParams.numDOFsNode;
            obj.nne      = cParams.numNodesElem;
            obj.nel      = cParams.numElements;
            obj.b        = cParams.wingspan;
            obj.rhoinf   = cParams.rhoInf;
            obj.vinf     = cParams.vInf;
            obj.c        = cParams.chord;
            obj.cl       = cParams.Cl;
            obj.lambda   = cParams.lambda;
            obj.TnD2     = cParams.nodalConnec;
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