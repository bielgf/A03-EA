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
        ndof
        data
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
            obj.geometricdiscretization();
        end
        
        function computeSectionSolver(obj)
            s.x_prim = obj.x_prim;
            s.Tn     = obj.Tn;
            s.mD1    = obj.mD1;
            s.Tm     = obj.Tm;
            s.open   = obj.open;
            s.E      = obj.E;
            s.G      = obj.G;
            s.ni     = obj.ni;
            s.M_x_prim = obj.M_x_prim;
            s.M_y_prim = obj.M_y_prim;
            s.M_z_prim = obj.M_z_prim;
            s.S_x_prim = obj.S_x_prim;
            s.S_y_prim = obj.S_y_prim;
            s.data   = obj.data; % delete..
            sec      = SectionSolver(s);
            [obj.sigma,obj.s_norm,obj.tau_s,obj.s_shear,obj.tau_t,obj.s_tor,obj.x_s_prim] = sec.compute();
            obj.xnod = sec.xnod;
            obj.mD2  = sec.mD2;
            obj.fe   = sec.fe;
            obj.me   = sec.me;
            obj.ndof = sec.ndof;
        end

        function computeBeamSolver(obj)
            obj.beamsolver();
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.d        = cParams.d;
            obj.h1       = cParams.h1;
            obj.h2       = cParams.h2;
            obj.N1       = cParams.N1;
            obj.N2       = cParams.N2;
            obj.N3       = cParams.N3;
            obj.mD1      = cParams.mD1;
            obj.open     = cParams.open;
            obj.M_x_prim = cParams.M_x_prim;
            obj.M_y_prim = cParams.M_y_prim;
            obj.M_z_prim = cParams.M_z_prim;
            obj.S_x_prim = cParams.S_x_prim;
            obj.S_y_prim = cParams.S_y_prim;
            obj.E        = cParams.E;
            obj.G        = cParams.G;
            obj.ni       = cParams.ni;
            obj.data     = cParams; % ELIMINAR-HO !!
            obj.TnD2     = cParams.TnD2;
            obj.TmD2     = cParams.TmD2;
            obj.TdD2     = cParams.TdD2;
            obj.pD2      = cParams.pD2;
            obj.be       = cParams.be;
            obj.Me       = cParams.Me;
            obj.g        = cParams.g;
            obj.ze       = cParams.ze;
            obj.xi_p     = cParams.xi_p;
        end

        function geometricdiscretization(obj)
            [obj.x_prim,obj.Tm,obj.Tn] = GeometricDiscret(obj.N1,obj.N2,obj.N3,obj.d,obj.h1,obj.h2);
        end

        function beamsolver(obj)
            [Kel] = stiffnessFunction(obj.data,obj.xnod',obj.TnD2,obj.mD2,obj.TmD2);
            [fel] = forceFunction(obj.data,obj.xnod',obj.TnD2,obj.fe,obj.me);
            
            obj.FD2 = [find(obj.xnod == obj.be), 1, -obj.Me*obj.g; find(obj.xnod == obj.be), 3, -obj.Me*obj.g*(((obj.d + obj.xi_p) - obj.x_s_prim) - obj.ze)];
            obj.data.ndof = obj.ndof;
            
            [obj.K,obj.F] = assemblyFunction(obj.data,obj.TdD2,Kel,fel);
            [up,vp] = applyBC(obj.data,obj.pD2);
            [obj.F] = pointLoads(obj.data,obj.F,obj.FD2);
            [obj.u,obj.r] = solveSystem(obj.data,obj.K,obj.F,up,vp);
            [xel,Sel,Mbel,Mtel] = internalforcesFunction(obj.data,obj.xnod',obj.TnD2,obj.TdD2,Kel,obj.u);
        end

    end

end