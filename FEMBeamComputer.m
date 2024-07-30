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
        nnod
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
            obj.sectionSolver();
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

        function sectionSolver(obj)
            % Function to get the section properties
            [obj.x_0_prim,obj.y_0_prim,obj.x_s_prim,obj.y_s_prim,obj.A_tot,obj.I_xx_prim,obj.I_yy_prim,obj.I_xy_prim,obj.J,obj.A_in] = sectionProperties(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,obj.open);
            % Function to get the normal stress distribution
            obj.computemD2();
            obj.computeXnod();
            [obj.sigma,obj.s_norm] = normalStressDistribution(obj.x_prim,obj.Tn,obj.x_0_prim,obj.y_0_prim,obj.I_xx_prim,obj.I_yy_prim,obj.I_xy_prim,obj.M_x_prim,obj.M_y_prim);
            % Function to get the tangential stress distribution due to shear
            [obj.tau_s,obj.s_shear] = tangentialStressDistributionShear(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,obj.x_0_prim,obj.y_0_prim,obj.I_xx_prim,obj.I_yy_prim,obj.I_xy_prim,obj.S_x_prim,obj.S_y_prim,obj.x_s_prim,obj.y_s_prim,obj.A_in,obj.open);
            % Function to get the  tangential stress distribution due to torsion
            [obj.tau_t,obj.s_tor] = tangentialStressDistributionTorsion(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,obj.M_z_prim,obj.J,obj.open,obj.A_in);
        end

        function computemD2(obj)
            obj.mD2 = [ % Young's Modulus, Shear Modulus, Bending Inertia, Torsional Inertia 
       obj.E   obj.G      obj.I_xx_prim       obj.J 
       ];
        end

        function computeXnod(obj)
            [obj.xnod,obj.fe,obj.me] = GetForceMomentElement(obj.data,obj.x_s_prim);
            obj.nnod = size(obj.xnod,2);
            obj.ndof = obj.nnod*obj.ni;

            obj.FD2 = [
        find(obj.xnod == obj.be)       1       -obj.Me*obj.g
        find(obj.xnod == obj.be)       3       -obj.Me*obj.g*(((obj.d + obj.xi_p) - obj.x_s_prim) - obj.ze)
        ];
        end

        function beamsolver(obj)
            [Kel] = stiffnessFunction(obj.data,obj.xnod',obj.TnD2,obj.mD2,obj.TmD2);
            [fel] = forceFunction(obj.data,obj.xnod',obj.TnD2,obj.fe,obj.me);
            obj.data.ndof = obj.ndof;
            [obj.K,obj.F] = assemblyFunction(obj.data,obj.TdD2,Kel,fel);
            [up,vp] = applyBC(obj.data,obj.pD2);
            [obj.F] = pointLoads(obj.data,obj.F,obj.FD2);
            [obj.u,obj.r] = solveSystem(obj.data,obj.K,obj.F,up,vp);
            [xel,Sel,Mbel,Mtel] = internalforcesFunction(obj.data,obj.xnod',obj.TnD2,obj.TdD2,Kel,obj.u);
        end

    end

end