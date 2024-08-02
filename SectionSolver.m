classdef SectionSolver < handle

    properties (Access = private)
        x_prim % xDomain, xSection, xGlobal, xLocal.. ?
        Tn % connec
        mD1
        Tm
        open
        E
        G
        ni
        nel
        b
        rhoinf
        vinf
        c
        cl
        lambda
        g
        d
        xi_p
        za
        zm
        M_x_prim
        M_y_prim
        M_z_prim
        S_x_prim
        S_y_prim
    end

    properties (Access = public)
        mD2
        xnod
        fe
        me
        ndof
    end

    methods (Access = public)
        function obj = SectionSolver(cParams)
            obj.init(cParams);
        end

        function [sigma,s_norm,tau_s,s_shear,tau_t,s_tor,x_s_prim] = compute(obj)
            % Function to get the section properties
            [x_0_prim,y_0_prim,x_s_prim,y_s_prim,~,I_xx_prim,I_yy_prim,I_xy_prim,J,A_in] = sectionProperties(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,obj.open);
            % Function to get the normal stress distribution
            obj.computemD2(I_xx_prim,J);
            obj.computeXnod(x_s_prim);
            [sigma,s_norm] = normalStressDistribution(obj.x_prim,obj.Tn,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,obj.M_x_prim,obj.M_y_prim);
            % Function to get the tangential stress distribution due to shear
            [tau_s,s_shear] = tangentialStressDistributionShear(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,obj.S_x_prim,obj.S_y_prim,x_s_prim,y_s_prim,A_in,obj.open);
            % Function to get the  tangential stress distribution due to torsion
            [tau_t,s_tor] = tangentialStressDistributionTorsion(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,obj.M_z_prim,J,obj.open,A_in);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.x_prim = cParams.x_prim;
            obj.Tn     = cParams.Tn;
            obj.mD1    = cParams.mD1;
            obj.Tm     = cParams.Tm;
            obj.open   = cParams.open;
            obj.E      = cParams.E;
            obj.G      = cParams.G;
            obj.ni     = cParams.ni;
            obj.M_x_prim = cParams.M_x_prim;
            obj.M_y_prim = cParams.M_y_prim;
            obj.M_z_prim = cParams.M_z_prim;
            obj.S_x_prim = cParams.S_x_prim;
            obj.S_y_prim = cParams.S_y_prim;
            obj.nel      = cParams.nel;
            obj.b        = cParams.b;
            obj.rhoinf   = cParams.rhoinf;
            obj.vinf     = cParams.vinf;
            obj.c        = cParams.c;
            obj.cl       = cParams.cl;
            obj.lambda   = cParams.lambda;
            obj.g        = cParams.g;
            obj.d        = cParams.d;
            obj.xi_p     = cParams.xi_p;
            obj.za       = cParams.za;
            obj.zm       = cParams.zm;
        end

        function computemD2(obj,I_xx_prim,J)
            obj.mD2 = [ % Young's Modulus, Shear Modulus, Bending Inertia, Torsional Inertia
                obj.E   obj.G      I_xx_prim       J
                ];
        end

        function computeXnod(obj,x_s_prim)
            [obj.xnod,obj.fe,obj.me] = GetForceMomentElement(obj.nel,obj.b,obj.rhoinf,obj.vinf,obj.c,obj.cl,obj.lambda,obj.g,obj.d,obj.xi_p,obj.za,obj.zm,x_s_prim);
            nnod = size(obj.xnod,2);
            obj.ndof = nnod*obj.ni;
        end
    end
end