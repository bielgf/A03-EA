classdef MAIN < handle
    
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
    end

    properties (Access = public)
        sigma
        s_norm
        tau_s
        s_shear
        tau_t
        s_tor
    end

    methods (Access = public)
        
        function obj = MAIN(cParams)
            obj.init(cParams);
        end

        function computeGeoDiscret(obj)
            obj.geometricdiscretization();
        end
        
        function computeSectionSolver(obj)
            obj.sectionsolver();
        end

        function computeBeamSolver(obj)
            obj.beamsolver();
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.d           = cParams.d;
            obj.h1          = cParams.h1;
            obj.h2          = cParams.h2;
            obj.N1          = cParams.N1;
            obj.N2          = cParams.N2;
            obj.N3          = cParams.N3;
            obj.mD1         = cParams.mD1;
            obj.open        = cParams.open;
            obj.M_x_prim    = cParams.M_x_prim;
            obj.M_y_prim    = cParams.M_y_prim;
            obj.M_z_prim    = cParams.M_z_prim;
            obj.S_x_prim    = cParams.S_x_prim;
            obj.S_y_prim    = cParams.S_y_prim;
        end

        function geometricdiscretization(obj)
            [obj.x_prim,obj.Tm,obj.Tn] = GeometricDiscret(obj.N1,obj.N2,obj.N3,obj.d,obj.h1,obj.h2);
        end

        function sectionsolver(obj)
            % Function to get the section properties
            [obj.x_0_prim,obj.y_0_prim,obj.x_s_prim,obj.y_s_prim,obj.A_tot,obj.I_xx_prim,obj.I_yy_prim,obj.I_xy_prim,obj.J,obj.A_in] = sectionProperties(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,obj.open);
            % Function to get the normal stress distribution
            [obj.sigma,obj.s_norm] = normalStressDistribution(obj.x_prim,obj.Tn,obj.x_0_prim,obj.y_0_prim,obj.I_xx_prim,obj.I_yy_prim,obj.I_xy_prim,obj.M_x_prim,obj.M_y_prim);
            % Function to get the tangential stress distribution due to shear
            [obj.tau_s,obj.s_shear] = tangentialStressDistributionShear(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,obj.x_0_prim,obj.y_0_prim,obj.I_xx_prim,obj.I_yy_prim,obj.I_xy_prim,obj.S_x_prim,obj.S_y_prim,obj.x_s_prim,obj.y_s_prim,obj.A_in,obj.open);
            % Function to get the  tangential stress distribution due to torsion
            [obj.tau_t,obj.s_tor] = tangentialStressDistributionTorsion(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,obj.M_z_prim,obj.J,obj.open,obj.A_in);
        end

        function beamsolver(obj)
            [Kel] = stiffnessFunction(data,xnod',TnD2,mD2,TmD2);
            [fel] = forceFunction(data,xnod',TnD2,fe,me);
            [K,F] = assemblyFunction(data,TdD2,Kel,fel);
            [up,vp] = applyBC(data,pD2);
            [F] = pointLoads(data,F,FD2);
            [u,r] = solveSystem(data,K,F,up,vp);
            [xel,Sel,Mbel,Mtel] = internalforcesFunction(data,xnod',TnD2,TdD2,Kel,u);
        end

    end

end