classdef SectionSolver < handle

    properties (Access = private)
        x_prim  % xDomain, xSection, xGlobal, xLocal.. ?
        Tn      % connec
        mD1
        Tm
        open
        M_x_prim
        M_y_prim
        M_z_prim
        S_x_prim
        S_y_prim
    end

    methods (Access = public)
        function obj = SectionSolver(cParams)
            obj.init(cParams);
        end

        function [sigma,s_norm,tau_s,s_shear,tau_t,s_tor,x_s_prim,I_xx_prim,J] = compute(obj)
            [x_0_prim,y_0_prim,x_s_prim,y_s_prim,~,I_xx_prim,I_yy_prim,I_xy_prim,J,A_in] = sectionProperties(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,obj.open);
            [sigma,s_norm] = normalStressDistribution(obj.x_prim,obj.Tn,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,obj.M_x_prim,obj.M_y_prim);
            [tau_s,s_shear] = tangentialStressDistributionShear(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,obj.S_x_prim,obj.S_y_prim,x_s_prim,y_s_prim,A_in,obj.open);
            [tau_t,s_tor] = tangentialStressDistributionTorsion(obj.x_prim,obj.Tn,obj.mD1,obj.Tm,obj.M_z_prim,J,obj.open,A_in);
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.x_prim   = cParams.x_prim;
            obj.Tn       = cParams.Tn;
            obj.mD1      = cParams.mD1;
            obj.Tm       = cParams.Tm;
            obj.open     = cParams.open;
            obj.M_x_prim = cParams.M_x_prim;
            obj.M_y_prim = cParams.M_y_prim;
            obj.M_z_prim = cParams.M_z_prim;
            obj.S_x_prim = cParams.S_x_prim;
            obj.S_y_prim = cParams.S_y_prim;
        end

    end
end