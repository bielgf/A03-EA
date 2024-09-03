classdef SectionSolver < handle

    properties (Access = private)
        x_prim  % xDomain, xSection, xGlobal, xLocal.. ?
        Tn
        materialProp
        Tm
        open
        xBendMoment
        yBendMoment
        zBendMoment
        xShearForce
        yShearForce
    end

    methods (Access = public)
        function obj = SectionSolver(cParams)
            obj.init(cParams);
        end

        function [sigma,s_norm,tau_s,s_shear,tau_t,s_tor,x_s_prim,I_xx_prim,J] = compute(obj)
            [x_0_prim,y_0_prim,x_s_prim,y_s_prim,~,I_xx_prim,I_yy_prim,I_xy_prim,J,A_in] = sectionProperties(obj.x_prim,obj.Tn,obj.materialProp,obj.Tm,obj.open);
            [sigma,s_norm] = normalStressDistribution(obj.x_prim,obj.Tn,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,obj.xBendMoment,obj.yBendMoment);
            [tau_s,s_shear] = tangentialStressDistributionShear(obj.x_prim,obj.Tn,obj.materialProp,obj.Tm,x_0_prim,y_0_prim,I_xx_prim,I_yy_prim,I_xy_prim,obj.xShearForce,obj.yShearForce,x_s_prim,y_s_prim,A_in,obj.open);
            [tau_t,s_tor] = tangentialStressDistributionTorsion(obj.x_prim,obj.Tn,obj.materialProp,obj.Tm,obj.zBendMoment,J,obj.open,A_in);
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.x_prim   = cParams.x_prim;
            obj.Tn       = cParams.Tn;
            obj.materialProp      = cParams.materialProp;
            obj.Tm       = cParams.Tm;
            obj.open     = cParams.open;
            obj.xBendMoment = cParams.xBendMoment;
            obj.yBendMoment = cParams.yBendMoment;
            obj.zBendMoment = cParams.zBendMoment;
            obj.xShearForce = cParams.xShearForce;
            obj.yShearForce = cParams.yShearForce;
        end

    end
end