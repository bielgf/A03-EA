classdef SectionSolver < handle

    properties (Access = private)
        xSection
        nodalSecConnec
        materialProp
        materialSecConnec
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

        function [sigma,s_norm,tau_s,s_shear,tau_t,s_tor,xShearCenter,xSecInertia,J] = compute(obj)
            % xS = obj.xSection;
            % ...
            [x_0_prim,y_0_prim,xShearCenter,y_s_prim,~,xSecInertia,I_yy_prim,I_xy_prim,J,A_in] = sectionProperties(obj.xSection,obj.nodalSecConnec,obj.materialProp,obj.materialSecConnec,obj.open);
            %secProps = sectionProperties(xS,...)
            [sigma,s_norm] = normalStressDistribution(obj.xSection,obj.nodalSecConnec,x_0_prim,y_0_prim,xSecInertia,I_yy_prim,I_xy_prim,obj.xBendMoment,obj.yBendMoment);
            [tau_s,s_shear] = tangentialStressDistributionShear(obj.xSection,obj.nodalSecConnec,obj.materialProp,obj.materialSecConnec,x_0_prim,y_0_prim,xSecInertia,I_yy_prim,I_xy_prim,obj.xShearForce,obj.yShearForce,xShearCenter,y_s_prim,A_in,obj.open);
            [tau_t,s_tor] = tangentialStressDistributionTorsion(obj.xSection,obj.nodalSecConnec,obj.materialProp,obj.materialSecConnec,obj.zBendMoment,J,obj.open,A_in);
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.xSection          = cParams.xSection;
            obj.nodalSecConnec    = cParams.nodalSecConnec;
            obj.materialProp      = cParams.materialProp;
            obj.materialSecConnec = cParams.materialSecConnec;
            obj.open              = cParams.open;
            obj.xBendMoment       = cParams.xBendMoment;
            obj.yBendMoment       = cParams.yBendMoment;
            obj.zBendMoment       = cParams.zBendMoment;
            obj.xShearForce       = cParams.xShearForce;
            obj.yShearForce       = cParams.yShearForce;
        end

    end
end