classdef SectionSolver < handle

    properties (Access = private)
        geomParams
    end

    methods (Access = public)

        function obj = SectionSolver(cParams)
            obj.init(cParams);
        end

        function [xShearCenter,xSecInertia,J] = compute(obj)
            xS    = obj.geomParams.xSection;
            nSecC = obj.geomParams.nodalSecConnec;
            mProp = obj.geomParams.materialProp;
            mC    = obj.geomParams.materialSecConnec;
            open  = obj.geomParams.open;
            [~,~,xShearCenter,~,~,xSecInertia,~,~,J,~] = sectionProperties(xS,nSecC,mProp,mC,open);
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.geomParams = cParams;
        end

    end
end