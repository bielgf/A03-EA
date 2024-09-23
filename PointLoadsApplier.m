classdef PointLoadsApplier < handle
    
    properties (Access = private)
        nDOFsNode
        externalForce
    end

    properties
        F
    end

    methods (Access = public)
        
        function obj = PointLoadsApplier(cParams)
            obj.init(cParams);
        end

        function computePointLoads(obj)
            obj.compute;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nDOFsNode     = cParams.nDofsNode;
            obj.externalForce = cParams.extF;
            obj.F             = cParams.F;
        end

        function compute(obj)
            ni = obj.nDOFsNode;
            Fc  = obj.externalForce;

            for i = 1:size(Fc,1)
                I = nod2dof(ni,Fc(i,1),Fc(i,2));
                obj.F(I) = Fc(i,3);
            end

        end

    end
    
end