classdef PointLoadsApplier < handle
    
    properties (Access = private)
        beamParams
        externalForce
    end

    properties (Access = public)
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
            obj.beamParams    = cParams.beamP;
            obj.externalForce = cParams.extF;
            obj.F             = cParams.F;
        end

        function compute(obj)
            ni = obj.beamParams.numDOFsNode;
            Fc  = obj.externalForce;

            for i = 1:size(Fc,1)
                I = nod2dof(ni,Fc(i,1),Fc(i,2));
                obj.F(I) = Fc(i,3);
            end

        end

    end
    
end