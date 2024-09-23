classdef BoundaryConditionsApplier < handle
    
    properties (Access = private)
        beamParams
        connec
    end

    properties (Access = public)
        up
        vp
    end

    methods (Access = public)

        function obj = BoundaryConditionsApplier(cParams)
            obj.init(cParams);
        end

        function computeBC(obj)
            obj.compute();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.beamParams = cParams.beamP;
            obj.connec     = cParams.con;
        end

        function compute(obj)
            ni = obj.beamParams.numDOFsNode;
            p  = obj.connec.fixedNodes;

            np = size(p,1);
            obj.vp = zeros(np,1);
            obj.up = zeros(np,1);
            for i = 1:np
                obj.vp(i) = nod2dof(ni,p(i,1),p(i,2));
                obj.up(i) = p(i,3);
            end
        end

    end

end