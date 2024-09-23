classdef BoundaryConditionsApplier < handle
    
    properties (Access = private)
        beamParams
        connec
    end

    methods (Access = public)

        function obj = BoundaryConditionsApplier(cParams)
            obj.init(cParams);
        end

        function [up,vp] = computeBC(obj)
            [up,vp] = obj.compute();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.beamParams = cParams.beamP;
            obj.connec     = cParams.con;
        end

        function [up,vp] = compute(obj)
            ni = obj.beamParams.numDOFsNode;
            p  = obj.connec.fixedNodes;

            np = size(p,1);
            vp = zeros(np,1);
            up = zeros(np,1);
            for i = 1:np
                vp(i) = nod2dof(ni,p(i,1),p(i,2));
                up(i) = p(i,3);
            end
        end

    end

end