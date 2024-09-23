classdef SystemSolver < handle
    
    properties (Access = private)
        nDOFTotal
        K
        F
        up
        vp
    end

    methods (Access = public)

        function obj = SystemSolver(cParams)
            obj.init(cParams);
        end

        function [u,r] = computeSystem(obj)
            [u,r] = obj.compute();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nDOFTotal  = cParams.nDofTotal;
            obj.K          = cParams.K;
            obj.F          = cParams.F;
            obj.up         = cParams.up;
            obj.vp         = cParams.vp;
        end

        function [u,r] = compute(obj)
            ndof = obj.nDOFTotal;
            k    = obj.K;
            f    = obj.F;
            upoint   = obj.up;
            vpoint   = obj.vp;

            vf = setdiff((1:ndof)',vpoint);
            u = zeros(ndof,1);
            u(vpoint) = upoint;
            LHS = k(vf,vf);
            RHS = f(vf) - k(vf,vpoint)*u(vpoint);

            method = 'Direct';                      % method = Direct or Iterative
            solver = Solver.create(LHS,RHS,method);
            solver.compute();
            u(vf) = solver.x;
            r = k(vpoint,:)*u-f(vpoint);
        end

    end

end