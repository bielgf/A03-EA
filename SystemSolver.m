classdef SystemSolver < handle
    
    properties (Access = private)
        beamParams % -> nDofTotal
        K
        F
        up
        vp
    end

    properties (Access = public)
        u
        r
    end

    methods (Access = public)

        function obj = SystemSolver(cParams)
            obj.init(cParams);
        end

        function computeSystem(obj)
            obj.compute();
        end

    end

    methods (Access = public)

        function init(obj,cParams)
            obj.beamParams = cParams.beamP;
            obj.K          = cParams.K;
            obj.F          = cParams.F;
            obj.up         = cParams.up;
            obj.vp         = cParams.vp;
        end

        function compute(obj)
            ndof = obj.beamParams.totalDOFs;
            k    = obj.K;
            f    = obj.F;
            upoint   = obj.up;
            vpoint   = obj.vp;

            vf = setdiff((1:ndof)',vpoint);
            obj.u = zeros(ndof,1);
            obj.u(vpoint) = upoint;
            LHS = k(vf,vf);
            RHS = f(vf) - k(vf,vpoint)*obj.u(vpoint); % private function

            method = 'Direct';      % method = Direct or Iterative
            solver = Solver.create(LHS,RHS,method);
            solver.compute();
            obj.u(vf) = solver.x;
            obj.r = k(vpoint,:)*obj.u-f(vpoint); % private function
        end

    end

end