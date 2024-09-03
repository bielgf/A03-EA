classdef ForceMomentElemCompute < handle
    
    properties (Access = private)
        numDOFsNode
        nel
        b
        rhoinf
        vinf
        c
        cl
        lambda
        g
        d
        xi_p
        za
        zm
        x_s_prim
    end

    methods (Access = public)
        
        function obj = ForceMomentElemCompute(cParams)
            obj.init(cParams);
        end

        function [xnod,fe,me,ndof] = compute(obj)
            [xnod,fe,me] = GetForceMomentElement(obj.nel,obj.b,obj.rhoinf,obj.vinf,obj.c,obj.cl,obj.lambda,obj.g,obj.d,obj.xi_p,obj.za,obj.zm,obj.x_s_prim);
            nnod = size(xnod,2);
            ndof = nnod*obj.numDOFsNode;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.numDOFsNode       = cParams.numDOFsNode;
            obj.nel      = cParams.nel;
            obj.b        = cParams.b;
            obj.rhoinf   = cParams.rhoinf;
            obj.vinf     = cParams.vinf;
            obj.c        = cParams.c;
            obj.cl       = cParams.cl;
            obj.lambda   = cParams.lambda;
            obj.g        = cParams.g;
            obj.d        = cParams.d;
            obj.xi_p     = cParams.xi_p;
            obj.za       = cParams.za;
            obj.zm       = cParams.zm;
            obj.x_s_prim = cParams.x_s_prim;
        end

    end

end