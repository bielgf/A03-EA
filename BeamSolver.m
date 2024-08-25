classdef BeamSolver < handle
    
    properties (Access = private)
        nne
        ni
        nel
        xnod
        TnD2
        TdD2
        mD2
        TmD2
        pD2
        fe
        me
        ndof
        FD2
        be
        Me
        g
        d
        xi_p
        x_s_prim
        ze
    end

    methods (Access = public)
        
        function obj = BeamSolver(cParams)
            obj.init(cParams);
        end

        function [K,F,u,r] = compute(obj)
            [Kel] = stiffnessFunction(obj.nne,obj.ni,obj.nel,obj.xnod',obj.TnD2,obj.mD2,obj.TmD2);
            [fel] = forceFunction(obj.nne,obj.ni,obj.nel,obj.xnod',obj.TnD2,obj.fe,obj.me);
            obj.computeFD2();
            [K,F] = assemblyFunction(obj.ndof,obj.nel,obj.nne,obj.ni,obj.TdD2,Kel,fel);
            [up,vp] = applyBC(obj.ni,obj.pD2);
            [F] = pointLoads(obj.ni,F,obj.FD2);
            [u,r] = solveSystem(obj.ndof,K,F,up,vp);
            [xel,Sel,Mbel,Mtel] = internalforcesFunction(obj.nel,obj.ni,obj.nne,obj.xnod',obj.TnD2,obj.TdD2,Kel,u);
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.nne = cParams.nne;
            obj.ni  = cParams.ni;
            obj.nel = cParams.nel;
            obj.xnod = cParams.xnod;
            obj.TnD2 = cParams.TnD2;
            obj.TdD2 = cParams.TdD2;
            obj.mD2 = cParams.mD2;
            obj.TmD2 = cParams.TmD2;
            obj.pD2 = cParams.pD2;
            obj.fe = cParams.fe;
            obj.me = cParams.me;
            obj.ndof = cParams.ndof;
            obj.be = cParams.be;
            obj.Me = cParams.Me;
            obj.g = cParams.g;
            obj.d = cParams.d;
            obj.xi_p = cParams.xi_p;
            obj.x_s_prim = cParams.x_s_prim;
            obj.ze = cParams.ze;
        end

        function computeFD2(obj)
            obj.FD2 = [find(obj.xnod == obj.be), 1, -obj.Me*obj.g; find(obj.xnod == obj.be), 3, -obj.Me*obj.g*(((obj.d + obj.xi_p) - obj.x_s_prim) - obj.ze)];
        end

    end

end