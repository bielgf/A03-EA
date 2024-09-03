classdef BeamSolver < handle
    
    properties (Access = private)
        numNodesElem
        numDOFsNode
        numElements
        xnod
        nodalConnec
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
        beamWidth
        xi_p
        x_s_prim
        ze
    end

    methods (Access = public)
        
        function obj = BeamSolver(cParams)
            obj.init(cParams);
        end

        function [K,F,u,r] = compute(obj)
            [Kel] = stiffnessFunction(obj.numNodesElem,obj.numDOFsNode,obj.numElements,obj.xnod',obj.nodalConnec,obj.mD2,obj.TmD2);
            [fel] = forceFunction(obj.numNodesElem,obj.numDOFsNode,obj.numElements,obj.xnod',obj.nodalConnec,obj.fe,obj.me);
            obj.computeFD2();
            [K,F] = assemblyFunction(obj.ndof,obj.numElements,obj.numNodesElem,obj.numDOFsNode,obj.TdD2,Kel,fel);
            [up,vp] = applyBC(obj.numDOFsNode,obj.pD2);
            [F] = pointLoads(obj.numDOFsNode,F,obj.FD2);
            [u,r] = solveSystem(obj.ndof,K,F,up,vp);
            [xel,Sel,Mbel,Mtel] = internalforcesFunction(obj.numElements,obj.numDOFsNode,obj.numNodesElem,obj.xnod',obj.nodalConnec,obj.TdD2,Kel,u);
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.numNodesElem = cParams.numNodesElem;
            obj.numDOFsNode  = cParams.numDOFsNode;
            obj.numElements = cParams.numElements;
            obj.xnod = cParams.xnod;
            obj.nodalConnec = cParams.nodalConnec;
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
            obj.beamWidth = cParams.beamWidth;
            obj.xi_p = cParams.xi_p;
            obj.x_s_prim = cParams.x_s_prim;
            obj.ze = cParams.ze;
        end

        function computeFD2(obj)
            obj.FD2 = [find(obj.xnod == obj.be), 1, -obj.Me*obj.g; find(obj.xnod == obj.be), 3, -obj.Me*obj.g*(((obj.beamWidth + obj.xi_p) - obj.x_s_prim) - obj.ze)];
        end

    end

end