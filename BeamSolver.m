classdef BeamSolver < handle
    
    properties (Access = private)
        numNodesElem
        numDOFsNode
        numElements
        xnod
        nodalConnec
        dofsConnec
        mD2
        materialConnec
        fixedNodes
        fe
        me
        ndof
        FD2
        xEngine
        engineMass
        g
        beamWidth
        chiP
        x_s_prim
        zEngine
    end

    methods (Access = public)
        
        function obj = BeamSolver(cParams)
            obj.init(cParams);
        end

        function [K,F,u,r] = compute(obj)
            [Kel] = stiffnessFunction(obj.numNodesElem,obj.numDOFsNode,obj.numElements,obj.xnod',obj.nodalConnec,obj.mD2,obj.materialConnec);
            [fel] = forceFunction(obj.numNodesElem,obj.numDOFsNode,obj.numElements,obj.xnod',obj.nodalConnec,obj.fe,obj.me);
            obj.computeFD2();
            [K,F] = assemblyFunction(obj.ndof,obj.numElements,obj.numNodesElem,obj.numDOFsNode,obj.dofsConnec,Kel,fel);
            [up,vp] = applyBC(obj.numDOFsNode,obj.fixedNodes);
            [F] = pointLoads(obj.numDOFsNode,F,obj.FD2);
            [u,r] = solveSystem(obj.ndof,K,F,up,vp);
            [xel,Sel,Mbel,Mtel] = internalforcesFunction(obj.numElements,obj.numDOFsNode,obj.numNodesElem,obj.xnod',obj.nodalConnec,obj.dofsConnec,Kel,u);
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.numNodesElem = cParams.numNodesElem;
            obj.numDOFsNode  = cParams.numDOFsNode;
            obj.numElements = cParams.numElements;
            obj.xnod = cParams.xnod;
            obj.nodalConnec = cParams.nodalConnec;
            obj.dofsConnec = cParams.dofsConnec;
            obj.mD2 = cParams.mD2;
            obj.materialConnec = cParams.materialConnec;
            obj.fixedNodes = cParams.fixedNodes;
            obj.fe = cParams.fe;
            obj.me = cParams.me;
            obj.ndof = cParams.ndof;
            obj.xEngine = cParams.xEngine;
            obj.engineMass = cParams.engineMass;
            obj.g = cParams.g;
            obj.beamWidth = cParams.beamWidth;
            obj.chiP = cParams.chiP;
            obj.x_s_prim = cParams.x_s_prim;
            obj.zEngine = cParams.zEngine;
        end

        function computeFD2(obj)
            obj.FD2 = [find(obj.xnod == obj.xEngine), 1, -obj.engineMass*obj.g;
                       find(obj.xnod == obj.xEngine), 3, -obj.engineMass*obj.g*(((obj.beamWidth + obj.chiP) - obj.x_s_prim) - obj.zEngine)];
        end

    end

end