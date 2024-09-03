classdef BeamSolver < handle
    
    properties (Access = private)
        numNodesElem
        numDOFsNode
        numElements
        xGlobal
        nodalConnec
        dofsConnec
        beamProp
        materialConnec
        fixedNodes
        forceElem
        momentElem
        totalDOFs
        FD2
        xEngine
        engineMass
        g
        beamWidth
        chiP
        xShearCenter
        zEngine
    end

    methods (Access = public)
        
        function obj = BeamSolver(cParams)
            obj.init(cParams);
        end

        function [K,F,u,r] = compute(obj)
            [Kel] = stiffnessFunction(obj.numNodesElem,obj.numDOFsNode,obj.numElements,obj.xGlobal',obj.nodalConnec,obj.beamProp,obj.materialConnec);
            [fel] = forceFunction(obj.numNodesElem,obj.numDOFsNode,obj.numElements,obj.xGlobal',obj.nodalConnec,obj.forceElem,obj.momentElem);
            obj.computeFD2();
            [K,F] = assemblyFunction(obj.totalDOFs,obj.numElements,obj.numNodesElem,obj.numDOFsNode,obj.dofsConnec,Kel,fel);
            [up,vp] = applyBC(obj.numDOFsNode,obj.fixedNodes);
            [F] = pointLoads(obj.numDOFsNode,F,obj.FD2);
            [u,r] = solveSystem(obj.totalDOFs,K,F,up,vp);
            [xel,Sel,Mbel,Mtel] = internalforcesFunction(obj.numElements,obj.numDOFsNode,obj.numNodesElem,obj.xGlobal',obj.nodalConnec,obj.dofsConnec,Kel,u);
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.numNodesElem = cParams.numNodesElem;
            obj.numDOFsNode  = cParams.numDOFsNode;
            obj.numElements = cParams.numElements;
            obj.xGlobal = cParams.xGlobal;
            obj.nodalConnec = cParams.nodalConnec;
            obj.dofsConnec = cParams.dofsConnec;
            obj.beamProp = cParams.beamProp;
            obj.materialConnec = cParams.materialConnec;
            obj.fixedNodes = cParams.fixedNodes;
            obj.forceElem = cParams.forceElem;
            obj.momentElem = cParams.momentElem;
            obj.totalDOFs = cParams.totalDOFs;
            obj.xEngine = cParams.xEngine;
            obj.engineMass = cParams.engineMass;
            obj.g = cParams.g;
            obj.beamWidth = cParams.beamWidth;
            obj.chiP = cParams.chiP;
            obj.xShearCenter = cParams.xShearCenter;
            obj.zEngine = cParams.zEngine;
        end

        function computeFD2(obj)
            obj.FD2 = [find(obj.xGlobal == obj.xEngine), 1, -obj.engineMass*obj.g;
                       find(obj.xGlobal == obj.xEngine), 3, -obj.engineMass*obj.g*(((obj.beamWidth + obj.chiP) - obj.xShearCenter) - obj.zEngine)];
        end

    end

end