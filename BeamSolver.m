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
        externalForce
        totalDOFs
    end

    methods (Access = public)
        
        function obj = BeamSolver(cParams)
            obj.init(cParams);
        end

        function [K,F,u,r] = compute(obj)
            [Kel] = stiffnessFunction(obj.numNodesElem,obj.numDOFsNode,obj.numElements,obj.xGlobal',obj.nodalConnec,obj.beamProp,obj.materialConnec); % passar a class
            [fel] = forceFunction(obj.numNodesElem,obj.numDOFsNode,obj.numElements,obj.xGlobal',obj.nodalConnec,obj.forceElem,obj.momentElem); % passar a class
            [K,F] = assemblyFunction(obj.totalDOFs,obj.numElements,obj.numNodesElem,obj.numDOFsNode,obj.dofsConnec,Kel,fel); % passar a class
            [up,vp] = applyBC(obj.numDOFsNode,obj.fixedNodes); % passar a class
            [F] = pointLoads(obj.numDOFsNode,F,obj.externalForce); % passar a class
            [u,r] = solveSystem(obj.totalDOFs,K,F,up,vp);
            [xel,Sel,Mbel,Mtel] = internalforcesFunction(obj.numElements,obj.numDOFsNode,obj.numNodesElem,obj.xGlobal',obj.nodalConnec,obj.dofsConnec,Kel,u); % passar a class
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
            obj.externalForce = cParams.externalForce;
            obj.totalDOFs = cParams.totalDOFs;
        end

    end

end