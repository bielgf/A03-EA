classdef BeamSolver < handle
    
    properties (Access = private)
        beamParams
        connec
        beamProp
        externalForce

        numNodesElem
        numDOFsNode
        numElements
        xGlobal
        nodalConnec
        dofsConnec
        forceElem
        momentElem
        totalDOFs
        fixedNodes
    end

    methods (Access = public)
        
        function obj = BeamSolver(cParams)
            obj.init(cParams);
        end

        function [K,F,u,r] = compute(obj)
            s.beamP  = obj.beamParams;
            s.con    = obj.connec;
            s.beamPr = obj.beamProp;
            s.extF   = obj.externalForce;

            sF = StiffnessFunctionClass(s);
            sF.computeStiffnessMatrix();
            s.Kel = sF.Kel;

            fF = ForceFunctionClass(s);
            fF.computeForceVector();
            s.fel = fF.fel;

            aF = AssemblyFunctionClass(s);
            aF.computeAssembly();
            K = aF.K;
            F = aF.F;

            bc = applyBoundaryConditionsClass(s);
            bc.computeBC();
            up = bc.up;
            vp = bc.vp;




%             [Kel] = stiffnessFunction(obj.numNodesElem,obj.numDOFsNode,obj.numElements,obj.xGlobal',obj.nodalConnec,obj.beamProp,obj.materialConnec); % passar a class
%             [fel] = forceFunction(obj.numNodesElem,obj.numDOFsNode,obj.numElements,obj.xGlobal',obj.nodalConnec,obj.forceElem,obj.momentElem); % passar a class
%             [K,F] = assemblyFunction(obj.totalDOFs,obj.numElements,obj.numNodesElem,obj.numDOFsNode,obj.dofsConnec,Kel,fel); % passar a class
%             [up,vp] = applyBC(obj.numDOFsNode,obj.fixedNodes); % passar a class
            [F] = pointLoads(obj.numDOFsNode,F,obj.externalForce); % passar a class
            [u,r] = solveSystem(obj.totalDOFs,K,F,up,vp);
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.beamParams    = cParams.beamParams;
            obj.connec        = cParams.connec;
            obj.beamProp      = cParams.beamProp;
            obj.externalForce = cParams.externalForce;

            obj.numNodesElem = cParams.numNodesElem;
            obj.numDOFsNode = cParams.numDOFsNode;
            obj.numElements = cParams.numElements;
            obj.xGlobal = cParams.xGlobal;
            obj.nodalConnec = cParams.nodalConnec;
            obj.dofsConnec = cParams.dofsConnec;
            obj.forceElem = cParams.forceElem;
            obj.momentElem = cParams.momentElem;
            obj.totalDOFs = cParams.totalDOFs;
            obj.fixedNodes = cParams.fixedNodes;
        end

    end

end