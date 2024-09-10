classdef ForceMomentElemCompute < handle
    
    properties (Access = private)
        numDOFsNode
        numElements
        wingspan
        rhoInf
        vInf
        chord
        Cl
        lambda
        g
        beamWidth
        chiP
        aeroCenter
        centerOfMass
        xShearCenter
    end

    methods (Access = public)
        
        function obj = ForceMomentElemCompute(cParams)
            obj.init(cParams);
        end

        function [xGlobal,forceElem,momentElem,totalDOFs] = compute(obj)
            %nElem = obj.numElements;
            [xGlobal,forceElem,momentElem] = GetForceMomentElement(obj.numElements,obj.wingspan,obj.rhoInf,obj.vInf,obj.chord,obj.Cl,obj.lambda,obj.g,obj.beamWidth,obj.chiP,obj.aeroCenter,obj.centerOfMass,obj.xShearCenter);
            nnod = size(xGlobal,2);
            totalDOFs = nnod*obj.numDOFsNode;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.numDOFsNode  = cParams.numDOFsNode;
            obj.numElements  = cParams.numElements;
            obj.wingspan     = cParams.wingspan;
            obj.rhoInf       = cParams.rhoInf;
            obj.vInf         = cParams.vInf;
            obj.chord        = cParams.chord;
            obj.Cl           = cParams.Cl;
            obj.lambda       = cParams.lambda;
            obj.g            = cParams.g;
            obj.beamWidth    = cParams.beamWidth;
            obj.chiP         = cParams.chiP;
            obj.aeroCenter   = cParams.aeroCenter;
            obj.centerOfMass = cParams.centerOfMass;
            obj.xShearCenter = cParams.xShearCenter;
        end

    end

end